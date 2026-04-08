# ============================================================================
#  w10c_eof_pca_circulation_indices.R
#  Nechako Watershed Drought — EOF/PCA CIRCULATION INDICES
#
#  PURPOSE
#  -------
#  Replace and augment the simple composite description in Section 4.7 of the
#  manuscript with quantitative PC time series as circulation indices.
#  Implements the full framework :
#    Phase 1 — Prepare North Pacific anomaly fields (cosine-weighted, standardised)
#    Phase 2 — EOF/PCA decomposition of Z500, SLP, SST (all-month + JJA + DJF)
#    Phase 3 — PC–SPEI correlation, variance partitioning, lag analysis
#    Phase 4 — 2022–2025 event quantification (ranking, z-scores, GEV)
#    Phase 5 — Figures & tables for Sections 4.7.1, 4.7.2 of the manuscript
#
#  OUTPUTS  (all written to OUT_DIR/eof_analysis/)
#  -------
#    w10c_eof_state.rds        — full EOF results (loadings, PCs, variances)
#    eof_pc_timeseries.csv     — monthly PC time series for all fields
#    Table4_PC_correlations.csv — Table 4: PC vs SPEI + teleconnection indices
#    Table5_regression.csv     — Table 5: SPEI ~ PC1 + PC2 + PC3 results
#    ms2_eof_blanks.csv        — manuscript fill-in values (§4.7.1, §4.7.2)
#    FigS1_scree_diagrams.pdf/png
#    FigS2_eof_spatial_patterns.pdf/png
#    Fig8_pc_timeseries_SPEI.pdf/png
#    Fig9_correlation_heatmap.pdf/png
#
#  SCRIPT CHAIN
#  -----------
#  Run AFTER : w9_atmospheric_diagnostics.R (anomaly NetCDFs + shared state)
#              w1–w3 (SPEI basin CSVs must exist)
#              w10b_further_atm_diag_index_skill.R (optional — for PDO/PNA cache)
#  Run BEFORE: (none — terminal analysis step; feeds manuscript directly)
#
#  PACKAGE REQUIREMENTS
#  -------------------
#  Required  : terra, ggplot2, patchwork, scales, dplyr, tidyr, lubridate,
#               readr, sf, rnaturalearth, rnaturalearthdata, lmtest, sandwich
#  Optional  : evd (GEV fitting; empirical fallback used if absent)
#
#  REFERENCES
#  ----------
#  Bretherton et al. (1999) effective sample size formula.
#  Jolliffe & Cadima (2016) for EOF/PCA methodology.
#  Newey & West (1987) for autocorrelation-robust regression SEs.
#
#  CHANGELOG
#  ---------
#  Section 4: added diagnostic warning when eof_results is empty so the
#        root cause (missing w9 NetCDFs or failed prcomp) surfaces immediately
#        rather than propagating silently to Section 11.
#  Section 11: replaced  description = NA_character_  with
#        description = rep(NA_character_, length(ms_blanks))  so that
#        data.frame() does not crash with "arguments imply differing number
#        of rows: 0, 1" when ms_blanks is empty (all EOF combinations failed).
#        Also added an explicit empty-guard that writes a zero-row CSV and
#        prints a clear warning instead of aborting the script.
# ============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(readr)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(lmtest)
  library(sandwich)
})

# ============================================================================
#  SECTION 1: CONFIGURATION
#  Loads scalar constants from w9 shared state; adds EOF-specific constants.
# ============================================================================
cat("============================================================\n")
cat(" SECTION 1: Configuration\n")
cat("============================================================\n")

WD_PATH  <- "D:/Nechako_Drought/Nechako/"
DATA_DIR <- "D:/Nechako_Drought/Nechako/monthly_data_direct"
ATM_DIR  <- file.path(WD_PATH, "atmospheric_diagnostics")    # w9 output dir

# ── Load w9 shared state (inherits all configuration scalars) ─────────────────
rds_w9 <- file.path(ATM_DIR, "w9_shared_state.rds")
if (!file.exists(rds_w9))
  stop("w9 shared state not found: ", rds_w9,
       "\n  Run w9_atmospheric_diagnostics.R first.")
w9 <- readRDS(rds_w9)
list2env(w9, envir = environment())   # unpack all w9 scalars into local scope
cat("  w9 shared state loaded.\n")
OUT_DIR  <- file.path(ATM_DIR, "eof_analysis")               # w10c outputs
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Index CSV directories (for SPEI) ─────────────────────────────────────────
source(file.path(WD_PATH, "DROUGHT_ANALYSIS_utils.R"))      # load_basin_avg_csv
source(file.path(WD_PATH, "utils_teleconnection_addon.R"))  # load_teleconnection

# ── North Pacific EOF domain (0–360 longitude convention, matching ERA5 raw) ──
# 140°E to 240°E (= 120°W), 20°N–65°N
# Rationale: captures PDO, PNA, and ENSO teleconnection patterns affecting BC.
EOF_LON_MIN <- 140    # 140°E
EOF_LON_MAX <- 240    # 240°E = 120°W
EOF_LAT_MIN <-  20    # 20°N
EOF_LAT_MAX <-  65    # 65°N
# In wrapped (-180:+180) convention for map plotting:
EOF_LON_MIN_WGS <- -180  # wraps from 180°E onward
EOF_LON_MAX_WGS <- -120  # 120°W in -180:+180

# ── EOF analysis settings ─────────────────────────────────────────────────────
N_PCS_MAX  <- 10L    # max PCs to extract; displayed in scree; retained by Kaiser
# Seasons to analyse separately
SEASONS <- list(
  all = NULL,             # all 12 months
  JJA = c(6L, 7L, 8L),   # drought season
  DJF = c(12L, 1L, 2L)   # teleconnection season
)

# ── 2022–2025 event window ────────────────────────────────────────────────────
# Inherited from w9: DROUGHT_FOCUS_START, DROUGHT_FOCUS_END

# ── Correlation / regression settings ────────────────────────────────────────
ALPHA_CORR  <- 0.05       # significance level for correlations
MAX_LAG_MO  <- 2L         # maximum lag (months) for PC–index lag analysis
NW_BANDWIDTH <- NULL      # NULL = automatic Newey-West bandwidth selection
SPEI_SCALES_MAIN <- c(1L, 2L, 3L)   # SPEI (Penman-Monteith) accumulation scales
SPI_SCALES_MAIN  <- c(1L, 2L, 3L)   # SPI accumulation scales

# ── Figure DPI (use w9 constant; fallback if not in shared state) ─────────────
FIG_DPI  <- if (exists("FIGURE_DPI")) FIGURE_DPI else 300L
FIG_W    <- if (exists("FIGURE_WIDTH_STD"))  FIGURE_WIDTH_STD  else 12
FIG_W_WD <- if (exists("FIGURE_WIDTH_WIDE")) FIGURE_WIDTH_WIDE else 16

cat("  EOF domain: [", EOF_LON_MIN, "–", EOF_LON_MAX, "°lon,",
    EOF_LAT_MIN, "–", EOF_LAT_MAX, "°lat] (0-360 convention)\n")
cat("  N PCs max:", N_PCS_MAX, " | Seasons:", paste(names(SEASONS), collapse=", "), "\n\n")

# ============================================================================
#  SECTION 2: HELPER FUNCTIONS
# ============================================================================
cat("============================================================\n")
cat(" SECTION 2: Defining helper functions\n")
cat("============================================================\n")

# ── 2-A: Longitude wrap ───────────────────────────────────────────────────────
#' Wrap a data.frame x column from 0:360 convention to -180:180 for ggplot2.
wrap_lon <- function(x_vec) {
  ifelse(x_vec > 180, x_vec - 360, x_vec)
}

# ── 2-B: Effective sample size (Bretherton et al. 1999) ──────────────────────
#' Estimate effective sample size n_eff = n * (1 - r1) / (1 + r1)
#' where r1 is the average lag-1 autocorrelation of the two series.
#' A floor of 3 is enforced to avoid degenerate tests.
effective_n <- function(x, y) {
  n  <- sum(!is.na(x) & !is.na(y))
  r1 <- 0.5 * (cor(x[-length(x)], x[-1], use = "complete.obs") +
                 cor(y[-length(y)], y[-1], use = "complete.obs"))
  r1 <- max(min(r1, 0.99), -0.99)   # clamp to avoid div/0
  n_eff <- n * (1 - r1) / (1 + r1)
  max(round(n_eff), 3L)
}

# ── 2-C: Correlation with effective-n significance test ──────────────────────
#' Pearson correlation between x and y with p-value adjusted for autocorrelation
#' using the effective sample size approach of Bretherton et al. (1999).
cor_test_effn <- function(x, y) {
  ok  <- !is.na(x) & !is.na(y)
  r   <- cor(x[ok], y[ok])
  n   <- sum(ok)
  n_e <- effective_n(x[ok], y[ok])
  t_stat <- r * sqrt((n_e - 2) / (1 - r^2))
  p_val  <- 2 * pt(-abs(t_stat), df = n_e - 2)
  list(r = r, p = p_val, n = n, n_eff = n_e, t = t_stat)
}

# ── 2-D: Cosine-latitude weighting for PCA ───────────────────────────────────
#' Apply sqrt(cos(lat)) weighting to each column of a time × space matrix.
#' Columns with weight = 0 (poles) or NA are excluded.
apply_cosine_weights <- function(mat, lat_vec) {
  w <- sqrt(pmax(cos(lat_vec * pi / 180), 0))
  ok <- w > 0 & !is.na(w)
  list(mat_w = sweep(mat[, ok, drop = FALSE], 2, w[ok], `*`),
       cell_idx = which(ok),
       weights  = w[ok])
}

# ── 2-E: Core EOF / PCA function ─────────────────────────────────────────────
#' Compute EOF analysis on a SpatRaster cropped to the EOF domain.
#'
#' Algorithm:
#'   1. Remove all-NA cells (e.g. SST land mask).
#'   2. Standardise each valid cell's time series to zero mean, unit variance.
#'   3. Apply sqrt(cos(lat)) weights to the standardised columns.
#'   4. Run prcomp() (SVD-based, numerically stable).
#'   5. Return EOFs (spatial), PCs (temporal), and eigenvalue diagnostics.
#'
#' @param rast_anom  SpatRaster: full-domain anomaly raster (time in layers)
#' @param dates      Date vector aligned to rast_anom layers
#' @param month_sel  Integer vector of months to subset (NULL = all months)
#' @param n_pcs      Maximum number of PCs to retain
#' @param field_name Character label for progress messages
#' @return Named list: scores [n_time × n_pcs], loadings [n_valid × n_pcs],
#'         var_pct [vector], dates_sub, cell_idx, lat_vec, lon_vec,
#'         tmpl_rast (template raster for spatial reconstruction)
compute_eof <- function(rast_anom, dates, month_sel = NULL,
                        n_pcs = N_PCS_MAX, field_name = "") {
  
  # ── Crop to North Pacific domain ──────────────────────────────────────────
  domain_ext <- terra::ext(EOF_LON_MIN, EOF_LON_MAX, EOF_LAT_MIN, EOF_LAT_MAX)
  r_domain   <- terra::crop(rast_anom, domain_ext)
  cat(sprintf("    [%s] Cropped to EOF domain: %d cells × %d layers\n",
              field_name, terra::ncell(r_domain), terra::nlyr(r_domain)))
  
  # ── Subset to selected months ─────────────────────────────────────────────
  if (!is.null(month_sel)) {
    m_vec  <- as.integer(format(dates, "%m"))
    t_idx  <- which(m_vec %in% month_sel)
    r_sub  <- terra::subset(r_domain, t_idx)
    d_sub  <- dates[t_idx]
  } else {
    r_sub  <- r_domain
    t_idx  <- seq_len(length(dates))
    d_sub  <- dates
  }
  n_time <- terra::nlyr(r_sub)
  cat(sprintf("    [%s] Season subset: %d time steps\n", field_name, n_time))
  
  # ── Extract matrix: rows = time, cols = space ─────────────────────────────
  # terra::values() returns [n_cells × n_layers]; transpose for [n_time × n_cells]
  mat_raw <- t(terra::values(r_sub))   # [n_time × n_cells]
  
  # ── Identify fully valid (non-NA throughout time) cells ───────────────────
  na_per_cell <- colSums(is.na(mat_raw))
  valid_cells <- which(na_per_cell == 0)
  if (length(valid_cells) < 10)
    stop(sprintf("[%s] Fewer than 10 valid cells after NA removal.", field_name))
  mat_valid <- mat_raw[, valid_cells, drop = FALSE]
  cat(sprintf("    [%s] Valid cells: %d / %d\n",
              field_name, length(valid_cells), terra::ncell(r_domain)))
  
  # ── Get latitude coordinates for each valid cell ──────────────────────────
  xy_all   <- terra::xyFromCell(r_domain, seq_len(terra::ncell(r_domain)))
  lat_valid <- xy_all[valid_cells, "y"]
  lon_valid <- xy_all[valid_cells, "x"]
  
  # ── Standardise each cell (zero mean, unit variance) ─────────────────────
  cell_mean <- colMeans(mat_valid)
  cell_sd   <- apply(mat_valid, 2, sd)
  # Guard against zero-variance cells (should not occur on anomaly data)
  zero_var  <- cell_sd < .Machine$double.eps * 100
  if (any(zero_var)) {
    mat_valid  <- mat_valid[, !zero_var, drop = FALSE]
    lat_valid  <- lat_valid[!zero_var]
    lon_valid  <- lon_valid[!zero_var]
    valid_cells <- valid_cells[!zero_var]
    cell_mean  <- cell_mean[!zero_var]
    cell_sd    <- cell_sd[!zero_var]
    cat(sprintf("    [%s] Removed %d zero-variance cells.\n",
                field_name, sum(zero_var)))
  }
  mat_std <- sweep(sweep(mat_valid, 2, cell_mean, `-`), 2, cell_sd, `/`)
  
  # ── Apply cosine-latitude weighting ──────────────────────────────────────
  wt_res   <- apply_cosine_weights(mat_std, lat_valid)
  mat_w    <- wt_res$mat_w         # [n_time × n_valid_weighted]
  lat_used <- lat_valid[wt_res$cell_idx]
  lon_used <- lon_valid[wt_res$cell_idx]
  vc_used  <- valid_cells[wt_res$cell_idx]
  n_pcs_eff <- min(n_pcs, ncol(mat_w) - 1L, n_time - 1L)
  
  # ── SVD via prcomp ────────────────────────────────────────────────────────
  cat(sprintf("    [%s] Running prcomp (n_time=%d, n_space=%d, n_pcs=%d)...\n",
              field_name, n_time, ncol(mat_w), n_pcs_eff))
  pca_res  <- prcomp(mat_w, center = FALSE, scale. = FALSE,
                     rank. = n_pcs_eff)
  
  # ── Eigenvalues & variance explained ─────────────────────────────────────
  # prcomp stores standard deviations; eigenvalues = sdev^2
  # ── Eigenvalues & variance explained ─────────────────────────────────────
  # prcomp stores standard deviations; eigenvalues = sdev^2.
  # IMPORTANT: with a wide matrix (n_space >> n_time), prcomp's internal
  # cross-product path may return all min(n_time-1, n_space-1) eigenvalues
  # regardless of the rank. parameter.  We:
  #   (a) keep the full eigenvalue vector for the Kaiser criterion (which must
  #       be assessed against the full eigenvalue spectrum, not a truncated one)
  #   (b) truncate to n_pcs_eff for all stored outputs so downstream data.frame()
  #       calls never encounter a length mismatch (1 × n_pcs vs 912).
  eigenvals_full <- pca_res$sdev^2               # all eigenvalues prcomp computed
  n_ev           <- length(eigenvals_full)        # may be > n_pcs_eff
  
  # Kaiser criterion: computed on the FULL spectrum before truncation so that
  # kaiser_n reflects the true number of meaningful modes (not just the top 10).
  kaiser_n <- sum(eigenvals_full > mean(eigenvals_full))
  
  # Truncate to the n_pcs_eff components that were actually extracted
  n_keep    <- min(n_pcs_eff, n_ev)
  eigenvals <- eigenvals_full[seq_len(n_keep)]   # [n_pcs_eff]
  
  total_var_approx <- ncol(mat_w)   # each col has ~unit variance (trace approximation)
  var_pct   <- (eigenvals / total_var_approx) * 100   # [n_pcs_eff]
  cum_var   <- cumsum(var_pct)                         # [n_pcs_eff]
  
  # ── PC scores (temporal) — unweighted for physical interpretability ───────
  # pca_res$x gives [n_time × n_pcs] scores in the WEIGHTED space.
  # Store these directly: they are proportional to the projection of the
  # (weighted) anomaly field onto the EOF patterns.
  scores <- pca_res$x   # [n_time × n_pcs]
  
  # Normalise PC scores to unit variance for comparability across fields
  scores_sd <- apply(scores, 2, sd)
  scores_norm <- sweep(scores, 2, scores_sd, `/`)  # unit-variance PCs
  
  # ── EOF spatial patterns (loadings in weighted space) ─────────────────────
  # pca_res$rotation: [n_space_weighted × n_pcs] — these are the EOFs in
  # the cosine-weighted standardised space.  For spatial plotting we map
  # them back to the full raster grid.
  loadings <- pca_res$rotation   # [n_weighted_cells × n_pcs]
  
  # ── Sign convention: PC1 positive during 2022-2025 drought ───────────────
  # Check mean PC score during drought period; flip sign if negative, so
  # positive PC1 consistently denotes the drought-associated pattern.
  focus_months <- which(d_sub >= as.Date(sprintf("%d-01-01", DROUGHT_FOCUS_START)) &
                          d_sub <= as.Date(sprintf("%d-12-31", DROUGHT_FOCUS_END)))
  if (length(focus_months) >= 3) {
    for (k in seq_len(n_pcs_eff)) {
      if (mean(scores_norm[focus_months, k], na.rm = TRUE) < 0) {
        scores_norm[, k] <- -scores_norm[, k]
        loadings[, k]    <- -loadings[, k]
      }
    }
  }
  
  # ── Kaiser criterion: PCs with eigenvalue > 1 ────────────────────────────
  # When data are standardised, eigenvalues are relative to ncol(mat_w).
  # "Eigenvalue > 1" here is interpreted as "explains > 100%/ncol variance",
  # equivalent to the Kaiser criterion.  We scale eigenvalues to the number
  # of variables used.
  # NOTE: the strict Kaiser criterion requires the full covariance matrix.
  # With SVD on truncated rank, we report an approximation and advise
  # inspection of the scree diagram.
  cat(sprintf("    [%s] Kaiser criterion (eigenval > mean, full spectrum): %d PCs\n",
              field_name, kaiser_n))
  cat(sprintf("    [%s] Variance explained (top %d PCs): %s%%\n",
              field_name, n_keep,
              paste(round(var_pct[1:min(5, n_keep)], 1), collapse = ", ")))
  
  # ── Template raster for spatial reconstruction ────────────────────────────
  # A single-layer template of the domain at the EOF resolution
  tmpl <- terra::subset(r_domain, 1L)
  terra::values(tmpl) <- NA_real_
  
  list(
    field        = field_name,
    scores       = scores_norm,          # [n_time × n_pcs] normalised
    scores_raw   = scores,               # [n_time × n_pcs] raw
    loadings     = loadings,             # [n_weighted_cells × n_pcs]
    eigenvals    = eigenvals,            # [n_pcs_eff] — truncated to extracted PCs
    eigenvals_full = eigenvals_full,     # [all]       — full spectrum for diagnostics
    var_pct      = var_pct,              # [n_pcs_eff]
    cum_var      = cum_var,              # [n_pcs_eff]
    kaiser_n     = kaiser_n,             # scalar — based on full eigenvalue spectrum
    n_pcs        = n_keep,              # = n_pcs_eff (number of extracted PCs)
    dates_sub    = d_sub,
    t_idx        = t_idx,               # index into full date vector
    vc_used      = vc_used,             # cell indices in r_domain
    lat_used     = lat_used,
    lon_used     = lon_used,
    weights      = wt_res$weights,
    tmpl_rast    = tmpl
  )
}

# ── 2-F: Reconstruct EOF spatial pattern as SpatRaster ───────────────────────
#' Map EOF loadings back onto the domain raster template.
eof_to_raster <- function(eof_res, pc_k = 1L) {
  tmpl <- eof_res$tmpl_rast
  vals <- rep(NA_real_, terra::ncell(tmpl))
  vals[eof_res$vc_used] <- eof_res$loadings[, pc_k]
  terra::values(tmpl) <- vals
  tmpl
}

# ── 2-G: PC–SPEI lag correlation table ──────────────────────────────────────
#' Pearson r (with effective-n p-value) between a PC and a SPEI series
#' at lags 0:max_lag_mo.  Returns a data.frame.
pc_spei_lag_cor <- function(pc_ts, spei_ts, pc_dates, spei_dates,
                            max_lag = MAX_LAG_MO, label_pc = "PC1") {
  # Align on common dates at lag = 0
  merged0 <- data.frame(date = pc_dates, pc = pc_ts) |>
    inner_join(data.frame(date = spei_dates, spei = spei_ts), by = "date") |>
    arrange(date)
  
  bind_rows(lapply(0:max_lag, function(lag) {
    if (lag == 0) {
      x <- merged0$pc; y <- merged0$spei
    } else {
      # PC leads SPEI by `lag` months: PC at t, SPEI at t + lag
      n  <- nrow(merged0)
      x  <- merged0$pc[1:(n - lag)]
      y  <- merged0$spei[(lag + 1):n]
    }
    ok <- !is.na(x) & !is.na(y)
    if (sum(ok) < 20) return(NULL)
    ct <- cor_test_effn(x[ok], y[ok])
    data.frame(pc = label_pc, lag_mo = lag, r = ct$r,
               p_effn = ct$p, n = ct$n, n_eff = ct$n_eff)
  }))
}

# ── 2-H: Newey–West regression of SPEI on multiple PCs ───────────────────────
#' Multiple regression SPEI ~ PC1 + PC2 + ... with Newey-West robust SEs.
#' Returns a tidy data.frame of coefficients + diagnostics.
nw_regression <- function(spei_vec, pc_mat, pc_names,
                          bandwidth = NW_BANDWIDTH) {
  df  <- as.data.frame(pc_mat)
  colnames(df) <- pc_names
  df$spei <- spei_vec
  df <- na.omit(df)
  if (nrow(df) < (length(pc_names) + 2))
    return(NULL)
  
  frm   <- as.formula(paste("spei ~", paste(pc_names, collapse = " + ")))
  fit   <- lm(frm, data = df)
  nw_se <- coeftest(fit, vcov = sandwich::NeweyWest(fit, lag = bandwidth,
                                                    prewhite = FALSE))
  
  tidy_nw <- data.frame(
    term      = rownames(nw_se),
    estimate  = nw_se[, "Estimate"],
    std_error = nw_se[, "Std. Error"],
    t_stat    = nw_se[, "t value"],
    p_nw      = nw_se[, "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
  
  s   <- summary(fit)
  adj_r2 <- s$adj.r.squared
  
  list(coefs = tidy_nw, adj_r2 = adj_r2, fit = fit)
}

# ── 2-I: Historical ranking of event-period running mean ─────────────────────
#' Compute N-year running mean of a monthly time series and rank the focus
#' period window against all overlapping windows of the same length.
rank_event_period <- function(monthly_ts, monthly_dates,
                              focus_start, focus_end,
                              window_years = NULL) {
  # Default window = duration of focus period
  if (is.null(window_years))
    window_years <- as.integer(focus_end - focus_start + 1L)
  w_months <- window_years * 12L
  n <- length(monthly_ts)
  
  # N-year running means (non-overlapping approach: mean over each window)
  # Use a rolling window aligned to calendar year starts
  years_all <- as.integer(format(monthly_dates, "%Y"))
  unique_yrs <- sort(unique(years_all))
  
  means_list <- list()
  for (i in seq_len(length(unique_yrs) - window_years + 1)) {
    yr_start <- unique_yrs[i]
    yr_end   <- yr_start + window_years - 1L
    idx_w    <- which(years_all >= yr_start & years_all <= yr_end)
    if (length(idx_w) < (w_months * 0.9)) next  # skip if < 90% data
    means_list[[length(means_list) + 1]] <- data.frame(
      yr_start   = yr_start,
      yr_end     = yr_end,
      mean_pc    = mean(monthly_ts[idx_w], na.rm = TRUE),
      n_months   = length(idx_w)
    )
  }
  all_means <- bind_rows(means_list)
  
  # Focus period value
  focus_idx <- which(years_all >= focus_start & years_all <= focus_end)
  focus_val <- mean(monthly_ts[focus_idx], na.rm = TRUE)
  
  # Rank (1 = highest)
  all_means$rank <- rank(-all_means$mean_pc, ties.method = "first")
  focus_row  <- all_means[all_means$yr_start == focus_start, ]
  focus_rank <- if (nrow(focus_row) > 0) focus_row$rank[1] else NA_integer_
  
  # Z-score relative to full distribution
  focus_z <- (focus_val - mean(all_means$mean_pc, na.rm = TRUE)) /
    sd(all_means$mean_pc, na.rm = TRUE)
  
  list(all_means  = all_means,
       focus_val  = focus_val,
       focus_rank = focus_rank,
       focus_z    = focus_z,
       n_windows  = nrow(all_means),
       window_yr  = window_years)
}

# ── 2-J: GEV return period (empirical fallback if evd unavailable) ────────────
#' Fit a GEV to annual block maxima of |PC values| and estimate return period
#' for the focus-period mean.  Uses evd::fgev if available; otherwise falls
#' back to the empirical Weibull plotting position.
gev_return_period <- function(monthly_ts, monthly_dates,
                              focus_start, focus_end) {
  years_all <- as.integer(format(monthly_dates, "%Y"))
  
  # Annual maxima of ABSOLUTE value (most extreme year in each direction)
  ann_max <- tapply(monthly_ts, years_all, function(x) max(x, na.rm = TRUE))
  ann_max <- sort(as.numeric(ann_max), decreasing = TRUE)
  
  focus_idx <- which(years_all >= focus_start & years_all <= focus_end)
  focus_val <- mean(monthly_ts[focus_idx], na.rm = TRUE)
  
  if (requireNamespace("evd", quietly = TRUE)) {
    fit <- tryCatch(
      evd::fgev(ann_max),
      error = function(e) NULL)
    if (!is.null(fit)) {
      # Return period from GEV: T = 1 / (1 - F(x))
      T_gev <- tryCatch({
        F_x <- evd::pgev(focus_val,
                         loc   = fit$estimate["loc"],
                         scale = fit$estimate["scale"],
                         shape = fit$estimate["shape"])
        1 / (1 - F_x)
      }, error = function(e) NA_real_)
      return(list(method = "GEV (evd::fgev)", T_years = T_gev))
    }
  }
  
  # Empirical Weibull plotting position (rank-based)
  n <- length(ann_max)
  ranks_emp <- rank(-ann_max, ties.method = "first")   # 1 = largest
  T_emp <- (n + 1) / ranks_emp
  closest <- which.min(abs(ann_max - focus_val))
  list(method   = "Empirical (Weibull)",
       T_years  = T_emp[closest])
}

cat("  All helper functions defined.\n\n")

# ============================================================================
#  SECTION 3: LOAD ERA5 ANOMALIES + SPEI DATA
# ============================================================================
cat("============================================================\n")
cat(" SECTION 3: Loading ERA5 anomalies and SPEI indices\n")
cat("============================================================\n")

# ── ERA5 anomaly NetCDFs (from w9) ────────────────────────────────────────────
req_nc <- function(path) {
  if (!file.exists(path))
    stop("Anomaly file not found:\n  ", path,
         "\n  Run w9_atmospheric_diagnostics.R first.")
  cat(sprintf("  [load] %s\n", basename(path)))
  terra::rast(path)
}
z500_anom <- req_nc(file.path(DATA_DIR, sprintf(
  "z500_monthly_anomaly_clim%d_%d.nc", CLIM_START, CLIM_END)))
slp_anom  <- req_nc(file.path(DATA_DIR, sprintf(
  "slp_monthly_anomaly_clim%d_%d.nc",  CLIM_START, CLIM_END)))
sst_anom  <- req_nc(file.path(DATA_DIR, sprintf(
  "sst_monthly_anomaly_clim%d_%d.nc",  CLIM_START, CLIM_END)))

# ── Longitude convention normalisation ────────────────────────────────────────
#
#  ERA5 raw files use 0–360 longitudes internally in terra (confirmed: DOM_XMIN
#  = 100, DOM_XMAX = 250 in w9/w10a Section 6).  However, when anomaly fields
#  are written with terra::writeCDF() and re-read, terra (via the underlying
#  netCDF library) may re-encode the x-axis as –180:+180 (confirmed by the
#  ms1_atmospheric_blanks.csv output, which shows Z500 ridge lon = –176.75 °).
#
#  The EOF domain is defined in the 0–360 convention (EOF_LON_MIN = 140,
#  EOF_LON_MAX = 240).  terra::crop() fails with "[crop] extents do not overlap"
#  when the raster is in –180:+180 and the crop extent is in 0–360, because
#  the 140–240 rectangle lies entirely outside the –180:+180 extent.
#
#  Fix: detect the convention from the x-extent of the first loaded raster and,
#  if –180:+180, rotate all three anomaly rasters to 0–360 with terra::rotate().
#  terra::rotate() reorders columns so that longitude 0° is at the left edge,
#  converting a –180:+180 raster to 0–360.  All subsequent EOF crop calls then
#  work as coded.
#
#  This rotation is applied ONLY inside w10c (local copies); the original
#  anomaly NetCDF files on disk are not modified.

.lon_max <- terra::ext(z500_anom)[2]   # xmax of first raster layer
if (.lon_max <= 181) {
  # Raster is in –180:+180 convention → rotate to 0–360
  cat("  [INFO] Anomaly rasters are in –180:+180 convention (xmax =",
      round(.lon_max, 1), "). Rotating to 0–360 for EOF domain crop.\n")
  z500_anom <- terra::rotate(z500_anom)
  slp_anom  <- terra::rotate(slp_anom)
  sst_anom  <- terra::rotate(sst_anom)
  cat(sprintf("  [INFO] After rotation: Z500 x-extent [%.1f, %.1f]\n",
              terra::ext(z500_anom)[1], terra::ext(z500_anom)[2]))
} else {
  cat("  [INFO] Anomaly rasters already in 0–360 convention (xmax =",
      round(.lon_max, 1), "). No rotation needed.\n")
}

# ── Full date spine ────────────────────────────────────────────────────────────
all_dates <- seq.Date(as.Date(sprintf("%d-01-01", START_YEAR)),
                      as.Date(sprintf("%d-12-01", END_YEAR)),
                      by = "month")
n_total <- length(all_dates)
stopifnot(terra::nlyr(z500_anom) == n_total)
cat(sprintf("  Date spine: %s – %s (%d months)\n",
            min(all_dates), max(all_dates), n_total))

# ── SPEI (Penman-Monteith) basin-averaged time series ─────────────────────────
cat("\n  Loading SPEI (PM) basin means...\n")
spei_list <- lapply(SPEI_SCALES_MAIN, function(sc) {
  df <- tryCatch(load_basin_avg_csv("spei", sc),
                 error = function(e) NULL,
                 warning = function(w) suppressWarnings(load_basin_avg_csv("spei", sc)))
  if (is.null(df) || nrow(df) == 0) {
    cat(sprintf("  ⚠ SPEI-%d: not found, skipped\n", sc)); return(NULL)
  }
  df$scale <- sc
  df$index <- sprintf("SPEI%d", sc)
  df
})
names(spei_list) <- sprintf("SPEI%d", SPEI_SCALES_MAIN)
spei_list <- Filter(Negate(is.null), spei_list)
cat(sprintf("  SPEI scales loaded: %s\n", paste(names(spei_list), collapse = ", ")))

# ── SPI basin-averaged time series ────────────────────────────────────────────
cat("\n  Loading SPI basin means...\n")
spi_list <- lapply(SPI_SCALES_MAIN, function(sc) {
  df <- tryCatch(load_basin_avg_csv("spi", sc),
                 error = function(e) NULL,
                 warning = function(w) suppressWarnings(load_basin_avg_csv("spi", sc)))
  if (is.null(df) || nrow(df) == 0) {
    cat(sprintf("  ⚠ SPI-%d: not found, skipped\n", sc)); return(NULL)
  }
  df$scale <- sc
  df$index <- sprintf("SPI%d", sc)
  df
})
names(spi_list) <- sprintf("SPI%d", SPI_SCALES_MAIN)
spi_list <- Filter(Negate(is.null), spi_list)
cat(sprintf("  SPI scales loaded: %s\n", paste(names(spi_list), collapse = ", ")))

# ── Combined drought index list (SPEI-PM + SPI) ───────────────────────────────
index_list <- c(spei_list, spi_list)
cat(sprintf("  Drought indices available: %s\n",
            paste(names(index_list), collapse = ", ")))

# ── Standard teleconnection indices (from cache) ──────────────────────────────
cat("\n  Loading standard teleconnection indices...\n")
tele_names <- c("pdo", "pna", "nino34", "ao")
tele_list  <- lapply(setNames(tele_names, toupper(tele_names)), function(nm) {
  tryCatch(
    load_teleconnection(nm, start_year = START_YEAR, end_year = END_YEAR),
    error = function(e) { cat(sprintf("  ⚠ %s: %s\n", nm, e$message)); NULL })
})
tele_list <- Filter(Negate(is.null), tele_list)
cat(sprintf("  Indices loaded: %s\n\n", paste(names(tele_list), collapse = ", ")))

# ============================================================================
#  SECTION 4: EOF/PCA ANALYSIS
#  Runs for each field × season combination; stores results in eof_results.
# ============================================================================
cat("============================================================\n")
cat(" SECTION 4: EOF / PCA analysis\n")
cat("============================================================\n")

eof_results <- list()   # keyed as "{field}_{season}"

fields <- list(
  Z500 = z500_anom,
  SLP  = slp_anom,
  SST  = sst_anom
)

for (fld in names(fields)) {
  for (seas_nm in names(SEASONS)) {
    key <- sprintf("%s_%s", fld, seas_nm)
    cat(sprintf("\n  --- EOF: %s | Season: %s ---\n", fld, seas_nm))
    
    result <- tryCatch(
      compute_eof(
        rast_anom  = fields[[fld]],
        dates      = all_dates,
        month_sel  = SEASONS[[seas_nm]],
        n_pcs      = N_PCS_MAX,
        field_name = key),
      error = function(e) {
        cat(sprintf("  ✗ EOF failed for %s: %s\n", key, e$message))
        NULL
      })
    
    if (!is.null(result)) {
      eof_results[[key]] <- result
      cat(sprintf("  ✓ %s complete: %d PCs, top-5 var = [%s]%%\n",
                  key,
                  result$n_pcs,
                  paste(round(result$var_pct[1:min(5, result$n_pcs)], 1),
                        collapse = ", ")))
    }
  }
}

cat(sprintf("\n  EOF analysis complete: %d field×season combinations.\n\n",
            length(eof_results)))

# [FIX] Diagnostic guard: warn explicitly when all EOF combinations failed so
# the root cause (missing w9 NetCDFs, layer-count mismatch, etc.) is visible
# immediately rather than only manifesting as a cryptic error in Section 11.
if (length(eof_results) == 0) {
  warning(
    "eof_results is empty — all EOF/PCA combinations failed.\n",
    "  Likely causes:\n",
    "    1. Anomaly NetCDFs not found in DATA_DIR (run w9 first).\n",
    "    2. Layer count mismatch between NetCDF and all_dates.\n",
    "    3. Fewer than 10 valid (non-NA) cells in the EOF domain.\n",
    "  Sections 6–9 will produce no output. Section 11 will write an\n",
    "  empty ms2_eof_blanks.csv. Check messages above for ✗ EOF failed lines.",
    call. = FALSE
  )
}

# ============================================================================
#  SECTION 5: SAVE EOF STATE
# ============================================================================
cat("============================================================\n")
cat(" SECTION 5: Saving EOF state RDS\n")
cat("============================================================\n")

# ── PC time series CSV: all-month analyses for the three fields ───────────────
pc_ts_rows <- list()
for (fld in names(fields)) {
  key    <- sprintf("%s_all", fld)
  if (!key %in% names(eof_results)) next
  res    <- eof_results[[key]]
  n_p    <- min(res$n_pcs, 5L)   # export up to 5 PCs per field
  for (k in seq_len(n_p)) {
    pc_ts_rows[[length(pc_ts_rows) + 1]] <- data.frame(
      date     = res$dates_sub,
      field    = fld,
      pc       = sprintf("PC%d", k),
      value    = res$scores[, k],
      var_pct  = res$var_pct[k],
      season   = "all"
    )
  }
}
pc_ts_df <- bind_rows(pc_ts_rows)
write_csv(pc_ts_df, file.path(OUT_DIR, "eof_pc_timeseries.csv"))
cat(sprintf("  PC time series CSV: %d rows → eof_pc_timeseries.csv\n",
            nrow(pc_ts_df)))

# ── Full EOF state RDS ────────────────────────────────────────────────────────
w10c_state <- list(
  eof_results  = eof_results,
  pc_ts_df     = pc_ts_df,
  all_dates    = all_dates,
  EOF_domain   = list(lon_min = EOF_LON_MIN, lon_max = EOF_LON_MAX,
                      lat_min = EOF_LAT_MIN, lat_max = EOF_LAT_MAX),
  config       = list(N_PCS_MAX = N_PCS_MAX, SEASONS = SEASONS,
                      CLIM_START = CLIM_START, CLIM_END = CLIM_END,
                      DROUGHT_FOCUS_START = DROUGHT_FOCUS_START,
                      DROUGHT_FOCUS_END   = DROUGHT_FOCUS_END)
)
rds_w10c <- file.path(OUT_DIR, "w10c_eof_state.rds")
saveRDS(w10c_state, rds_w10c)
cat(sprintf("  EOF state RDS: %s\n\n", basename(rds_w10c)))

# ============================================================================
#  SECTION 6: PC–SPEI CORRELATION AND VARIANCE PARTITIONING
#  Implements framework Steps 3.1–3.4 (new Section 4.7.1)
# ============================================================================
cat("============================================================\n")
cat(" SECTION 6: PC–SPEI correlation and regression\n")
cat("============================================================\n")

corr_rows    <- list()   # for Table 4
reg_rows     <- list()   # for Table 5
ms_blanks    <- list()   # manuscript text fill-in values

# Primary scales for variance-partitioning regression (one per index type)
PRIMARY_SCALES <- c(
  SPEI = sprintf("SPEI%d", if (3L %in% SPEI_SCALES_MAIN) 3L else SPEI_SCALES_MAIN[length(SPEI_SCALES_MAIN)]),
  SPI  = sprintf("SPI%d",  if (3L %in% SPI_SCALES_MAIN)  3L else SPI_SCALES_MAIN[length(SPI_SCALES_MAIN)])
)
PRIMARY_SCALES <- PRIMARY_SCALES[PRIMARY_SCALES %in% names(index_list)]

for (fld in names(fields)) {
  key_all <- sprintf("%s_all", fld)
  if (!key_all %in% names(eof_results)) next
  res <- eof_results[[key_all]]
  n_p <- min(res$kaiser_n + 1L, res$n_pcs, 4L)  # top PCs to correlate
  
  for (index_nm in names(index_list)) {
    index_df  <- index_list[[index_nm]]
    
    for (k in seq_len(n_p)) {
      pc_vals   <- res$scores[, k]
      pc_dates  <- res$dates_sub
      
      # Lag correlation table
      lc <- tryCatch(
        pc_spei_lag_cor(pc_vals, index_df$value, pc_dates, index_df$date,
                        max_lag   = MAX_LAG_MO,
                        label_pc  = sprintf("%s-PC%d", fld, k)),
        error = function(e) NULL)
      if (!is.null(lc)) corr_rows <- c(corr_rows, list(
        lc %>% mutate(field = fld, index_nm = index_nm)))
    }
  }
  
  # ── Multiple regression: Index ~ PC1 + PC2 + PC3 (for each primary scale) ──
  for (prim_key in PRIMARY_SCALES) {
    if (!prim_key %in% names(index_list)) next
    index_df <- index_list[[prim_key]]
    n_p_reg <- min(n_p, 3L)
    pc_labels <- sprintf("%s_PC%d", fld, seq_len(n_p_reg))
    
    # Align PC scores to index dates (common dates)
    pc_wide <- data.frame(date = res$dates_sub,
                          do.call(cbind, lapply(seq_len(n_p_reg), function(k)
                            res$scores[, k])))
    colnames(pc_wide)[-1] <- pc_labels
    merged_reg <- index_df %>%
      select(date, spei = value) %>%
      inner_join(pc_wide, by = "date") %>%
      arrange(date)
    
    if (nrow(merged_reg) < 30) next
    
    reg_res <- tryCatch(
      nw_regression(merged_reg$spei,
                    as.matrix(merged_reg[, pc_labels]),
                    pc_labels),
      error = function(e) { cat("  ⚠ NW regression failed:", e$message, "\n"); NULL })
    
    if (!is.null(reg_res)) {
      reg_row <- reg_res$coefs %>%
        mutate(field = fld, index_nm = prim_key,
               adj_r2 = reg_res$adj_r2)
      reg_rows <- c(reg_rows, list(reg_row))
      
      cat(sprintf("  [%s ~ %s] adj-R² = %.3f\n",
                  fld, prim_key, reg_res$adj_r2))
      
      # Manuscript blank: field-level adj-R²
      ms_blanks[[sprintf("reg_adjR2_%s_%s", fld, prim_key)]] <-
        round(reg_res$adj_r2 * 100, 1)
    }
  }
}

# Bind and save correlation table
if (length(corr_rows)) {
  corr_df <- bind_rows(corr_rows)
  write_csv(corr_df, file.path(OUT_DIR, "Table4_PC_correlations.csv"))
  cat(sprintf("  Table 4 saved: %d correlation rows\n", nrow(corr_df)))
} else {
  corr_df <- NULL
  cat("  ⚠ No correlation rows computed.\n")
}

# Bind and save regression table
if (length(reg_rows)) {
  reg_df <- bind_rows(reg_rows)
  write_csv(reg_df, file.path(OUT_DIR, "Table5_regression.csv"))
  cat(sprintf("  Table 5 saved: %d regression term rows\n", nrow(reg_df)))
} else {
  reg_df <- NULL
  cat("  ⚠ No regression rows computed.\n")
}

# ============================================================================
#  SECTION 7: COMPARISON WITH STANDARD TELECONNECTION INDICES
#  Validates physical interpretability of EOFs (framework Step 3.5)
# ============================================================================
cat("\n============================================================\n")
cat(" SECTION 7: PC vs. standard teleconnection index comparison\n")
cat("============================================================\n")

tele_corr_rows <- list()

for (fld in names(fields)) {
  key_all <- sprintf("%s_all", fld)
  if (!key_all %in% names(eof_results)) next
  res <- eof_results[[key_all]]
  n_p <- min(3L, res$n_pcs)
  
  for (tele_nm in names(tele_list)) {
    tele_df   <- tele_list[[tele_nm]]
    # Rename value column if needed
    if (!"value" %in% names(tele_df)) {
      val_col <- setdiff(names(tele_df), "date")
      tele_df <- tele_df %>% rename(value = all_of(val_col[1]))
    }
    
    for (k in seq_len(n_p)) {
      pc_df <- data.frame(date = res$dates_sub, pc = res$scores[, k])
      mg    <- inner_join(pc_df, tele_df %>% select(date, value), by = "date")
      if (nrow(mg) < 30) next
      ct    <- cor_test_effn(mg$pc, mg$value)
      tele_corr_rows[[length(tele_corr_rows) + 1]] <- data.frame(
        field  = fld, pc = sprintf("PC%d", k),
        tele   = tele_nm,
        r      = round(ct$r, 3),
        p_effn = round(ct$p, 4),
        n_eff  = ct$n_eff
      )
    }
  }
}

if (length(tele_corr_rows)) {
  tele_corr_df <- bind_rows(tele_corr_rows)
  write_csv(tele_corr_df,
            file.path(OUT_DIR, "Table4_PC_vs_teleconnections.csv"))
  cat("  Table 4 (teleconnections) saved.\n")
  
  # Best matching index per PC (for physical interpretability)
  best_match <- tele_corr_df %>%
    group_by(field, pc) %>%
    slice_max(abs(r), n = 1, with_ties = FALSE) %>%
    ungroup()
  cat("\n  Best-matching standard index per PC:\n")
  print(best_match %>% select(field, pc, tele, r, p_effn))
  
  # Populate manuscript blanks
  for (i in seq_len(nrow(best_match))) {
    bm  <- best_match[i, ]
    key <- sprintf("pattern_cor_%s_%s_%s", bm$field, bm$pc, bm$tele)
    ms_blanks[[key]] <- bm$r
  }
}

# ============================================================================
#  SECTION 8: 2022–2025 EVENT QUANTIFICATION
#  Implements framework Step 4 / new Section 4.7.2
# ============================================================================
cat("\n============================================================\n")
cat(" SECTION 8: 2022–2025 event quantification\n")
cat("============================================================\n")

event_rows <- list()

for (fld in names(fields)) {
  key_all <- sprintf("%s_all", fld)
  if (!key_all %in% names(eof_results)) next
  res   <- eof_results[[key_all]]
  pc1   <- res$scores[, 1]
  dates <- res$dates_sub
  
  # ── 8-A: 2022-2025 mean, z-score, and historical ranking ─────────────────
  rk <- rank_event_period(pc1, dates,
                          focus_start   = DROUGHT_FOCUS_START,
                          focus_end     = DROUGHT_FOCUS_END,
                          window_years  = as.integer(DROUGHT_FOCUS_END -
                                                       DROUGHT_FOCUS_START + 1L))
  cat(sprintf("  [%s PC1] 2022–%d mean = %.2f | Rank %d of %d windows | z = %.2f\n",
              fld, DROUGHT_FOCUS_END,
              rk$focus_val, rk$focus_rank, rk$n_windows, rk$focus_z))
  
  # ── 8-B: GEV return period ───────────────────────────────────────────────
  gev_res <- tryCatch(
    gev_return_period(pc1, dates, DROUGHT_FOCUS_START, DROUGHT_FOCUS_END),
    error = function(e) NULL)
  T_years <- if (!is.null(gev_res)) gev_res$T_years else NA_real_
  method  <- if (!is.null(gev_res)) gev_res$method  else "failed"
  cat(sprintf("  [%s PC1] Return period (%s): %.0f years\n",
              fld, method, T_years))
  
  # ── Manuscript blanks ─────────────────────────────────────────────────────
  ms_blanks[[sprintf("event_mean_%s_PC1", fld)]]        <- round(rk$focus_val, 2)
  ms_blanks[[sprintf("event_zscore_%s_PC1", fld)]]      <- round(rk$focus_z, 2)
  ms_blanks[[sprintf("event_rank_%s_PC1", fld)]]        <- rk$focus_rank
  ms_blanks[[sprintf("event_n_windows_%s_PC1", fld)]]   <- rk$n_windows
  ms_blanks[[sprintf("event_return_yr_%s_PC1", fld)]]   <- round(T_years, 0)
  ms_blanks[[sprintf("varexp_%s_PC1", fld)]]            <-
    round(eof_results[[key_all]]$var_pct[1], 1)
  ms_blanks[[sprintf("kaiser_n_%s_all", fld)]]          <-
    eof_results[[key_all]]$kaiser_n
  
  event_rows[[length(event_rows) + 1]] <- data.frame(
    field      = fld,
    focus_mean = rk$focus_val,
    focus_z    = rk$focus_z,
    rank       = rk$focus_rank,
    n_windows  = rk$n_windows,
    T_return   = T_years,
    gev_method = method
  )
}

event_df <- bind_rows(event_rows)
write_csv(event_df, file.path(OUT_DIR, "event_quantification_2022_2025.csv"))
cat("\n  Event quantification CSV saved.\n")

# ============================================================================
#  SECTION 9: TREND ANALYSIS ON PC TIME SERIES
#  Implements framework Step 3.6 (parallel to thermodynamic fraction trends)
# ============================================================================
cat("\n============================================================\n")
cat(" SECTION 9: Trend analysis on PC time series\n")
cat("============================================================\n")

trend_rows <- list()

for (fld in names(fields)) {
  key_all <- sprintf("%s_all", fld)
  if (!key_all %in% names(eof_results)) next
  res <- eof_results[[key_all]]
  
  for (k in seq_len(min(3L, res$n_pcs))) {
    pc_vals <- res$scores[, k]
    dates_k <- res$dates_sub
    t_idx   <- as.numeric(dates_k - min(dates_k)) / 365.25   # years since start
    
    ok  <- !is.na(pc_vals)
    fit <- lm(pc_vals[ok] ~ t_idx[ok])
    nw  <- tryCatch(
      coeftest(fit, vcov = sandwich::NeweyWest(fit, prewhite = FALSE)),
      error = function(e) summary(fit)$coefficients)
    
    slope    <- nw[2, "Estimate"]
    se_slope <- nw[2, "Std. Error"]
    p_slope  <- nw[2, "Pr(>|t|)"]
    
    cat(sprintf("  [%s PC%d] trend = %.4f/yr ± %.4f (p=%.3f)%s\n",
                fld, k, slope, se_slope, p_slope,
                ifelse(p_slope < ALPHA_CORR, " *", "")))
    
    trend_rows[[length(trend_rows) + 1]] <- data.frame(
      field = fld, pc = sprintf("PC%d", k),
      slope_per_yr = slope, se = se_slope, p_nw = p_slope,
      significant  = p_slope < ALPHA_CORR
    )
  }
}

if (length(trend_rows)) {
  trend_df <- bind_rows(trend_rows)
  write_csv(trend_df, file.path(OUT_DIR, "pc_trend_analysis.csv"))
  cat("  PC trend table saved.\n")
}

# ============================================================================
#  SECTION 10: VISUALISATION
#  10-A: Scree diagrams
#  10-B: EOF spatial pattern maps
#  10-C: PC time series with SPEI overlay and 2022–2025 highlight
#  10-D: Correlation heatmap
# ============================================================================
cat("\n============================================================\n")
cat(" SECTION 10: Generating figures\n")
cat("============================================================\n")

# ── Shared colour palette (consistent with w9 / w10a) ────────────────────────
PAL_EOF  <- "RdBu"
COL_PC   <- "#1F4E79"   # dark navy for PC line
COL_SPEI <- "#922B21"   # dark red for SPEI line
COL_FOCUS <- "#E67E22"  # orange for 2022-2025 shading

# ── Shared time-series theme ─────────────────────────────────────────────────
theme_ts <- function(base_size = 10) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(colour = "grey92", linewidth = 0.4),
      axis.line.x       = element_line(colour = "grey40", linewidth = 0.4),
      plot.title        = element_text(face = "bold", size = base_size,
                                       hjust = 0),
      plot.subtitle     = element_text(size = base_size - 1, colour = "grey40",
                                       hjust = 0),
      legend.position   = "bottom",
      legend.key.size   = unit(0.4, "cm"),
      legend.text       = element_text(size = base_size - 1)
    )
}

# ── Helper: save pdf + png ────────────────────────────────────────────────────
save_fig <- function(p, stem, w, h) {
  pdf_path <- file.path(OUT_DIR, paste0(stem, ".pdf"))
  png_path <- file.path(OUT_DIR, paste0(stem, ".png"))
  tryCatch({
    ggplot2::ggsave(pdf_path, p, width = w, height = h, units = "in",
                    device = "pdf")
    cat(sprintf("  ✓ %s.pdf\n", stem))
  }, error = function(e) cat(sprintf("  ✗ %s.pdf: %s\n", stem, e$message)))
  tryCatch({
    ggplot2::ggsave(png_path, p, width = w, height = h, units = "in",
                    dpi = FIG_DPI, device = "png")
    cat(sprintf("  ✓ %s.png\n", stem))
  }, error = function(e) cat(sprintf("  ✗ %s.png: %s\n", stem, e$message)))
}

# ── 10-A: SCREE DIAGRAMS ─────────────────────────────────────────────────────
cat("\n  10-A: Scree diagrams\n")

scree_rows <- list()
for (fld in names(fields)) {
  key_all <- sprintf("%s_all", fld)
  if (!key_all %in% names(eof_results)) next
  res <- eof_results[[key_all]]
  scree_rows[[fld]] <- data.frame(
    field   = fld,
    pc      = seq_len(res$n_pcs),
    eigenval = res$eigenvals,
    var_pct  = res$var_pct,
    cum_var  = res$cum_var,
    kaiser   = res$eigenvals > mean(res$eigenvals)  # modified Kaiser
  )
}
scree_df <- bind_rows(scree_rows)

if (nrow(scree_df) == 0) {
  cat("  ⚠ No EOF results available for scree diagram — skipping.\n")
} else {
  
  p_scree <- ggplot(scree_df, aes(x = pc, y = var_pct)) +
    geom_col(aes(fill = kaiser), width = 0.7, colour = "white", linewidth = 0.3) +
    geom_line(aes(y = cum_var), colour = "grey30", linewidth = 0.7,
              linetype = "dashed") +
    geom_point(aes(y = cum_var), colour = "grey30", size = 1.5) +
    geom_hline(yintercept = 80, colour = "#E74C3C", linewidth = 0.5,
               linetype = "dotted", alpha = 0.7) +
    scale_fill_manual(values = c("TRUE" = "#1F4E79", "FALSE" = "#AED6F1"),
                      labels = c("TRUE" = "Retained (Kaiser)", "FALSE" = "Rejected"),
                      name   = NULL) +
    facet_wrap(~ field, ncol = 3) +
    scale_x_continuous(breaks = seq(1, N_PCS_MAX, 2)) +
    scale_y_continuous(
      name     = "Variance explained (%)",
      sec.axis = dup_axis(name = "Cumulative variance (%)")) +
    labs(title    = "Scree diagram — North Pacific EOF analysis",
         subtitle  = paste0("Domain: ", EOF_LON_MIN, "°E–", EOF_LON_MAX-360, "°W ",
                            "(", EOF_LON_MIN, "–", EOF_LON_MAX, "° in 0-360), ",
                            EOF_LAT_MIN, "°–", EOF_LAT_MAX, "°N  |  ",
                            "Anomalies relative to WMO ", CLIM_START, "–", CLIM_END,
                            "  |  Red dotted: 80% cumulative variance"),
         x = "PC number") +
    theme_ts(10)
  
  save_fig(p_scree, "FigS1_scree_diagrams", w = 12, h = 5)
  
} # end scree guard

# ── 10-B: EOF SPATIAL PATTERN MAPS ──────────────────────────────────────────
cat("\n  10-B: EOF spatial pattern maps\n")

# Shared map infrastructure (reuse borders from w9 if still in environment)
if (!exists("borders")) {
  borders <- ne_countries(scale = "medium", returnclass = "sf")
}

make_eof_panel <- function(eof_res, pc_k, fld_label, var_pct_k) {
  r_eof <- eof_to_raster(eof_res, pc_k)
  
  # Wrap lon to -180:180 for ggplot
  df_eof <- as.data.frame(r_eof, xy = TRUE, na.rm = TRUE)
  names(df_eof)[3] <- "loading"
  df_eof$x <- wrap_lon(df_eof$x)
  
  # Symmetric colour limits (98th pctile of absolute value)
  lim <- quantile(abs(df_eof$loading), 0.98, na.rm = TRUE)
  
  ggplot2::ggplot() +
    ggplot2::geom_raster(data = df_eof,
                         aes(x = x, y = y, fill = loading)) +
    ggplot2::scale_fill_distiller(
      palette   = PAL_EOF,
      direction = -1,
      limits    = c(-lim, lim),
      oob       = scales::squish,
      name      = "Loading",
      guide     = guide_colorbar(barwidth = 6, barheight = 0.4,
                                 title.position = "top")) +
    ggplot2::geom_sf(data = borders, fill = NA, colour = "grey30",
                     linewidth = 0.2, inherit.aes = FALSE) +
    ggplot2::coord_sf(
      xlim   = c(140 - 360, -120),   # 140°E = -220 wraps to +140; use -220:-120 in WGS84
      ylim   = c(EOF_LAT_MIN, EOF_LAT_MAX),
      expand = FALSE, crs = sf::st_crs(4326)) +
    ggplot2::labs(
      title    = sprintf("(%s) %s EOF%d",
                         letters[(which(names(fields) == sub("_.*","",fld_label)) - 1) *
                                   2 + pc_k],
                         fld_label, pc_k),
      subtitle = sprintf("%.1f%% variance", var_pct_k)) +
    ggplot2::theme_void(base_size = 8.5) +
    ggplot2::theme(
      plot.title    = element_text(size = 8.5, face = "bold", hjust = 0,
                                   margin = margin(b = 2)),
      plot.subtitle = element_text(size = 7, colour = "grey35"),
      legend.position  = "bottom",
      legend.key.size  = unit(0.35, "cm"),
      legend.text      = element_text(size = 6.5),
      legend.title     = element_text(size = 7, face = "bold"),
      panel.border     = element_rect(colour = "grey50", fill = NA,
                                      linewidth = 0.3))
}

panels_eof <- list()
for (fld in names(fields)) {
  key_all <- sprintf("%s_all", fld)
  if (!key_all %in% names(eof_results)) next
  res <- eof_results[[key_all]]
  for (k in 1:2) {
    if (k > res$n_pcs) next
    p_k <- tryCatch(
      make_eof_panel(res, k, fld, res$var_pct[k]),
      error = function(e) { cat(sprintf("  ✗ EOF map %s PC%d: %s\n", fld, k, e$message)); NULL })
    if (!is.null(p_k)) panels_eof[[sprintf("%s_PC%d", fld, k)]] <- p_k
  }
}

if (length(panels_eof) >= 2) {
  fig_eof <- wrap_plots(panels_eof, ncol = 2, guides = "collect") &
    theme(legend.position = "bottom")
  
  fig_eof <- fig_eof +
    plot_annotation(
      title    = "EOF spatial patterns — North Pacific circulation and SST",
      subtitle = sprintf(
        "Domain: %d°E–%d°W, %d°–%d°N | Anomalies vs WMO %d–%d | Cosine-lat weighted | All months",
        EOF_LON_MIN, 360 - EOF_LON_MAX,
        EOF_LAT_MIN, EOF_LAT_MAX, CLIM_START, CLIM_END),
      theme = theme(plot.title    = element_text(size = 10, face = "bold",
                                                 hjust = 0.5),
                    plot.subtitle = element_text(size = 7.5, colour = "grey35",
                                                 hjust = 0.5)))
  
  save_fig(fig_eof, "FigS2_eof_spatial_patterns", w = 10, h = 3 * ceiling(length(panels_eof) / 2))
}

# ── 10-C: PC TIME SERIES WITH SPEI OVERLAY ───────────────────────────────────
cat("\n  10-C: PC time series with SPEI overlay\n")

pc_ts_panels <- list()

# Use SPEI-3 as the primary overlay series (fallback to first available SPEI)
spei_overlay_nm <- if ("SPEI3" %in% names(index_list)) "SPEI3" else
  names(index_list)[startsWith(names(index_list), "SPEI")][1]
spei_overlay    <- if (!is.na(spei_overlay_nm) && spei_overlay_nm %in% names(index_list)) {
  index_list[[spei_overlay_nm]]
} else NULL

for (fld in names(fields)) {
  key_all <- sprintf("%s_all", fld)
  if (!key_all %in% names(eof_results)) next
  res      <- eof_results[[key_all]]
  pc1_df   <- data.frame(date = res$dates_sub, pc1 = res$scores[, 1])
  var1_pct <- round(res$var_pct[1], 1)
  
  # Drought focus period shading rectangle
  shade_df <- data.frame(
    xmin = as.Date(sprintf("%d-01-01", DROUGHT_FOCUS_START)),
    xmax = as.Date(sprintf("%d-12-31", DROUGHT_FOCUS_END)),
    ymin = -Inf, ymax = Inf)
  
  p_pc <- ggplot() +
    # 2022-2025 shading
    geom_rect(data = shade_df,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = COL_FOCUS, alpha = 0.12, inherit.aes = FALSE) +
    # Zero line
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.4) +
    # PC1 line
    geom_line(data = pc1_df, aes(x = date, y = pc1,
                                 colour = "PC1"), linewidth = 0.7) +
    scale_y_continuous(
      name = sprintf("%s PC1 (normalised)", fld),
      sec.axis = sec_axis(~ .,
                          name = sprintf("%s (standardised)", spei_overlay_nm)))
  
  if (!is.null(spei_overlay)) {
    mg_sp <- inner_join(pc1_df, spei_overlay %>% select(date, spei = value),
                        by = "date")
    p_pc <- p_pc +
      geom_line(data = mg_sp, aes(x = date, y = spei, colour = "SPEI"),
                linewidth = 0.7, linetype = "dashed", alpha = 0.85)
  }
  
  p_pc <- p_pc +
    scale_colour_manual(
      values = c("PC1" = COL_PC, "SPEI" = COL_SPEI),
      labels = c("PC1" = sprintf("%s PC1 (%.1f%% var.)", fld, var1_pct),
                 "SPEI" = spei_overlay_nm),
      name   = NULL) +
    scale_x_date(date_breaks = "10 years", date_labels = "%Y") +
    # Annotate drought period
    annotate("text", x = as.Date(sprintf("%d-07-01", DROUGHT_FOCUS_START)),
             y = Inf, label = sprintf("%d–%d\nDrought",
                                      DROUGHT_FOCUS_START, DROUGHT_FOCUS_END),
             vjust = 1.8, hjust = 0.5, size = 2.8,
             colour = COL_FOCUS, fontface = "bold") +
    labs(
      title    = sprintf("(%s)  %s PC1 and %s — %d–%d",
                         letters[which(names(fields) == fld)],
                         fld, spei_overlay_nm, START_YEAR, END_YEAR),
      subtitle = sprintf(
        "PC1 explains %.1f%% of North Pacific %s variance | %d–%d highlighted",
        var1_pct, fld, DROUGHT_FOCUS_START, DROUGHT_FOCUS_END),
      x = NULL) +
    theme_ts(9.5)
  
  pc_ts_panels[[fld]] <- p_pc
}

if (length(pc_ts_panels) >= 1) {
  fig_pc_ts <- wrap_plots(pc_ts_panels, ncol = 1) +
    plot_annotation(
      title    = sprintf(
        "North Pacific EOF PC1 time series and %s — Nechako Basin drought context",
        spei_overlay_nm),
      subtitle = paste0(
        "EOF domain: ", EOF_LON_MIN, "°E–", 360 - EOF_LON_MAX, "°W, ",
        EOF_LAT_MIN, "°–", EOF_LAT_MAX, "°N  |  ",
        "Shading: ", DROUGHT_FOCUS_START, "–", DROUGHT_FOCUS_END, " drought period  |  ",
        "Dashed: ", spei_overlay_nm, " (standardised, right axis)"),
      theme = theme(
        plot.title    = element_text(size = 9.5, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 7.5, colour = "grey35", hjust = 0.5,
                                     margin = margin(b = 4))))
  
  save_fig(fig_pc_ts, "Fig8_pc_timeseries_SPEI",
           w = FIG_W, h = 3.5 * length(pc_ts_panels))
}

# ── 10-D: CORRELATION HEATMAP ────────────────────────────────────────────────
cat("\n  10-D: Correlation heatmap\n")

if (!is.null(corr_df)) {
  
  # ── Filter: SPEI-PM only (drop SPI), lag = 0, PC1-PC3 only ─────────────────
  SPEI_LABELS <- c("SPEI1", "SPEI2", "SPEI3")
  MS_PCS      <- c("PC1", "PC2", "PC3")
  
  heatmap_df <- corr_df %>%
    dplyr::filter(
      lag_mo   == 0,
      index_nm %in% SPEI_LABELS,
      pc       %in% paste0(field, "-", MS_PCS[match(
        sub(".*-", "", pc), MS_PCS)])
    ) %>%
    # pc column already contains e.g. "Z500-PC1" — use directly as label
    dplyr::mutate(
      pc_label  = pc,
      sig_label = ifelse(p_effn < 0.01, "**",
                         ifelse(p_effn < ALPHA_CORR, "*", "")),
      r_label   = sprintf("%.2f%s", r, sig_label),
      pc_label  = factor(pc_label,
                         levels = c("Z500-PC1", "Z500-PC2", "Z500-PC3",
                                    "SLP-PC1",  "SLP-PC2",  "SLP-PC3",
                                    "SST-PC1",  "SST-PC2",  "SST-PC3")),
      spei_label = dplyr::recode(index_nm,
                                 SPEI1 = "SPEI-1",
                                 SPEI2 = "SPEI-2",
                                 SPEI3 = "SPEI-3"),
      spei_label = factor(spei_label, levels = c("SPEI-1", "SPEI-2", "SPEI-3"))
    ) %>%
    dplyr::filter(!is.na(pc_label))
  
  # Vertical separator positions (between Z500/SLP and SLP/SST groups)
  vline_pos <- c(3.5, 6.5)
  
  p_heat <- ggplot2::ggplot(heatmap_df,
                            ggplot2::aes(x = pc_label, y = spei_label, fill = r)) +
    ggplot2::geom_vline(xintercept = vline_pos, colour = "grey55",
                        linewidth = 0.9, linetype = "solid") +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = r_label),
                       size = 3.8, colour = "grey10", fontface = "bold") +
    # Field group labels above tiles
    ggplot2::annotate("text", x = 2, y = 3.7, label = "Z500",
                      size = 3.8, fontface = "bold", colour = "grey20") +
    ggplot2::annotate("text", x = 5, y = 3.7, label = "SLP",
                      size = 3.8, fontface = "bold", colour = "grey20") +
    ggplot2::annotate("text", x = 8, y = 3.7, label = "SST",
                      size = 3.8, fontface = "bold", colour = "grey20") +
    ggplot2::scale_fill_distiller(
      palette   = "RdBu",
      direction = 1,
      limits    = c(-1, 1),
      oob       = scales::squish,
      name      = "Pearson r",
      guide     = ggplot2::guide_colorbar(
        barwidth = 10, barheight = 0.5,
        title.position = "top", title.hjust = 0.5)) +
    ggplot2::scale_x_discrete(name = NULL) +
    ggplot2::scale_y_discrete(
      name   = NULL,
      expand = ggplot2::expansion(add = c(0.5, 0.9))) +
    ggplot2::labs(
      title    = "North Pacific circulation PCs vs basin-average SPEI (PM) — Nechako Basin",
      subtitle = paste0(
        "Lag = 0 months  |  PC1\u2013PC3 for Z500, SLP, SST  |  ",
        "** p < 0.01, * p < 0.05 (effective sample size)  |  ",
        "Anomalies vs WMO ", CLIM_START, "\u2013", CLIM_END)) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(size = 11, face = "bold", hjust = 0,
                                              margin = ggplot2::margin(b = 4)),
      plot.subtitle   = ggplot2::element_text(size = 8.5, colour = "grey35", hjust = 0,
                                              margin = ggplot2::margin(b = 8)),
      axis.text.x     = ggplot2::element_text(size = 10, angle = 35, hjust = 1),
      axis.text.y     = ggplot2::element_text(size = 11, face = "bold"),
      panel.grid      = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.title    = ggplot2::element_text(size = 9, face = "bold"),
      legend.text     = ggplot2::element_text(size = 9),
      plot.margin     = ggplot2::margin(t = 10, r = 10, b = 8, l = 8)
    )
  
  FIG9_W <- 9.5
  FIG9_H <- 4.2
  
  tryCatch({
    ggplot2::ggsave(file.path(OUT_DIR, "Fig9_correlation_heatmap.pdf"),
                    p_heat, width = FIG9_W, height = FIG9_H,
                    units = "in", device = "pdf")
    cat("  \u2713 Fig9_correlation_heatmap.pdf\n")
  }, error = function(e) cat(sprintf("  \u2717 Fig9_correlation_heatmap.pdf: %s\n", e$message)))
  
  tryCatch({
    ggplot2::ggsave(file.path(OUT_DIR, "Fig9_correlation_heatmap.png"),
                    p_heat, width = FIG9_W, height = FIG9_H,
                    units = "in", dpi = 600L, device = "png")
    cat("  \u2713 Fig9_correlation_heatmap.png (600 DPI)\n")
  }, error = function(e) cat(sprintf("  \u2717 Fig9_correlation_heatmap.png: %s\n", e$message)))
}

# ============================================================================
#  SECTION 11: MANUSCRIPT TEXT BLANKS
# ============================================================================
cat("\n============================================================\n")
cat(" SECTION 11: Manuscript text blanks CSV\n")
cat("============================================================\n")

# Add regression adj-R² values already in ms_blanks; add combined-field values
# Multi-field regression: collect adj-R² across fields
all_adjr2 <- unlist(ms_blanks[grep("^reg_adjR2_", names(ms_blanks))])
if (length(all_adjr2)) {
  ms_blanks[["reg_adjR2_max"]] <- round(max(all_adjr2, na.rm = TRUE), 1)
  ms_blanks[["reg_adjR2_min"]] <- round(min(all_adjr2, na.rm = TRUE), 1)
}

# [FIX] Guard against an empty ms_blanks list.
# When eof_results is empty (Section 4 failed for all fields), Sections 6–8
# skip entirely and ms_blanks stays as list().  The original code then crashed
# because data.frame() received columns of length 0, 0, and 1 (NA_character_).
# We write an empty CSV with the correct column structure and exit gracefully.
if (length(ms_blanks) == 0) {
  warning(
    "ms_blanks is empty — no manuscript fill-in values were produced.\n",
    "  This means all EOF combinations in Section 4 failed.\n",
    "  Writing zero-row ms2_eof_blanks.csv and skipping Section 11 output.",
    call. = FALSE
  )
  blanks_df <- data.frame(
    placeholder = character(0),
    value       = character(0),
    description = character(0)
  )
  write_csv(blanks_df, file.path(OUT_DIR, "ms2_eof_blanks.csv"))
  cat("  ms2_eof_blanks.csv saved (0 rows — no EOF results).\n")
} else {
  # [FIX] Use rep(NA_character_, length(ms_blanks)) instead of NA_character_
  # so all three columns have the same length and data.frame() does not crash.
  blanks_df <- data.frame(
    placeholder = names(ms_blanks),
    value       = unlist(ms_blanks),
    description = rep(NA_character_, length(ms_blanks)),
    row.names   = NULL
  )
  
  # Human-readable descriptions for each placeholder
  desc_map <- c(
    "varexp_Z500_PC1"              = "% variance explained by Z500 EOF1",
    "varexp_SLP_PC1"               = "% variance explained by SLP EOF1",
    "varexp_SST_PC1"               = "% variance explained by SST EOF1",
    "event_mean_Z500_PC1"          = "Mean Z500-PC1 during 2022-2025 (normalised σ)",
    "event_mean_SLP_PC1"           = "Mean SLP-PC1 during 2022-2025 (normalised σ)",
    "event_mean_SST_PC1"           = "Mean SST-PC1 during 2022-2025 (normalised σ)",
    "event_zscore_Z500_PC1"        = "Z-score of 2022-2025 Z500-PC1 mean vs all windows",
    "event_rank_Z500_PC1"          = "Rank of 2022-2025 Z500-PC1 mean (1=highest)",
    "event_n_windows_Z500_PC1"     = "Number of overlapping event-length windows in record",
    "event_return_yr_Z500_PC1"     = "Estimated return period of 2022-2025 Z500-PC1 (years)",
    "reg_adjR2_max"                = "Max adjusted R² across fields (SPEI ~ PC1+PC2+PC3, %)",
    "reg_adjR2_min"                = "Min adjusted R² across fields (SPEI ~ PC1+PC2+PC3, %)"
  )
  blanks_df$description <- desc_map[blanks_df$placeholder]
  
  write_csv(blanks_df, file.path(OUT_DIR, "ms2_eof_blanks.csv"))
  cat("  ms2_eof_blanks.csv saved.\n")
  
  cat("\n  KEY MANUSCRIPT FILL-IN VALUES:\n")
  cat("  ┌──────────────────────────────────────────────────────────────┐\n")
  for (nm in names(ms_blanks)) {
    cat(sprintf("  │  %-42s = %-12s │\n",
                nm, as.character(ms_blanks[[nm]])))
  }
  cat("  └──────────────────────────────────────────────────────────────┘\n")
}

# ============================================================================
#  DONE
# ============================================================================
cat("\n============================================================\n")
cat(" w10c COMPLETE\n")
cat("============================================================\n")
cat(sprintf("  Output directory: %s\n", OUT_DIR))
cat("  Files written:\n")
cat("    w10c_eof_state.rds            — full EOF results object\n")
cat("    eof_pc_timeseries.csv         — monthly PC time series\n")
cat("    Table4_PC_correlations.csv    — PC vs SPEI lag correlations\n")
cat("    Table4_PC_vs_teleconnections.csv — PC vs PDO/PNA/Niño-3.4/AO\n")
cat("    Table5_regression.csv         — SPEI ~ PC1+PC2+PC3 (Newey-West)\n")
cat("    event_quantification_2022_2025.csv\n")
cat("    pc_trend_analysis.csv\n")
cat("    ms2_eof_blanks.csv            — plug these into §4.7.1, §4.7.2\n")
cat("    FigS1_scree_diagrams.pdf/png\n")
cat("    FigS2_eof_spatial_patterns.pdf/png\n")
cat("    Fig8_pc_timeseries_SPEI.pdf/png\n")
cat("    Fig9_correlation_heatmap.pdf/png\n")
cat("\n")
cat("  SCRIPT CHAIN:\n")
cat("    w9 → w10a → w10b → w10c  (all independent after w9)\n")
cat("============================================================\n")