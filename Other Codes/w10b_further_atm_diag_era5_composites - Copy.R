# ==============================================================================
# w10a_era5_composites.R
# ERA5 DATA PREPARATION + B1 HISTORICAL COMPOSITES + B3 BOOTSTRAP SIGNIFICANCE
#
# Split from w10_further_atmospheric_diagnostics.R (was one monolithic file).
# B2 (climate index skill) is now in w10b_index_skill.R — it does not need ERA5.
#
# Modules:
#   B1  Historical composite analysis (Z500, SLP, SST)
#   B3  Bootstrap fingerprint significance test
#
# Outputs:
#   composite_results/   (B1 rasters, spatial correlations, severity ranking)
#   bootstrap_results/   (B3 p-value rasters, stippled maps, analogue list)
#
# Run AFTER: w9_atmospheric_diagnostics.R  (anomaly NetCDFs are reused if present)
# Run BEFORE: w10b_index_skill.R (independent — can run in parallel)
# ==============================================================================

library(terra)
library(lubridate)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(parallel)

setwd("D:/Nechako_Drought/Nechako/")

for (d in c("composite_results", "bootstrap_results"))
  if (!dir.exists(d)) dir.create(d)

DATA_DIR <- "monthly_data_direct"   # same as w9

cat("\n========== w10a: LOADING COMMON DATASETS ==========\n")

# ------------------------------------------------------------------------------
# DATES  (from precipitation file)
# ------------------------------------------------------------------------------
precip       <- rast(file.path(DATA_DIR, "total_precipitation_monthly.nc"))
dates_precip <- as.Date(time(precip))
if (is.null(dates_precip) || all(is.na(dates_precip)))
  dates_precip <- seq(as.Date("1950-01-01"), by = "month", length.out = nlyr(precip))
month_nums <- as.integer(format(dates_precip, "%m"))
year_nums  <- as.integer(format(dates_precip, "%Y"))
cat(sprintf("Date range: %s to %s (%d months)\n",
            min(dates_precip), max(dates_precip), length(dates_precip)))

# ------------------------------------------------------------------------------
# ERA5 ANOMALIES
# Reuse pre-computed anomaly NetCDFs from the w9 pipeline when available;
# fall back to computing from raw ERA5 if they are absent.
#
# PERF: compute_anomaly() now subtracts the climatology mean from all layers
#       of a given calendar month in one vectorised raster operation instead of
#       looping over individual layers (was O(n_years) raster writes per month).
# ------------------------------------------------------------------------------

compute_anomaly <- function(r, clim_start = 1991, clim_end = 2020) {
  dates_r  <- as.Date(time(r))
  month_r  <- as.integer(format(dates_r, "%m"))
  year_r   <- as.integer(format(dates_r, "%Y"))
  clim_idx <- which(year_r >= clim_start & year_r <= clim_end)
  
  # Build a list of (layer indices, anomaly raster) pairs — one per calendar month
  anom_parts <- vector("list", 12)
  for (m in 1:12) {
    m_all  <- which(month_r == m)
    m_clim <- intersect(m_all, clim_idx)
    if (!length(m_clim)) next
    clim_m          <- app(r[[m_clim]], mean, na.rm = TRUE)   # 1-layer climatology
    anom_parts[[m]] <- list(idx = m_all, anom = r[[m_all]] - clim_m)
    # terra broadcasts the single-layer clim_m across all layers of r[[m_all]]
  }
  
  # Reassemble in original chronological order
  valid       <- !vapply(anom_parts, is.null, logical(1))
  all_idx     <- unlist(lapply(anom_parts[valid], `[[`, "idx"))
  all_anom    <- do.call(c, lapply(anom_parts[valid], `[[`, "anom"))
  result      <- all_anom[[order(all_idx)]]
  terra::time(result) <- dates_r
  result
}

load_or_compute_anomaly <- function(raw_nc, anomaly_nc,
                                    unit_divisor = 1,
                                    unit_offset  = 0,
                                    clim_start   = 1991,
                                    clim_end     = 2020) {
  if (file.exists(anomaly_nc)) {
    cat(sprintf("  Reusing pre-computed anomaly: %s\n", basename(anomaly_nc)))
    return(rast(anomaly_nc))
  }
  cat(sprintf("  Computing anomaly from %s...\n", basename(raw_nc)))
  r <- rast(raw_nc)
  if (unit_divisor != 1) r <- r / unit_divisor
  if (unit_offset  != 0) r <- r + unit_offset
  compute_anomaly(r, clim_start, clim_end)
}

cat("\nLoading / computing ERA5 anomalies...\n")

# w9 saves anomalies here; reuse them to avoid redundant computation
z500_anom <- load_or_compute_anomaly(
  raw_nc      = file.path(DATA_DIR, "geopotential_500hPa_monthly.nc"),
  anomaly_nc  = file.path(DATA_DIR, "z500_monthly_anomaly_clim1991_2020.nc"),
  unit_divisor = 9.80665)

slp_anom <- load_or_compute_anomaly(
  raw_nc      = file.path(DATA_DIR, "mean_sea_level_pressure_monthly.nc"),
  anomaly_nc  = file.path(DATA_DIR, "slp_monthly_anomaly_clim1991_2020.nc"),
  unit_divisor = 100)

sst_anom <- load_or_compute_anomaly(
  raw_nc      = file.path(DATA_DIR, "sea_surface_temperature_monthly.nc"),
  anomaly_nc  = file.path(DATA_DIR, "sst_monthly_anomaly_clim1991_2020.nc"),
  unit_offset  = -273.15)

dates_era5 <- as.Date(time(z500_anom))
cat(sprintf("ERA5 period: %s to %s\n", min(dates_era5), max(dates_era5)))

# ------------------------------------------------------------------------------
# T0 DROUGHT MONTHS  (produced by w9)
# ------------------------------------------------------------------------------
cat("\nLoading T0 drought months...\n")
t0_file <- "T0_drought_months_all_indices.csv"
if (!file.exists(t0_file)) stop("T0_drought_months_all_indices.csv not found.")

t0 <- read_csv(t0_file, show_col_types = FALSE)
t0_primary <- t0 %>%
  filter(index %in% c("SPEI6","SPEI12","SPI6","SPI12")) %>%
  mutate(date_ym = as.Date(paste0(date, "-01")))

drought_consensus <- t0_primary %>%
  group_by(date) %>%
  summarise(n_indices  = n_distinct(index),
            mean_value = mean(index_value, na.rm = TRUE),
            year       = first(year),
            month_num  = first(month_num),
            season     = first(season), .groups = "drop") %>%
  filter(n_indices >= 2) %>%
  mutate(date_ym          = as.Date(paste0(date, "-01")),
         is_drought_period = year >= 2022 & year <= 2025)

cat(sprintf("Consensus drought months: %d total (%d in 2022-2025)\n",
            nrow(drought_consensus), sum(drought_consensus$is_drought_period)))

# ==============================================================================
# MODULE B1: HISTORICAL COMPOSITE ANALYSIS
# ==============================================================================
cat("\n========== MODULE B1: HISTORICAL COMPOSITES ==========\n")
out_b1 <- "composite_results"

all_drought_dates  <- drought_consensus$date_ym
drought2025_dates  <- drought_consensus$date_ym[ drought_consensus$is_drought_period]
historical_dates   <- drought_consensus$date_ym[!drought_consensus$is_drought_period]

# PDO proxy from SST box
sst_pdo_box <- crop(sst_anom, ext(-160, -130, 40, 55))
pdo_proxy   <- global(sst_pdo_box, fun = "mean", na.rm = TRUE)[, 1]
pdo_df      <- data.frame(date_ym = as.Date(time(sst_anom)), pdo_proxy = pdo_proxy)

drought_consensus <- drought_consensus %>%
  left_join(pdo_df, by = "date_ym") %>%
  mutate(pdo_phase = ifelse(pdo_proxy > 0, "PDO+", "PDO-"))

pdo_pos_dates <- drought_consensus$date_ym[drought_consensus$pdo_phase == "PDO+" &
                                             !drought_consensus$is_drought_period]
pdo_neg_dates <- drought_consensus$date_ym[drought_consensus$pdo_phase == "PDO-" &
                                             !drought_consensus$is_drought_period]
cat(sprintf("Historical drought months: PDO+ = %d, PDO- = %d\n",
            length(pdo_pos_dates), length(pdo_neg_dates)))

composite_months <- function(anom_rast, target_dates, label) {
  dates_r <- as.Date(time(anom_rast))
  idx <- which(dates_r %in% target_dates)
  if (!length(idx)) return(NULL)
  comp <- app(anom_rast[[idx]], mean, na.rm = TRUE)   # app() is faster than mean()
  names(comp) <- label
  comp
}

z500_all_hist  <- composite_months(z500_anom, historical_dates,  "Z500_all_historical")
z500_pdo_pos   <- composite_months(z500_anom, pdo_pos_dates,     "Z500_PDO+")
z500_pdo_neg   <- composite_months(z500_anom, pdo_neg_dates,     "Z500_PDO-")
z500_2022_2025 <- composite_months(z500_anom, drought2025_dates, "Z500_2022-2025")
slp_all_hist   <- composite_months(slp_anom,  historical_dates,  "SLP_all_historical")
slp_pdo_pos    <- composite_months(slp_anom,  pdo_pos_dates,     "SLP_PDO+")
slp_pdo_neg    <- composite_months(slp_anom,  pdo_neg_dates,     "SLP_PDO-")
slp_2022_2025  <- composite_months(slp_anom,  drought2025_dates, "SLP_2022-2025")
sst_all_hist   <- composite_months(sst_anom,  historical_dates,  "SST_all_historical")
sst_pdo_pos    <- composite_months(sst_anom,  pdo_pos_dates,     "SST_PDO+")
sst_pdo_neg    <- composite_months(sst_anom,  pdo_neg_dates,     "SST_PDO-")
sst_2022_2025  <- composite_months(sst_anom,  drought2025_dates, "SST_2022-2025")

spatial_cor <- function(r1, r2) {
  v1 <- as.vector(values(r1)); v2 <- as.vector(values(r2))
  ok <- !is.na(v1) & !is.na(v2)
  if (sum(ok) < 10) return(NA)
  cor(v1[ok], v2[ok])
}

r_z500_all    <- spatial_cor(z500_2022_2025, z500_all_hist)
r_z500_pdopos <- spatial_cor(z500_2022_2025, z500_pdo_pos)
r_z500_pdoneg <- spatial_cor(z500_2022_2025, z500_pdo_neg)
r_slp_all     <- spatial_cor(slp_2022_2025,  slp_all_hist)
r_slp_pdopos  <- spatial_cor(slp_2022_2025,  slp_pdo_pos)

cat("\nSpatial correlations (2022-2025 vs composites):\n")
cat(sprintf("  Z500 vs All Hist: %.3f\n", r_z500_all))
cat(sprintf("  Z500 vs PDO+    : %.3f\n", r_z500_pdopos))
cat(sprintf("  Z500 vs PDO-    : %.3f\n", r_z500_pdoneg))
cat(sprintf("  SLP  vs All Hist: %.3f\n", r_slp_all))
cat(sprintf("  SLP  vs PDO+    : %.3f\n", r_slp_pdopos))

severity_by_year <- drought_consensus %>%
  group_by(year) %>%
  summarise(n_drought_months = n(),
            mean_severity    = mean(abs(mean_value), na.rm = TRUE),
            composite_score  = n_drought_months * mean_severity, .groups = "drop") %>%
  arrange(desc(composite_score)) %>%
  mutate(rank = row_number())

write.csv(severity_by_year, file.path(out_b1,"drought_severity_ranking.csv"), row.names=FALSE)

writeRaster(z500_all_hist,  file.path(out_b1,"z500_composite_all_hist.tif"),  overwrite=TRUE)
writeRaster(z500_pdo_pos,   file.path(out_b1,"z500_composite_pdopos.tif"),    overwrite=TRUE)
writeRaster(z500_pdo_neg,   file.path(out_b1,"z500_composite_pdoneg.tif"),    overwrite=TRUE)
writeRaster(z500_2022_2025, file.path(out_b1,"z500_composite_2022_2025.tif"), overwrite=TRUE)
writeRaster(slp_all_hist,   file.path(out_b1,"slp_composite_all_hist.tif"),   overwrite=TRUE)
writeRaster(slp_pdo_pos,    file.path(out_b1,"slp_composite_pdopos.tif"),     overwrite=TRUE)
writeRaster(sst_all_hist,   file.path(out_b1,"sst_composite_all_hist.tif"),   overwrite=TRUE)
writeRaster(sst_pdo_pos,    file.path(out_b1,"sst_composite_pdopos.tif"),     overwrite=TRUE)
writeRaster(sst_2022_2025,  file.path(out_b1,"sst_composite_2022_2025.tif"),  overwrite=TRUE)

write.csv(data.frame(
  variable   = c("Z500","Z500","Z500","SLP","SLP"),
  comparison = c("All Historical","PDO+","PDO-","All Historical","PDO+"),
  spatial_r  = c(r_z500_all, r_z500_pdopos, r_z500_pdoneg, r_slp_all, r_slp_pdopos)
), file.path(out_b1,"composite_spatial_correlations.csv"), row.names=FALSE)

cat("✓ B1 composites saved to:", out_b1, "\n")

# ==============================================================================
# MODULE B3: BOOTSTRAP FINGERPRINT SIGNIFICANCE
# PERF: inner bootstrap loop parallelised across available cores using a PSOCK
#       cluster (Windows-safe).  Each worker accumulates its own exceedance
#       count matrix; results are summed in the main process.
# ==============================================================================
cat("\n========== MODULE B3: BOOTSTRAP SIGNIFICANCE ==========\n")
out_b3     <- "bootstrap_results"
N_BOOT     <- 5000L
BLOCK_SIZE <- 6L
N_DROUGHT  <- 48L
set.seed(42L)

bg_idx      <- which(year_nums >= 1950 & year_nums <= 2021)
drought_idx <- which(year_nums >= 2022 & year_nums <= 2025)

if (length(bg_idx) < 100 || length(drought_idx) != N_DROUGHT) {
  cat("Insufficient data for bootstrap – skipping B3.\n")
} else {
  
  # Extract background matrices once (pixels × time) — main memory cost
  z500_bg_mat <- values(z500_anom[[bg_idx]])
  slp_bg_mat  <- values(slp_anom[[bg_idx]])
  sst_bg_mat  <- values(sst_anom[[bg_idx]])
  
  obs_z500     <- app(z500_anom[[drought_idx]], mean, na.rm = TRUE)
  obs_slp      <- app(slp_anom[[drought_idx]],  mean, na.rm = TRUE)
  obs_sst      <- app(sst_anom[[drought_idx]],  mean, na.rm = TRUE)
  obs_z500_vec <- as.vector(values(obs_z500))
  obs_slp_vec  <- as.vector(values(obs_slp))
  obs_sst_vec  <- as.vector(values(obs_sst))
  
  n_bg           <- length(bg_idx)
  n_pix          <- nrow(z500_bg_mat)
  block_starts   <- seq_len(n_bg - BLOCK_SIZE + 1L)
  n_blocks_needed <- ceiling(N_DROUGHT / BLOCK_SIZE)
  
  # ── Parallel block bootstrap ──────────────────────────────────────────────
  n_cores  <- max(1L, min(detectCores() - 1L, 4L))   # cap at 4 to limit RAM
  chunk_sz <- ceiling(N_BOOT / n_cores)
  chunks   <- split(seq_len(N_BOOT), ceiling(seq_len(N_BOOT) / chunk_sz))
  
  cat(sprintf("  Running %d bootstrap reps across %d core(s)...\n", N_BOOT, n_cores))
  
  run_chunk <- function(boot_ids) {
    cnt_z <- cnt_s <- cnt_sst <- numeric(n_pix)
    set.seed(boot_ids[1])  # deterministic per chunk
    for (b in boot_ids) {
      starts    <- sample(block_starts, n_blocks_needed, replace = TRUE)
      drawn_idx <- pmin(unlist(lapply(starts, function(s) s:(s + BLOCK_SIZE - 1L))), n_bg)[seq_len(N_DROUGHT)]
      cnt_z   <- cnt_z   + (abs(rowMeans(z500_bg_mat[, drawn_idx, drop = FALSE], na.rm = TRUE)) >= abs(obs_z500_vec))
      cnt_s   <- cnt_s   + (abs(rowMeans(slp_bg_mat[,  drawn_idx, drop = FALSE], na.rm = TRUE)) >= abs(obs_slp_vec))
      cnt_sst <- cnt_sst + (abs(rowMeans(sst_bg_mat[,  drawn_idx, drop = FALSE], na.rm = TRUE)) >= abs(obs_sst_vec))
    }
    list(z = cnt_z, s = cnt_s, sst = cnt_sst)
  }
  
  if (n_cores > 1L) {
    cl <- makeCluster(n_cores, type = "PSOCK")
    clusterExport(cl, c("z500_bg_mat","slp_bg_mat","sst_bg_mat",
                        "block_starts","n_blocks_needed","N_DROUGHT","n_bg",
                        "BLOCK_SIZE","obs_z500_vec","obs_slp_vec","obs_sst_vec","n_pix"),
                  envir = environment())
    results <- parLapply(cl, chunks, run_chunk)
    stopCluster(cl)
  } else {
    results <- lapply(chunks, run_chunk)
  }
  
  count_exceed_z500 <- Reduce(`+`, lapply(results, `[[`, "z"))
  count_exceed_slp  <- Reduce(`+`, lapply(results, `[[`, "s"))
  count_exceed_sst  <- Reduce(`+`, lapply(results, `[[`, "sst"))
  
  pval_z500 <- count_exceed_z500 / N_BOOT
  pval_slp  <- count_exceed_slp  / N_BOOT
  pval_sst  <- count_exceed_sst  / N_BOOT
  
  fdr_correct <- function(pvals, alpha = 0.05) {
    valid <- !is.na(pvals)
    adj   <- rep(NA_real_, length(pvals))
    adj[valid] <- p.adjust(pvals[valid], method = "BH")
    adj < alpha
  }
  sig_z500 <- fdr_correct(pval_z500)
  sig_slp  <- fdr_correct(pval_slp)
  sig_sst  <- fdr_correct(pval_sst)
  
  cat(sprintf("Fraction significant (FDR q<0.05): Z500=%.1f%%, SLP=%.1f%%, SST=%.1f%%\n",
              100*mean(sig_z500,na.rm=TRUE), 100*mean(sig_slp,na.rm=TRUE),
              100*mean(sig_sst,na.rm=TRUE)))
  
  make_pval_rast <- function(pv, tmpl) { r <- tmpl; values(r) <- pv; r }
  
  writeRaster(make_pval_rast(pval_z500,           obs_z500),
              file.path(out_b3,"pval_z500_bootstrap.tif"), overwrite=TRUE)
  writeRaster(make_pval_rast(pval_slp,            obs_slp),
              file.path(out_b3,"pval_slp_bootstrap.tif"),  overwrite=TRUE)
  writeRaster(make_pval_rast(pval_sst,            obs_sst),
              file.path(out_b3,"pval_sst_bootstrap.tif"),  overwrite=TRUE)
  writeRaster(make_pval_rast(as.numeric(sig_z500), obs_z500),
              file.path(out_b3,"sig_z500_fdr005.tif"), overwrite=TRUE)
  writeRaster(make_pval_rast(as.numeric(sig_slp),  obs_slp),
              file.path(out_b3,"sig_slp_fdr005.tif"),  overwrite=TRUE)
  writeRaster(make_pval_rast(as.numeric(sig_sst),  obs_sst),
              file.path(out_b3,"sig_sst_fdr005.tif"),  overwrite=TRUE)
  
  write.csv(data.frame(
    variable            = c("Z500","SLP","SST"),
    pct_significant_raw = c(100*mean(pval_z500<0.05,na.rm=TRUE),
                            100*mean(pval_slp<0.05, na.rm=TRUE),
                            100*mean(pval_sst<0.05, na.rm=TRUE)),
    pct_significant_fdr = c(100*mean(sig_z500,na.rm=TRUE),
                            100*mean(sig_slp, na.rm=TRUE),
                            100*mean(sig_sst, na.rm=TRUE)),
    N_bootstrap = N_BOOT, block_size_months = BLOCK_SIZE, drought_months = N_DROUGHT
  ), file.path(out_b3,"bootstrap_significance_summary.csv"), row.names=FALSE)
  
  # Stippled significance maps
  plot_with_stippling <- function(anom_rast, sig_rast, variable, units, col_lim,
                                  col_palette = "RdBu", reverse = TRUE) {
    anom_df <- as.data.frame(anom_rast, xy = TRUE); colnames(anom_df)[3] <- "value"
    sig_df  <- as.data.frame(sig_rast,  xy = TRUE); colnames(sig_df)[3]  <- "sig"
    stip_df <- sig_df %>% filter(sig == 1)
    ggplot() +
      geom_raster(data = anom_df, aes(x = x, y = y, fill = value)) +
      geom_point(data = stip_df, aes(x = x, y = y), size = 0.3, colour = "black", alpha = 0.5) +
      scale_fill_distiller(palette = col_palette, direction = if (reverse) -1 else 1,
                           limits = c(-col_lim, col_lim), oob = scales::squish,
                           name = sprintf("%s\n(%s)", variable, units)) +
      coord_fixed() +
      labs(title    = sprintf("2022\u20132025 %s Anomaly", variable),
           subtitle = sprintf("Dots = FDR-corrected significant (q<0.05, N=%d bootstrap)", N_BOOT),
           x = "Longitude", y = "Latitude") +
      theme_minimal(base_size = 11)
  }
  
  sig_z500_rast <- make_pval_rast(as.numeric(sig_z500), obs_z500)
  sig_slp_rast  <- make_pval_rast(as.numeric(sig_slp),  obs_slp)
  sig_sst_rast  <- make_pval_rast(as.numeric(sig_sst),  obs_sst)
  
  pdf(file.path(out_b3,"fig_bootstrap_significance_maps.pdf"), width=14, height=12)
  gridExtra::grid.arrange(
    plot_with_stippling(obs_z500, sig_z500_rast, "Z500", "m",   col_lim=50),
    plot_with_stippling(obs_slp,  sig_slp_rast,  "SLP",  "hPa", col_lim=2),
    plot_with_stippling(obs_sst,  sig_sst_rast,  "SST",  "\u00b0C", col_lim=1.5),
    ncol = 1,
    top = "2022\u20132025 Drought Fingerprint: Bootstrap Significance (FDR q<0.05)")
  dev.off()
  
  # Historical analogues
  cat("\nComputing historical analogues...\n")
  analogue_r <- vapply(seq_along(bg_idx), function(i) {
    bg_z500 <- z500_bg_mat[, i]
    ok <- !is.na(bg_z500) & !is.na(obs_z500_vec)
    if (sum(ok) < 50) return(NA_real_)
    cor(bg_z500[ok], obs_z500_vec[ok])
  }, numeric(1))
  
  data.frame(date = dates_era5[bg_idx], year = year_nums[bg_idx],
             spatial_r = analogue_r) %>%
    arrange(desc(spatial_r)) %>%
    mutate(rank = row_number()) %>%
    write.csv(file.path(out_b3,"historical_analogues_z500.csv"), row.names=FALSE)
  
  cat("✓ B3 bootstrap analysis complete.\n")
}

cat("\n========== w10a COMPLETE ==========\n")
cat("  composite_results/  — B1 composites\n")
cat("  bootstrap_results/  — B3 significance maps\n")
cat("  Next: run w10b_index_skill.R (independent, no ERA5 needed)\n\n")