# ============================================================================
# NECHAKO BASIN — TOTAL STORAGE DEFICIT INDEX  (Extended Dual-Pipeline)
#
# PIPELINE 1  (BASIN)
#   Original approach: basin-wide total storage aggregated across all GeoLakes
#   lakes in the Nechako watershed.
#   Input:  BEST_ESTIMATE_basin_total_Mm3.csv
#           Required columns:  date, total_storage_Mm3
#
# PIPELINE 2  (NECHAKO_RESERVOIR)
#   Single-lake NTSDI computed for the Nechako Reservoir only.
#   (Nechako Reservoir / Ootsa Lake Reservoir — created by Kenney Dam, 1952;
#    surface elevation ~850–854 m ASL; area ~890 km²; BC, Canada)
#   Two data-source options (set NECHAKO_DATA_SOURCE below):
#
#   "GEOLAKES"   — filter the GeoLakes / HydroLAKES dataset to the
#                  Nechako Reservoir entry using its HydroLAKES_ID.
#                  Expects same format as Pipeline 1 but for a single lake:
#                  columns  date, storage_Mm3  (or similar, see GEOLAKES_*
#                  settings).
#
#   "LOCAL_CSV"  — import a CSV exported from a local water authority
#                  (e.g. BC Hydro reservoir operations data, BC River
#                  Forecast Centre BCWIS gauge exports).
#                  Two sub-options (LOCAL_DATA_TYPE):
#                    "STORAGE" — CSV already contains storage in Mm³
#                    "LEVEL"   — CSV contains water-surface elevation (m ASL);
#                                storage derived via a piecewise-linear
#                                hypsometric (stage–volume) look-up table
#                                (edit HYPSO_TABLE below with real BC Hydro
#                                stage–volume data).
#
# Based on:
#   Awange et al. (2016) Adv. Water Resources 94, 45–59
#   Yirdaw et al. (2008) J. Hydrology 356, 84–92
#   Farahmand & AghaKouchak (2015) Adv. Water Resources 76, 140–145
#
# OUTPUT (per pipeline that is run):
#   Nechako_TSDI_Pipeline1_Output.csv
#   Nechako_NTSDI_Pipeline2_Output.csv
#   Diagnostic plots for each active pipeline
#   Comparison plot (if both pipelines run)
# ============================================================================

rm(list = ls())

# ============================================================================
# PACKAGES
# ============================================================================

packages <- c("tidyverse", "lubridate", "zoo")

for (p in packages) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p)
    library(p, character.only = TRUE)
  }
}

# ============================================================================
# ── SECTION 1 ──  PIPELINE SELECTOR
# ============================================================================
#
# Set PIPELINE_MODE to one of:
#   "BASIN"             — run Pipeline 1 only
#   "NECHAKO_RESERVOIR" — run Pipeline 2 only
#   "BOTH"              — run both and produce a comparison
#
PIPELINE_MODE <- "BOTH"

# ============================================================================
# ── SECTION 2 ──  PIPELINE 1  SETTINGS  (basin-wide GeoLakes)
# ============================================================================

P1_INPUT_FILE <- "Lakes/Glolakes/nechako_extracted/BEST_ESTIMATE_basin_total_Mm3.csv"
# Required columns: date, total_storage_Mm3

# ============================================================================
# ── SECTION 3 ──  PIPELINE 2  SETTINGS  (Nechako Reservoir only)
# ============================================================================

# --- 3a. Data source ---------------------------------------------------------
#
# "GEOLAKES"   — use a GeoLakes/HydroLAKES file filtered to Nechako Reservoir
# "LOCAL_CSV"  — use a CSV from BC Hydro / BC River Forecast Centre
#
NECHAKO_DATA_SOURCE <- "GEOLAKES"

# --- 3b. GeoLakes option -----------------------------------------------------
#
# File that contains per-lake GeoLakes storage time series.
# Must have columns for lake identifier, date, and storage.
#
GEOLAKES_FILE <- "Lakes/Glolakes/nechako_extracted/LandsatPlusICESat2_long.csv"

# HydroLAKES ID for the Nechako Reservoir.
# To find it: open your GeoLakes shapefile / CSV and look for the Nechako
# Reservoir polygon (centroid near 53.75°N, 126.0°W).
# Replace the placeholder value below with the real ID.
# Column names in the GeoLakes long-format file
GEOLAKES_DATE_COL    <- "date"            # date column
GEOLAKES_STORAGE_COL <- "storage_Mm3"    # storage column (Mm³)
GEOLAKES_LAT_COL     <- "lat"            # latitude column
GEOLAKES_LON_COL     <- "lon"            # longitude column

# Nechako Reservoir identified as file_idx=2437 (highest storage ~208 Mm³ mean)
# lake_id and lake_name are NA in this file; filter by centroid coordinates instead.
NECHAKO_LAT  <- 54.2493   # confirmed from LandsatPlusICESat2_long.csv
NECHAKO_LON  <- -125.7896
NECHAKO_TOL  <- 0.01      # tight tolerance — exact centroid known

# ── Find the Nechako Reservoir lake_id ────────────────────────────────────────
# Run this snippet ONCE in the console to identify the correct ID, then set
# NECHAKO_HYDROID (as a string) to match.
#
#   df_tmp <- read_csv(GEOLAKES_FILE, show_col_types = FALSE)
#   df_tmp %>%
#     filter(!is.na(lake_name)) %>%
#     filter(grepl("nechako|ootsa|kenney", lake_name, ignore.case = TRUE)) %>%
#     dplyr::select(lake_id, lake_name, lat, lon) %>%
#     distinct()
#
# If the above returns nothing, try searching by coordinates:
#   df_tmp %>%
#     filter(between(lat, 53.0, 54.5), between(lon, -127.5, -124.5)) %>%
#     dplyr::select(lake_id, lake_name, lat, lon) %>%
#     distinct()
# ─────────────────────────────────────────────────────────────────────────────

# --- 3c. Local CSV option ----------------------------------------------------
#
# CSV exported from BC Hydro reservoir operations reports, BCWIS, or similar.
#
LOCAL_CSV_FILE    <- "nechako_reservoir_levels.csv"  # adjust path as needed
LOCAL_DATE_COL    <- "date"                          # date column name
LOCAL_DATE_FORMAT <- "%Y-%m-%d"                      # adjust if needed
# e.g. "%d/%m/%Y" for DD/MM/YYYY

# What does the local CSV provide?
#   "STORAGE" — a column already in Mm³  (skip hypsometric conversion)
#   "LEVEL"   — water-surface elevation in metres ASL (convert via HYPSO_TABLE)
LOCAL_DATA_TYPE   <- "LEVEL"

# Column name for storage OR level (depending on LOCAL_DATA_TYPE)
LOCAL_VALUE_COL   <- "level_m_asl"

# --- 3d. Hypsometric (stage–volume) look-up table ----------------------------
#
# Used only when LOCAL_DATA_TYPE == "LEVEL".
#
# Replace the placeholder rows below with actual BC Hydro stage–volume data.
# Source: BC Hydro Reservoir Operations, Nechako Reservoir Rule Curves
#   (contact BC Hydro Water Licence Operations or request from DRSP reports).
#
# The Nechako Reservoir operating range is approximately 850–854 m ASL
# (2790–2800 ft); max. depth 305 m; surface area ~890 km².
#
# Format: data.frame with columns  level_m  and  storage_Mm3
#
HYPSO_TABLE <- data.frame(
  level_m     = c(842.0, 844.0, 846.0, 848.0, 849.0,
                  850.0, 851.0, 852.0, 853.0, 854.0, 855.0),
  storage_Mm3 = c(11000, 13000, 15500, 18500, 20000,
                  21500, 22500, 23000, 23400, 23700, 24000)
  # ↑ PLACEHOLDER VALUES — replace with real BC Hydro stage–volume curve data
)

# ============================================================================
# ── SECTION 4 ──  SHARED ALGORITHM SETTINGS
# ============================================================================

# Palmer drought severity class used to anchor p and q
# –1 mild | –2 moderate | –3 severe | –4 extreme
DROUGHT_CLASS_C <- -3

# Optional rolling smoothing of storage before computing TSD
USE_SMOOTHING  <- FALSE
SMOOTH_WINDOW  <- 3   # months (odd integer recommended)

# Blom plotting-position constant (standard for nonparametric SPI-style PIT)
BLOM_A <- 0.44

# ============================================================================
# ── SECTION 5 ──  HELPER FUNCTIONS
# ============================================================================

# ---- 5a. Level → storage conversion via piecewise linear interpolation -----

level_to_storage <- function(level_vec, hypso = HYPSO_TABLE) {
  # Validate table
  stopifnot(
    is.data.frame(hypso),
    "level_m"     %in% names(hypso),
    "storage_Mm3" %in% names(hypso)
  )
  
  hypso <- hypso[order(hypso$level_m), ]   # ensure ascending order
  
  out_of_range <- level_vec < min(hypso$level_m) |
    level_vec > max(hypso$level_m)
  if (any(out_of_range, na.rm = TRUE)) {
    n_oob <- sum(out_of_range, na.rm = TRUE)
    warning(
      n_oob, " level value(s) fall outside the hypsometric table range [",
      min(hypso$level_m), ", ", max(hypso$level_m), " m]. ",
      "Storage will be extrapolated (use with caution) or will return NA."
    )
  }
  
  # approx() performs piecewise-linear interpolation; rule = 1 → NA outside range
  approx(
    x      = hypso$level_m,
    y      = hypso$storage_Mm3,
    xout   = level_vec,
    method = "linear",
    rule   = 1       # NA outside table range; change to 2 to extrapolate
  )$y
}

# ---- 5b. Core TSDI / NTSDI computation  (Awange et al. 2016) ---------------
#
# Arguments:
#   df_in           — data frame with columns: date, storage_Mm3
#   drought_class_c — Palmer class (negative number)
#   use_smoothing   — logical
#   smooth_window   — integer (months)
#   blom_a          — Blom constant
#   label           — character string used in console messages and plots

compute_tsdi <- function(df_in,
                         drought_class_c = DROUGHT_CLASS_C,
                         use_smoothing   = USE_SMOOTHING,
                         smooth_window   = SMOOTH_WINDOW,
                         blom_a          = BLOM_A,
                         label           = "TSDI") {
  
  stopifnot(
    "date"        %in% names(df_in),
    "storage_Mm3" %in% names(df_in)
  )
  
  df <- df_in %>%
    mutate(
      date  = as.Date(date),
      year  = year(date),
      month = month(date)
    ) %>%
    arrange(date)
  
  # ── Optional smoothing ────────────────────────────────────────────────────
  if (use_smoothing) {
    df <- df %>%
      mutate(
        storage_Mm3 = zoo::rollmean(
          storage_Mm3,
          k     = smooth_window,
          fill  = NA,
          align = "center"
        )
      )
    cat(label, ": applied", smooth_window, "-month rolling mean smoothing.\n")
  }
  
  # ── Monthly climatology (mean / max / min per calendar month) ─────────────
  # Follows Awange et al. (2016) Eq. 1
  monthly_stats <- df %>%
    group_by(month) %>%
    summarise(
      clim_mean = mean(storage_Mm3, na.rm = TRUE),
      clim_max  = max( storage_Mm3, na.rm = TRUE),
      clim_min  = min( storage_Mm3, na.rm = TRUE),
      .groups   = "drop"
    )
  
  df <- left_join(df, monthly_stats, by = "month")
  
  # ── Total Storage Deficit  (TSD %) — Eq. 1 ───────────────────────────────
  df <- df %>%
    mutate(
      TSD = (storage_Mm3 - clim_mean) /
        (clim_max    - clim_min  ) * 100
    )
  
  # ── Cumulative TSD ────────────────────────────────────────────────────────
  df <- df %>%
    mutate(cumulative_TSD = cumsum(replace_na(TSD, 0)))
  
  # ── Identify dominant drought period for calibration ─────────────────────
  # Strategy: find the continuously-negative TSD run whose *endpoint* has the
  # most negative cumulative TSD value (deepest deficit).  This is more robust
  # than picking the longest run, which can be very short when the record is
  # dominated by a single wet or dry phase and produces p outside (0,1).
  is_dry <- df$TSD < 0
  r      <- rle(is_dry)
  
  ends   <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  
  dry_runs <- which(r$values == TRUE)
  
  if (length(dry_runs) == 0) {
    stop(label, ": no drought periods detected (no consecutive negative TSD months).")
  }
  
  # Score each dry run by the cumulative TSD at its endpoint (most negative = worst)
  run_scores <- sapply(dry_runs, function(ri) df$cumulative_TSD[ends[ri]])
  deepest_run <- dry_runs[which.min(run_scores)]
  dry_start   <- starts[deepest_run]
  dry_end     <- ends[deepest_run]
  
  drought_df <- df[dry_start:dry_end, ] %>%
    mutate(drought_time = seq_len(n()))
  
  if (nrow(drought_df) < 6) {
    # Fallback: use the longest dry run if the deepest is implausibly short
    longest_run <- dry_runs[which.max(r$lengths[dry_runs])]
    dry_start   <- starts[longest_run]
    dry_end     <- ends[longest_run]
    drought_df  <- df[dry_start:dry_end, ] %>%
      mutate(drought_time = seq_len(n()))
    message(label, ": deepest-deficit episode < 6 months; falling back to longest dry run.")
  }
  
  cat(
    "\n", label, " — dominant drought episode:\n",
    " ", format(min(drought_df$date)), "to", format(max(drought_df$date)),
    "(", nrow(drought_df), "months )\n"
  )
  
  # ── Drought monograph linear regression ──────────────────────────────────
  # Compute episode-local cumulative TSD here, AFTER any fallback selection,
  # so it always applies to the final drought_df (deepest or longest run).
  drought_df <- drought_df %>%
    mutate(local_cum_TSD = cumsum(TSD))
  
  fit <- lm(local_cum_TSD ~ drought_time, data = drought_df)
  m   <- coef(fit)[2]
  b   <- coef(fit)[1]
  
  cat(
    label, " — drought monograph regression:\n",
    "   slope (m) =", round(m, 4), "\n",
    "   intercept (b) =", round(b, 4), "\n"
  )
  
  denom <- m + b
  if (abs(denom) < 1e-6) {
    stop(
      label, ": m + b is near zero (", round(denom, 6), "). ",
      "Cannot compute p and q. ",
      "Ensure the dominant dry period has a clearly negative cumulative trend."
    )
  }
  
  # ── Coefficients p and q  (Eq. 3) ────────────────────────────────────────
  p <- 1 - (m / denom)
  q <- -drought_class_c / denom
  
  cat(
    label, " — derived coefficients:\n",
    "   p =", round(p, 6), "\n",
    "   q =", round(q, 6), "\n"
  )
  
  if (p <= 0 || p >= 1) {
    warning(
      label, ": p = ", round(p, 4), " is outside (0, 1). ",
      "Check the regression period and drought-class setting."
    )
  }
  
  # ── Recursive TSDI  (Eq. 2) ──────────────────────────────────────────────
  n    <- nrow(df)
  TSDI <- rep(NA_real_, n)
  
  # Seed: use first non-NA TSD value (Yirdaw et al. 2008 initial condition)
  first_valid <- which(!is.na(df$TSD))[1]
  if (!is.na(first_valid)) {
    TSDI[first_valid] <- 0.02 * df$TSD[first_valid]
  }
  
  for (i in seq_len(n)[-seq_len(max(1L, first_valid - 1L))][-1L]) {
    if (!is.na(df$TSD[i])) {
      prev <- if (!is.na(TSDI[i - 1])) TSDI[i - 1] else 0
      TSDI[i] <- p * prev + q * df$TSD[i]
    } else {
      # Gap: carry forward the last valid index value unchanged
      TSDI[i] <- TSDI[i - 1]
    }
  }
  
  df[[label]] <- TSDI
  
  # ── PIT normalization → approximately N(0,1) ─────────────────────────────
  # Nonparametric Blom plotting-position (Farahmand & AghaKouchak 2015)
  tsdi_valid <- df %>%
    filter(!is.na(.data[[label]])) %>%
    mutate(
      rank_tsdi      = rank(.data[[label]], ties.method = "average"),
      n_valid        = n(),
      u_blom         = (rank_tsdi - blom_a) / (n_valid + 1 - 2 * blom_a),
      norm_col_value = qnorm(u_blom)
    ) %>%
    rename(!!paste0(label, "_norm") := norm_col_value) %>%
    dplyr::select(date, rank_tsdi, n_valid, u_blom, !!paste0(label, "_norm"))
  
  df <- left_join(df, tsdi_valid, by = "date")
  
  # Return results and calibrated coefficients as a named list
  list(
    data      = df,
    p         = p,
    q         = q,
    m         = m,
    b         = b,
    drought_start = min(drought_df$date),
    drought_end   = max(drought_df$date)
  )
}

# ---- 5c. Diagnostic plots for a single pipeline ----------------------------

plot_pipeline <- function(res, label, basin_label) {
  
  df         <- res$data
  index_col  <- label
  norm_col   <- paste0(label, "_norm")
  p_val      <- res$p
  q_val      <- res$q
  
  # Drought band for annotation
  dband <- data.frame(
    xmin = res$drought_start,
    xmax = res$drought_end,
    ymin = -Inf, ymax = Inf
  )
  
  # ── Plot 1 — Storage ───────────────────────────────────────────────────
  p1 <- ggplot(df, aes(date, storage_Mm3)) +
    geom_rect(
      data = dband,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "lightblue", alpha = 0.3, inherit.aes = FALSE
    ) +
    geom_line(linewidth = 0.8, colour = "#2c7bb6") +
    labs(
      title    = paste(basin_label, "— Lake Storage"),
      subtitle = "Shaded band = dominant drought episode used for calibration",
      x        = NULL,
      y        = "Storage (Mm³)"
    ) +
    theme_bw(base_size = 12)
  
  print(p1)
  
  # ── Plot 2 — TSD (%) ───────────────────────────────────────────────────
  p2 <- ggplot(df, aes(date, TSD)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_ribbon(aes(ymin = pmin(TSD, 0), ymax = 0),
                fill = "#d7191c", alpha = 0.3) +
    geom_ribbon(aes(ymin = 0, ymax = pmax(TSD, 0)),
                fill = "#1a9641", alpha = 0.3) +
    geom_line(linewidth = 0.6) +
    labs(
      title    = paste(basin_label, "— Total Storage Deficit (TSD %)"),
      subtitle = "Red shading = deficit  |  Green shading = surplus",
      x        = NULL,
      y        = "TSD (%)"
    ) +
    theme_bw(base_size = 12)
  
  print(p2)
  
  # ── Plot 3 — Raw index ─────────────────────────────────────────────────
  # Palmer-style thresholds (approximate equivalents)
  thresh <- data.frame(
    yint  = c(-1, -2, -3, -4),
    label = c("Mild (-1)", "Moderate (-2)", "Severe (-3)", "Extreme (-4)")
  )
  
  p3 <- ggplot(df, aes(date, .data[[index_col]])) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_hline(data = thresh,
               aes(yintercept = yint, colour = label),
               linetype = "dotted", linewidth = 0.7) +
    geom_line(linewidth = 0.8) +
    scale_colour_manual(
      values = c("Mild (-1)"     = "#fee08b",
                 "Moderate (-2)" = "#fc8d59",
                 "Severe (-3)"   = "#d73027",
                 "Extreme (-4)"  = "#7b0000"),
      name = "Palmer class"
    ) +
    labs(
      title    = paste(basin_label, "—", label, "(Palmer-style thresholds)"),
      subtitle = paste0("p = ", round(p_val, 4),
                        "  |  q = ", round(q_val, 4),
                        "  |  C = ", DROUGHT_CLASS_C),
      x        = NULL,
      y        = label
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
  
  print(p3)
  
  # ── Plot 4 — Standardised index ────────────────────────────────────────
  spi_thresh <- data.frame(
    yint  = c(-0.5, -1.0, -1.5, -2.0),
    label = c("Near-normal (-0.5)",
              "Moderate (-1.0)",
              "Severe (-1.5)",
              "Extreme (-2.0)")
  )
  
  p4 <- ggplot(df, aes(date, .data[[norm_col]])) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_hline(data = spi_thresh,
               aes(yintercept = yint, colour = label),
               linetype = "dotted", linewidth = 0.7) +
    geom_line(linewidth = 0.8) +
    scale_colour_manual(
      values = c("Near-normal (-0.5)" = "#fee08b",
                 "Moderate (-1.0)"    = "#fc8d59",
                 "Severe (-1.5)"      = "#d73027",
                 "Extreme (-2.0)"     = "#7b0000"),
      name = "SPI-style class"
    ) +
    labs(
      title    = paste(basin_label, "—", norm_col,
                       "(Standardised / SPI-style)"),
      subtitle = "Nonparametric PIT via Blom plotting position → N(0,1)",
      x        = NULL,
      y        = paste0(norm_col, "  [standard normal units]")
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
  
  print(p4)
}

# ---- 5d. Summary console output -------------------------------------------

print_summary <- function(res, label, out_file) {
  
  df       <- res$data
  norm_col <- paste0(label, "_norm")
  
  cat("\n============================================================\n")
  cat(label, "ANALYSIS COMPLETE\n")
  cat("============================================================\n")
  cat("Time range :", format(min(df$date)), "to", format(max(df$date)), "\n")
  cat("Records    :", nrow(df), "\n")
  cat("\nCalibration:\n")
  cat("  p              =", round(res$p, 6), "\n")
  cat("  q              =", round(res$q, 6), "\n")
  cat("  Drought period :", format(res$drought_start),
      "to", format(res$drought_end), "\n")
  cat("\nRaw", label, "range    :",
      round(min(df[[label]],   na.rm = TRUE), 3), "to",
      round(max(df[[label]],   na.rm = TRUE), 3), "\n")
  cat("Standardised range :",
      round(min(df[[norm_col]], na.rm = TRUE), 3), "to",
      round(max(df[[norm_col]], na.rm = TRUE), 3), "\n")
  cat("Output CSV :", out_file, "\n")
}

# ============================================================================
# ── SECTION 6 ──  PIPELINE 1  —  Basin-wide GeoLakes TSDI
# ============================================================================

run_pipeline_1 <- function() {
  
  cat("\n\n========================================================\n")
  cat("PIPELINE 1 — Basin-wide TSDI  (GeoLakes aggregated)\n")
  cat("========================================================\n")
  
  df_raw <- read_csv(P1_INPUT_FILE, show_col_types = FALSE)
  
  cat("Columns found:", paste(names(df_raw), collapse = ", "), "\n")
  
  stopifnot(
    "date"              %in% names(df_raw),
    "total_storage_Mm3" %in% names(df_raw)
  )
  
  df_input <- df_raw %>%
    transmute(
      date        = as.Date(date),
      storage_Mm3 = total_storage_Mm3
    )
  
  res <- compute_tsdi(
    df_in           = df_input,
    drought_class_c = DROUGHT_CLASS_C,
    use_smoothing   = USE_SMOOTHING,
    smooth_window   = SMOOTH_WINDOW,
    blom_a          = BLOM_A,
    label           = "TSDI"
  )
  
  out_file <- "Nechako_TSDI_Pipeline1_Output.csv"
  write_csv(res$data, out_file)
  
  plot_pipeline(res, label = "TSDI", basin_label = "Nechako Basin (all lakes)")
  print_summary(res, label = "TSDI", out_file = out_file)
  
  invisible(res)
}

# ============================================================================
# ── SECTION 7 ──  PIPELINE 2  —  Nechako Reservoir NTSDI
# ============================================================================

run_pipeline_2 <- function() {
  
  cat("\n\n========================================================\n")
  cat("PIPELINE 2 — Nechako Reservoir NTSDI\n")
  cat("  Data source:", NECHAKO_DATA_SOURCE, "\n")
  cat("========================================================\n")
  
  # ── 7a. Load & prepare storage time series ────────────────────────────────
  
  if (NECHAKO_DATA_SOURCE == "GEOLAKES") {
    
    cat("Reading GeoLakes file:", GEOLAKES_FILE, "\n")
    df_geo <- read_csv(GEOLAKES_FILE, show_col_types = FALSE)
    cat("Columns found:", paste(names(df_geo), collapse = ", "), "\n")
    
    # Validate required columns
    missing_cols <- setdiff(
      c(GEOLAKES_DATE_COL, GEOLAKES_STORAGE_COL, GEOLAKES_LAT_COL, GEOLAKES_LON_COL),
      names(df_geo)
    )
    if (length(missing_cols) > 0) {
      stop("Missing columns: ", paste(missing_cols, collapse=", "),
           "\nAvailable: ", paste(names(df_geo), collapse=", "))
    }
    
    # Filter by centroid coordinates (lake_id / lake_name are NA in this file)
    cat(sprintf("  Filtering: lat=%.4f ± %.3f, lon=%.4f ± %.3f\n",
                NECHAKO_LAT, NECHAKO_TOL, NECHAKO_LON, NECHAKO_TOL))
    
    df_nechako <- df_geo %>%
      filter(abs(.data[[GEOLAKES_LAT_COL]] - NECHAKO_LAT) <= NECHAKO_TOL,
             abs(.data[[GEOLAKES_LON_COL]] - NECHAKO_LON) <= NECHAKO_TOL)
    
    if (nrow(df_nechako) == 0) {
      all_coords <- df_geo %>%
        dplyr::select(all_of(c(GEOLAKES_LAT_COL, GEOLAKES_LON_COL))) %>%
        distinct() %>% mutate(across(everything(), ~ round(.x, 4)))
      stop("No rows matched. All coordinates in file:\n",
           paste(capture.output(print(all_coords)), collapse="\n"),
           "\nUpdate NECHAKO_LAT / NECHAKO_LON in Section 3b.")
    }
    cat("  Matched", nrow(df_nechako), "rows\n")
    
    cat("Nechako Reservoir rows found:", nrow(df_nechako), "\n")
    
    df_input <- df_nechako %>%
      transmute(
        date        = as.Date(.data[[GEOLAKES_DATE_COL]]),
        storage_Mm3 = .data[[GEOLAKES_STORAGE_COL]]
      )
    
  } else if (NECHAKO_DATA_SOURCE == "LOCAL_CSV") {
    
    cat("Reading local authority CSV:", LOCAL_CSV_FILE, "\n")
    df_local <- read_csv(LOCAL_CSV_FILE, show_col_types = FALSE)
    cat("Columns found:", paste(names(df_local), collapse = ", "), "\n")
    
    stopifnot(
      LOCAL_DATE_COL  %in% names(df_local),
      LOCAL_VALUE_COL %in% names(df_local)
    )
    
    if (LOCAL_DATA_TYPE == "STORAGE") {
      
      cat("Using direct storage values (Mm³) from column:", LOCAL_VALUE_COL, "\n")
      
      df_input <- df_local %>%
        transmute(
          date        = as.Date(.data[[LOCAL_DATE_COL]], format = LOCAL_DATE_FORMAT),
          storage_Mm3 = .data[[LOCAL_VALUE_COL]]
        )
      
    } else if (LOCAL_DATA_TYPE == "LEVEL") {
      
      cat("Converting water-level (m ASL) → storage via hypsometric table.\n")
      cat("  Level column   :", LOCAL_VALUE_COL, "\n")
      cat("  Hypso. table range:",
          min(HYPSO_TABLE$level_m), "–", max(HYPSO_TABLE$level_m), "m ASL\n")
      
      df_input <- df_local %>%
        transmute(
          date        = as.Date(.data[[LOCAL_DATE_COL]], format = LOCAL_DATE_FORMAT),
          level_m_asl = .data[[LOCAL_VALUE_COL]]
        ) %>%
        mutate(
          storage_Mm3 = level_to_storage(level_m_asl)
        )
      
      n_na_storage <- sum(is.na(df_input$storage_Mm3))
      if (n_na_storage > 0) {
        cat(
          "  WARNING:", n_na_storage,
          "records produced NA storage (level outside hypsometric table range).",
          "Update HYPSO_TABLE or adjust rule= in level_to_storage().\n"
        )
      }
      
    } else {
      stop("LOCAL_DATA_TYPE must be 'STORAGE' or 'LEVEL'. Got: '", LOCAL_DATA_TYPE, "'")
    }
    
  } else {
    stop("NECHAKO_DATA_SOURCE must be 'GEOLAKES' or 'LOCAL_CSV'. Got: '",
         NECHAKO_DATA_SOURCE, "'")
  }
  
  # Filter out NA storage rows before computing index
  n_before <- nrow(df_input)
  df_input <- df_input %>% filter(!is.na(storage_Mm3))
  n_after  <- nrow(df_input)
  if (n_before != n_after) {
    cat("  Removed", n_before - n_after, "rows with NA storage.\n")
  }
  
  # ── 7b. Run TSDI algorithm on single-lake storage ─────────────────────────
  res <- compute_tsdi(
    df_in           = df_input,
    drought_class_c = DROUGHT_CLASS_C,
    use_smoothing   = USE_SMOOTHING,
    smooth_window   = SMOOTH_WINDOW,
    blom_a          = BLOM_A,
    label           = "NTSDI"      # Nechako Total Storage Deficit Index
  )
  
  # Append data source metadata to output
  res$data <- res$data %>%
    mutate(data_source = NECHAKO_DATA_SOURCE)
  
  out_file <- "Nechako_NTSDI_Pipeline2_Output.csv"
  write_csv(res$data, out_file)
  
  plot_pipeline(
    res,
    label       = "NTSDI",
    basin_label = paste0(
      "Nechako Reservoir  [source: ",
      ifelse(NECHAKO_DATA_SOURCE == "LOCAL_CSV",
             paste0("Local CSV (", LOCAL_DATA_TYPE, ")"),
             "GeoLakes"),
      "]"
    )
  )
  print_summary(res, label = "NTSDI", out_file = out_file)
  
  invisible(res)
}

# ============================================================================
# ── SECTION 8 ──  COMPARISON PLOT  (only when PIPELINE_MODE == "BOTH")
# ============================================================================

plot_comparison <- function(res1, res2) {
  
  cat("\n\n== Generating comparison plot ==\n")
  
  comp <- inner_join(
    res1$data %>% dplyr::select(date, TSDI, TSDI_norm),
    res2$data %>% dplyr::select(date, NTSDI, NTSDI_norm),
    by = "date"
  )
  
  if (nrow(comp) == 0) {
    warning("No overlapping dates between Pipeline 1 and Pipeline 2. Skipping comparison plot.")
    return(invisible(NULL))
  }
  
  cat("Overlapping records:", nrow(comp), "\n")
  
  # ── Raw index comparison ──────────────────────────────────────────────────
  comp_long_raw <- comp %>%
    pivot_longer(
      cols      = c(TSDI, NTSDI),
      names_to  = "Index",
      values_to = "value"
    ) %>%
    mutate(
      Index = recode(Index,
                     TSDI  = "Pipeline 1: Basin-wide TSDI",
                     NTSDI = "Pipeline 2: Nechako Reservoir NTSDI"
      )
    )
  
  pc1 <- ggplot(comp_long_raw, aes(date, value, colour = Index)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_line(linewidth = 0.8, alpha = 0.85) +
    scale_colour_manual(
      values = c("Pipeline 1: Basin-wide TSDI"           = "#2c7bb6",
                 "Pipeline 2: Nechako Reservoir NTSDI"   = "#d7191c")
    ) +
    labs(
      title    = "Comparison — Basin TSDI vs. Nechako Reservoir NTSDI",
      subtitle = paste0(
        "Overlap: ", format(min(comp$date)), " – ", format(max(comp$date)),
        "  |  Data source (P2): ", NECHAKO_DATA_SOURCE
      ),
      x      = NULL,
      y      = "Index value",
      colour = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
  
  print(pc1)
  
  # ── Standardised index comparison ─────────────────────────────────────────
  comp_long_norm <- comp %>%
    pivot_longer(
      cols      = c(TSDI_norm, NTSDI_norm),
      names_to  = "Index",
      values_to = "value"
    ) %>%
    mutate(
      Index = recode(Index,
                     TSDI_norm  = "Pipeline 1: Basin-wide TSDI_norm",
                     NTSDI_norm = "Pipeline 2: Nechako Reservoir NTSDI_norm"
      )
    )
  
  pc2 <- ggplot(comp_long_norm, aes(date, value, colour = Index)) +
    geom_hline(yintercept = 0,    linetype = "dashed",  colour = "grey50") +
    geom_hline(yintercept = -1.0, linetype = "dotted",  colour = "#fc8d59") +
    geom_hline(yintercept = -1.5, linetype = "dotted",  colour = "#d73027") +
    geom_line(linewidth = 0.8, alpha = 0.85) +
    scale_colour_manual(
      values = c("Pipeline 1: Basin-wide TSDI_norm"            = "#2c7bb6",
                 "Pipeline 2: Nechako Reservoir NTSDI_norm"    = "#d7191c")
    ) +
    labs(
      title    = "Comparison — Standardised Basin TSDI_norm vs. NTSDI_norm",
      subtitle = "Dotted lines: SPI-style moderate (–1.0) and severe (–1.5) thresholds",
      x        = NULL,
      y        = "Standardised index [N(0,1)]",
      colour   = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
  
  print(pc2)
  
  # ── Scatter / correlation ─────────────────────────────────────────────────
  corr_val <- round(cor(comp$TSDI_norm, comp$NTSDI_norm, use = "complete.obs"), 3)
  
  pc3 <- ggplot(comp, aes(TSDI_norm, NTSDI_norm)) +
    geom_point(alpha = 0.4, size = 1.5, colour = "#404040") +
    geom_smooth(method = "lm", se = TRUE, colour = "#2c7bb6", linewidth = 0.9) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4, size = 4,
             label = paste0("r = ", corr_val)) +
    labs(
      title    = "Scatter — Basin TSDI_norm vs. Nechako Reservoir NTSDI_norm",
      subtitle = "Dashed line = 1:1; blue line = OLS fit",
      x        = "Basin TSDI_norm (Pipeline 1)",
      y        = "Nechako Reservoir NTSDI_norm (Pipeline 2)"
    ) +
    theme_bw(base_size = 12)
  
  print(pc3)
  
  cat("Pearson r (TSDI_norm vs NTSDI_norm):", corr_val, "\n")
}

# ============================================================================
# ── SECTION 9 ──  MAIN EXECUTION
# ============================================================================

res_p1 <- NULL
res_p2 <- NULL

if (PIPELINE_MODE %in% c("BASIN", "BOTH")) {
  res_p1 <- run_pipeline_1()
}

if (PIPELINE_MODE %in% c("NECHAKO_RESERVOIR", "BOTH")) {
  res_p2 <- run_pipeline_2()
}

if (PIPELINE_MODE == "BOTH" && !is.null(res_p1) && !is.null(res_p2)) {
  plot_comparison(res_p1, res_p2)
}

cat("\n\nDone.\n")