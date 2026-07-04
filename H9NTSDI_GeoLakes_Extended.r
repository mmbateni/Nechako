# ============================================================================
# NECHAKO BASIN ??? TOTAL STORAGE DEFICIT INDEX  (Extended Dual-Pipeline)
# VERSION 2 ??? Changes from original (H9NTSDI_GeoLakes_Extended.r):
#
#   1. Pipeline 2 now supports LOCAL_DATA ONLY.
#      The GEOLAKES branch has been removed; NECHAKO_DATA_SOURCE is fixed to
#      "LOCAL_DATA" and its sub-options (STORAGE / LEVEL / LEVEL_DAILY) remain.
#
#   2. compute_tsdi() ??? p-validity retry loop added.
#      When the primary calibration episode yields p ??? (0,1) (which causes the
#      recursive formula to diverge), the function now iterates over all dry
#      runs in descending order of deficit severity until one produces a valid
#      p.  This fixes the Pipeline 1 numerical explosion seen in the original
#      run (p = ???7.07 ??? raw TSDI ??? ??1e182).
#
#   3. New Section 3.5  ??? Kenney Dam / Nechako Reservoir regulation settings.
#
#   4. New function check_regulated_lakes() (Section 5e).
#      Reads the GeoLakes long file, filters to the Nechako watershed bounding
#      box, flags each lake as IMPOUNDED / LIKELY IMPOUNDED / DOWNSTREAM /
#      NATURAL, and produces a summary CSV + diagnostic PDF.
#
#   5. Section 9 now calls check_regulated_lakes() when pipeline mode includes
#      BASIN or BOTH.
#
# Based on:
#   Awange et al. (2016) Adv. Water Resources 94, 45???59
#   Yirdaw et al. (2008) J. Hydrology 356, 84???92
#   Farahmand & AghaKouchak (2015) Adv. Water Resources 76, 140???145
# ============================================================================
rm(list = ls())

# --- OUTPUT DIRECTORY SETUP ---
OUT_DIR <- "ntsdi_results"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
sink(file.path(OUT_DIR, "console_log.txt"), append = TRUE, split = TRUE)

# ============================================================================
# PACKAGES
# ============================================================================
packages <- c("tidyverse", "lubridate", "zoo", "ggplot2", "dplyr", "readxl")
for (p in packages) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p)
    library(p, character.only = TRUE)
  }
}

# ============================================================================
# HELPER ??? standard time-series plot (shared across both pipelines)
# ============================================================================
plot_standard_basin_ts <- function(df, index_name, color = "steelblue",
                                   drought_threshold = -0.5) {
  if (!all(c("date", "value") %in% names(df)))
    stop("df must have 'date' and 'value' columns.")
  df <- df[order(df$date), ]
  df <- df[!is.na(df$value), ]
  
  in_drought <- df$value < drought_threshold
  rle_res    <- rle(in_drought)
  # FIX (bug #4): the old condition '& rle_res$lengths >= 1' was dead code --
  # rle() run lengths are always >= 1 by definition, so it never filtered
  # anything. No minimum-run-length filter on drought-event counting is
  # intended here, so we just count every contiguous drought run directly.
  n_ev       <- sum(rle_res$values)
  
  y_range <- max(df$value, na.rm = TRUE) - min(df$value, na.rm = TRUE)
  y_pad   <- y_range * 0.10
  y_lo    <- min(-3.5, min(df$value, na.rm = TRUE) - y_pad)
  y_hi    <- max( 3.5, max(df$value, na.rm = TRUE) + y_pad)
  
  drought_bands <- list(
    ggplot2::annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf,  ymax=-2.0, fill="#7B0025", alpha=0.15),
    ggplot2::annotate("rect", xmin=-Inf, xmax=Inf, ymin=-2.0,  ymax=-1.5, fill="#D73027", alpha=0.15),
    ggplot2::annotate("rect", xmin=-Inf, xmax=Inf, ymin=-1.5,  ymax=-1.0, fill="#F46D43", alpha=0.15),
    ggplot2::annotate("rect", xmin=-Inf, xmax=Inf, ymin=-1.0,  ymax= 1.0, fill="#2E7D32", alpha=0.15),
    ggplot2::annotate("rect", xmin=-Inf, xmax=Inf, ymin= 1.0,  ymax= 1.5, fill="#90EE90", alpha=0.15),
    ggplot2::annotate("rect", xmin=-Inf, xmax=Inf, ymin= 1.5,  ymax= 2.0, fill="#66BB6A", alpha=0.15),
    ggplot2::annotate("rect", xmin=-Inf, xmax=Inf, ymin= 2.0,  ymax= Inf, fill="#4575B4", alpha=0.15)
  )
  
  ggplot2::ggplot(df, ggplot2::aes(x = date, y = value)) +
    drought_bands +
    ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = c(-1, 1), linetype = "dashed",
                        color = "gray50", linewidth = 0.3) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = ifelse(value < drought_threshold, value, 0), ymax = 0),
      fill = "#c0392b", alpha = 0.35) +
    ggplot2::geom_line(color = color, linewidth = 0.6) +
    ggplot2::scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
    ggplot2::scale_y_continuous(limits = c(y_lo, y_hi)) +
    ggplot2::labs(
      title    = sprintf("%s Basin-Averaged Time Series", toupper(index_name)),
      subtitle = sprintf("%d drought events detected (onset < %.1f)", n_ev, drought_threshold),
      x = "Year", y = sprintf("%s Index Value", toupper(index_name))
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "gray30"),
      axis.title    = ggplot2::element_text(face = "bold"),
      axis.text.x   = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "bottom"
    )
}

# ============================================================================
# HELPER — standalone basin time series PNG (shared across both pipelines)
# ============================================================================
# Mirrors plot_index_timeseries()/plot_ssi_timeseries() in the companion
# H1SSI_SSMI_ERALand.R, H7WSC_StreamFlow.R, and H8SWEI_Snow.R scripts, so all
# four families of standardized index (SSI/SSMI, SWEI/SSPI-1, TSDI/NTSDI)
# share the same onset convention, drought-shading rule, and visual language:
#   - blue/green/red severity background bands
#   - dashed reference lines at +/-1, dotted at +/-1.5 and +/-2
#   - drought shading ONLY on runs that cross the onset threshold (default
#     -0.5, matching SSI_DROUGHT_THRESHOLD elsewhere), filled all the way up
#     to zero rather than stopped at onset
#   - subtitle reporting the number of drought events at that onset
# Unlike plot_standard_basin_ts() (used inline in the multi-panel diagnostic
# PDF via plot_pipeline()), this writes its own standalone PNG to
# out_dir/figures/, matching how the other three scripts' plots are saved.
plot_index_timeseries <- function(dates_vec, values_vec, index_label,
                                  title_label, out_dir, onset = -0.5) {
  df <- data.frame(date = dates_vec, value = values_vec)
  df <- df[is.finite(df$value), ]
  if (nrow(df) == 0L) {
    cat(sprintf("  [PLOT] Skipping %s time series: no finite values\n", index_label))
    return(invisible(NULL))
  }
  
  # Drought event count (contiguous runs below onset threshold)
  below_onset <- df$value < onset
  below_onset[is.na(below_onset)] <- FALSE
  n_events <- sum(rle(below_onset)$values)
  
  # Shade only the runs that actually cross the onset threshold, but filled
  # all the way up to zero (not stopped at onset) -- same rule used in the
  # companion SSI/SSMI, SWEI/SSPI-1 scripts. Separate contiguous runs so the
  # ribbon fill doesn't bridge unrelated dry spells.
  df$below <- below_onset
  df$grp   <- cumsum(c(1L, diff(as.integer(df$below)) != 0L))
  
  y_lo <- min(-2.2, floor(min(df$value) * 10) / 10)
  y_hi <- max( 2.2, ceiling(max(df$value) * 10) / 10)
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = date, y = value)) +
    
    # ---- severity background bands ----
  ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1,    ymax = Inf,
                    fill = "#4472C4", alpha = 0.10) +
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0,    ymax = 1,
                      fill = "#70AD47", alpha = 0.08) +
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = -1,   ymax = 0,
                      fill = "#70AD47", alpha = 0.05) +
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -1,
                      fill = "#C00000", alpha = 0.08) +
    
    # ---- reference lines ----
  ggplot2::geom_hline(yintercept = 0,             color = "black", linewidth = 0.6) +
    ggplot2::geom_hline(yintercept = c(-1, 1),      color = "grey30", linetype = "dashed",
                        linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = c(-1.5, 1.5),  color = "grey55", linetype = "dotted",
                        linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = c(-2, 2),      color = "grey55", linetype = "dotted",
                        linewidth = 0.4) +
    
    # ---- drought shading: ONLY months that cross the onset threshold,
    #      filled all the way up to zero (not stopped at onset) ----
  ggplot2::geom_ribbon(data = subset(df, below),
                       ggplot2::aes(ymin = value, ymax = 0, group = grp),
                       fill = "#C0504D", alpha = 0.45) +
    
    # ---- the series itself ----
  ggplot2::geom_line(color = "#1F77B4", linewidth = 0.45) +
    
    ggplot2::scale_x_date(date_breaks = "5 years", date_labels = "%Y",
                          expand = ggplot2::expansion(mult = c(0.01, 0.015))) +
    ggplot2::scale_y_continuous(breaks = seq(-4, 4, 1)) +
    ggplot2::coord_cartesian(ylim = c(y_lo, y_hi)) +
    
    ggplot2::labs(
      title    = sprintf("%s Time Series  |  NECHAKO BASIN", title_label),
      subtitle = sprintf("%d drought events detected  (onset < %.2f)",
                         n_events, onset),
      x = "Year", y = sprintf("%s Index Value", index_label)
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      panel.grid.minor    = ggplot2::element_blank(),
      panel.grid.major.x  = ggplot2::element_line(color = "grey85"),
      panel.grid.major.y  = ggplot2::element_blank(),
      plot.title          = ggplot2::element_text(face = "bold", size = 15),
      plot.subtitle       = ggplot2::element_text(color = "grey40", size = 11),
      axis.text.x         = ggplot2::element_text(angle = 30, hjust = 1),
      plot.margin         = ggplot2::margin(10, 16, 10, 10)
    )
  
  figures_dir <- file.path(out_dir, "figures")
  if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)
  fname <- file.path(figures_dir, sprintf("%s_basin_timeseries.png", tolower(index_label)))
  ggplot2::ggsave(fname, p, width = 11, height = 5.5, dpi = 200, bg = "white")
  cat(sprintf("  [PLOT] %s time series saved: %s  (%d events)\n",
              index_label, basename(fname), n_events))
  invisible(p)
}

# ============================================================================
# ?????? SECTION 1 ??????  PIPELINE SELECTOR
# ============================================================================
VALID_PIPELINE_MODES <- c("BASIN", "NECHAKO_RESERVOIR", "BOTH")

# ---- To use in batch/scheduled jobs, set DEFAULT_PIPELINE_MODE below ----
DEFAULT_PIPELINE_MODE <- "BOTH"

# Normalizes any accepted form of the answer ("1"/"2"/"3" or the mode name,
# any case/whitespace) to one of VALID_PIPELINE_MODES, or NA if invalid.
normalize_pipeline_answer <- function(ans) {
  ans <- toupper(trimws(ans))
  ans <- switch(ans,
                "1" = "BASIN",
                "2" = "NECHAKO_RESERVOIR",
                "3" = "BOTH",
                ans)
  if (ans %in% VALID_PIPELINE_MODES) ans else NA_character_
}

prompt_pipeline_mode <- function() {
  
  menu_text <- paste0(
    "Which pipeline should be run?\n\n",
    "  1) BASIN\n",
    "  2) NECHAKO_RESERVOIR\n",
    "  3) BOTH\n"
  )
  
  # Priority 1: command-line argument (e.g. Rscript H9NTSDI_GeoLakes_Extended.r BOTH)
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1 && !is.na(normalize_pipeline_answer(args[1]))) {
    mode <- normalize_pipeline_answer(args[1])
    cat(sprintf("\u2192 Pipeline set from command-line argument: %s\n", mode))
    return(mode)
  }
  
  # Priority 2: RStudio GUI dialog (pops a window - works with the Run button,
  # not just source()). This is what fires when you click "Run"/"Source" in
  # the RStudio editor rather than typing source() into the console.
  if (interactive() && requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    raw_input <- rstudioapi::showPrompt(
      title   = "Pipeline Selection",
      message = paste0(menu_text, "\nEnter BASIN / NECHAKO_RESERVOIR / BOTH (or 1/2/3):"),
      default = DEFAULT_PIPELINE_MODE
    )
    if (is.null(raw_input)) {
      # User clicked Cancel -> use the default
      cat(sprintf("\u2192 Dialog cancelled: using default pipeline %s\n", DEFAULT_PIPELINE_MODE))
      return(DEFAULT_PIPELINE_MODE)
    }
    mode <- normalize_pipeline_answer(raw_input)
    while (is.na(mode)) {
      raw_input <- rstudioapi::showPrompt(
        title   = "Invalid input \u2014 try again",
        message = "Please enter BASIN, NECHAKO_RESERVOIR, BOTH, or 1/2/3:",
        default = DEFAULT_PIPELINE_MODE
      )
      if (is.null(raw_input)) {
        cat(sprintf("\u2192 Dialog cancelled: using default pipeline %s\n", DEFAULT_PIPELINE_MODE))
        return(DEFAULT_PIPELINE_MODE)
      }
      mode <- normalize_pipeline_answer(raw_input)
    }
    cat(sprintf("\u2192 Pipeline set from RStudio dialog: %s\n", mode))
    return(mode)
  }
  
  # Priority 3: plain readline() fallback (R console, non-RStudio interactive,
  # or RStudio without the rstudioapi package available)
  if (interactive()) {
    cat("\n", menu_text, sep = "")
    mode <- NA_character_
    while (is.na(mode)) {
      ans  <- readline(prompt = "Enter BASIN / NECHAKO_RESERVOIR / BOTH (or 1/2/3): ")
      mode <- normalize_pipeline_answer(ans)
      if (is.na(mode)) cat("Invalid selection '", ans, "'. Please try again.\n", sep = "")
    }
    return(mode)
  }
  
  # Priority 4: non-interactive fallback (Rscript piped from a scheduler,
  # no CLI argument supplied) - try stdin, else fall back to the default.
  cat("\n", menu_text, sep = "")
  cat("Enter BASIN / NECHAKO_RESERVOIR / BOTH (or 1/2/3): ")
  ans  <- tryCatch(readLines(con = "stdin", n = 1), error = function(e) character(0))
  mode <- if (length(ans) >= 1) normalize_pipeline_answer(ans[1]) else NA_character_
  if (is.na(mode)) {
    cat(sprintf("\u2192 No valid input received: using default pipeline %s\n", DEFAULT_PIPELINE_MODE))
    cat("  (Edit DEFAULT_PIPELINE_MODE at the top of the script to change this.)\n")
    mode <- DEFAULT_PIPELINE_MODE
  }
  mode
}

PIPELINE_MODE <- prompt_pipeline_mode()
cat("\nSelected PIPELINE_MODE:", PIPELINE_MODE, "\n\n")

# ============================================================================
# ?????? SECTION 2 ??????  PIPELINE 1  SETTINGS  (basin-wide GeoLakes)
# ============================================================================
P1_INPUT_FILE <- "Lakes/Glolakes/nechako_extracted/BEST_ESTIMATE_basin_total_Mm3.csv"
# Required columns: date, total_storage_Mm3

# ============================================================================
# ?????? SECTION 3 ??????  PIPELINE 2  SETTINGS  (Nechako Reservoir ??? LOCAL CSV ONLY)
# ============================================================================
# NOTE (v2): Pipeline 2 now uses LOCAL CSV exclusively.
#   The GEOLAKES branch has been removed.  Set LOCAL_DATA_TYPE and the
#   relevant column / file settings below.

NECHAKO_DATA_SOURCE <- "LOCAL_DATA"   # FIXED ??? do not change

# Path to your daily level / storage file (BC Hydro, BCWIS, RFC, etc.)
# NOTE: despite the .csv extension, this file is whitespace-delimited
# (Year, DOY, Level_m — 3 fields, no header), NOT comma-separated.
# read_table() is used deliberately in run_pipeline_2(); do not switch to read_csv().
LOCAL_DATA_FILE    <- "D:/Nechako_Drought/Nechako/Lakes/dailylevel.txt"

# Choose data type:
#   "STORAGE"      ??? CSV already contains storage in Mm??
#   "LEVEL"        ??? CSV has a WSE column (m ASL); storage via hypsometric LUT
#   "LEVEL_DAILY"  ??? CSV has Year, Day-of-Year, Level_m (no header row)
LOCAL_DATA_TYPE   <- "LEVEL_DAILY"

# Column names (ignored when LOCAL_DATA_TYPE == "LEVEL_DAILY")
LOCAL_VALUE_COL   <- "Level_m"
LOCAL_DATE_COL    <- "Date"
LOCAL_DATE_FORMAT <- "%Y-%m-%d"

# GeoLakes long-format file (used only for the lake regulation diagnostic;
# Pipeline 2 no longer reads this file at all).
GEOLAKES_FILE        <- "Lakes/Glolakes/nechako_extracted/LandsatPlusICESat2_long.csv"
GEOLAKES_DATE_COL    <- "date"
GEOLAKES_STORAGE_COL <- "storage_Mm3"
GEOLAKES_LAT_COL     <- "lat"
GEOLAKES_LON_COL     <- "lon"

# Nechako Reservoir GeoLakes centroid (used for regulation diagnostic only)
NECHAKO_LAT  <- 54.2493
NECHAKO_LON  <- -125.7896
NECHAKO_TOL  <- 0.01

# Storage curve (stage???volume look-up table for LEVEL / LEVEL_DAILY)
Storage_curve <- read_excel("./Lakes/Storage curve.xlsx")
col_elev <- grep("Elevation", names(Storage_curve), value = TRUE)[1]
col_stor <- grep("Storage",   names(Storage_curve), value = TRUE)[1]

if (is.na(col_elev) || is.na(col_stor))
  stop("Could not find Elevation/Storage columns. Found names: ",
       paste(names(Storage_curve), collapse = ", "))

level_m    <- Storage_curve[[col_elev]] * 0.3048
storage_m3 <- Storage_curve[[col_stor]] * 1000000

HYPSO_TABLE <- data.frame(level_m, storage_m3) %>% arrange(level_m)

# ============================================================================
# ?????? SECTION 3.5 ??????  KENNEY DAM / REGULATION CHECK SETTINGS  (NEW in v2)
# ============================================================================
# Kenney Dam approximate coordinates (BC, Canada)
KENNEY_DAM_LAT <- 53.733    # ~53??44'N
KENNEY_DAM_LON <- -127.217  # ~127??13'W

# Nechako watershed coarse bounding box for the lake regulation scan
NECHAKO_WS_LAT  <- c(52.5, 55.5)
NECHAKO_WS_LON  <- c(-128.5, -122.0)

# Nechako Reservoir impoundment zone
# The reservoir (Ootsa / Whitesail / Eutsuk / Knewstubb / Natalkuz) lies WEST
# of the dam between approximately 53.1???54.6??N, 125.0???127.8??W at ???840 m ASL.
RESERVOIR_ZONE <- list(
  lat = c(53.1, 54.6),
  lon = c(-127.8, -125.0),
  elev_min_m = 840          # reservoir operating range ??? ~840 m
)

# Known sub-basins of the Nechako Reservoir system (for annotation)
RESERVOIR_COMPONENTS <- data.frame(
  name = c("Ootsa Lake / Nechako Res.", "Whitesail Lake",
           "Eutsuk Lake", "Knewstubb Lake", "Natalkuz Lake"),
  lat  = c(53.60, 53.55, 53.35, 53.68, 53.78),
  lon  = c(-126.00, -126.95, -127.25, -126.30, -126.55)
)

# ============================================================================
# ?????? SECTION 4 ??????  SHARED ALGORITHM SETTINGS
# ============================================================================
DROUGHT_CLASS_C <- -3      # Palmer class: ???1 mild | ???2 moderate | ???3 severe | ???4 extreme
USE_SMOOTHING   <- FALSE
SMOOTH_WINDOW   <- 3       # months
BLOM_A          <- 0.44    # Blom plotting-position constant

# ============================================================================
# ?????? SECTION 5 ??????  HELPER FUNCTIONS
# ============================================================================

# ---- 5a. Level ??? storage (piecewise linear) ---------------------------------
level_to_storage <- function(level_vec, hypso = HYPSO_TABLE) {
  stopifnot(
    is.data.frame(hypso),
    "level_m"    %in% names(hypso),
    "storage_m3" %in% names(hypso)
  )
  hypso <- hypso[order(hypso$level_m), ]
  out_of_range <- level_vec < min(hypso$level_m) |
    level_vec > max(hypso$level_m)
  if (any(out_of_range, na.rm = TRUE)) {
    warning(sum(out_of_range, na.rm = TRUE),
            " level value(s) outside hypsometric table range [",
            min(hypso$level_m), ", ", max(hypso$level_m), " m].")
  }
  approx(x = hypso$level_m, y = hypso$storage_m3,
         xout = level_vec, method = "linear", rule = 1)$y
}

# ---- 5b. Core TSDI / NTSDI computation  (Awange et al. 2016) ---------------
# CHANGE v2: added p-validity retry loop.
#   When the primary drought episode yields p ??? (0,1), the function now
#   iterates over all dry runs sorted by deficit severity, trying each until
#   one gives 0 < p < 1.  This prevents the recursive formula from diverging.
#
# Root cause of the original Pipeline 1 failure
# -----------------------------------------------
# Dominant episode (fallback to longest run): 1991-07 to 1992-06 (6 months)
#   slope  m = ???34.054  intercept b = +29.833
#   denom  m + b = ???4.221  (barely negative; intercept nearly cancels slope)
#   p = 1 ??? m/(m+b) = 1 ??? (???34.054/???4.221) = 1 ??? 8.069 = ???7.069
# With p = ???7.07 the recursion TSDI[i] = p??TSDI[i-1] + q??TSD[i] diverges
# exponentially ??? raw TSDI ??? ??1e182.  The standardised range looks reasonable
# only because the rank-based (Blom) normalisation still works on any finite
# ordering of values, masking the underlying explosion.
compute_tsdi <- function(df_in,
                         drought_class_c = DROUGHT_CLASS_C,
                         use_smoothing   = USE_SMOOTHING,
                         smooth_window   = SMOOTH_WINDOW,
                         blom_a          = BLOM_A,
                         label           = "TSDI") {
  stopifnot("date"       %in% names(df_in),
            "storage_m3" %in% names(df_in))
  
  df <- df_in %>%
    mutate(date  = as.Date(date),
           year  = year(date),
           month = month(date)) %>%
    arrange(date)
  
  # Optional smoothing
  if (use_smoothing) {
    df <- df %>%
      mutate(storage_m3 = zoo::rollmean(storage_m3, k = smooth_window,
                                        fill = NA, align = "center"))
    cat(label, ": applied", smooth_window, "-month rolling mean.\n")
  }
  
  # Monthly climatology (Awange et al. 2016, Eq. 1)
  # NOTE: bimonthly GeoLakes data has ~6 obs/year, so some calendar months
  # have ALL NA storage.  suppressWarnings() silences max/min on empty groups;
  # the subsequent mutate replaces -Inf / Inf / NaN with NA so those months
  # produce NA TSD rather than exploding the calculation.
  monthly_stats <- df %>%
    group_by(month) %>%
    summarise(
      clim_mean = suppressWarnings(mean(storage_m3, na.rm = TRUE)),
      clim_max  = suppressWarnings(max( storage_m3, na.rm = TRUE)),
      clim_min  = suppressWarnings(min( storage_m3, na.rm = TRUE)),
      .groups   = "drop"
    ) %>%
    mutate(
      clim_mean = ifelse(is.nan(clim_mean)       | is.infinite(clim_mean), NA_real_, clim_mean),
      clim_max  = ifelse(is.infinite(clim_max),  NA_real_, clim_max),
      clim_min  = ifelse(is.infinite(clim_min),  NA_real_, clim_min)
    )
  df <- left_join(df, monthly_stats, by = "month")
  
  # TSD (%)
  df <- df %>%
    mutate(TSD = (storage_m3 - clim_mean) / (clim_max - clim_min) * 100)
  
  # Cumulative TSD
  df <- df %>%
    mutate(cumulative_TSD = cumsum(replace_na(TSD, 0)))
  
  # Identify all dry runs
  df_valid <- df[!is.na(df$TSD), ]
  is_dry   <- df_valid$TSD < 0
  r        <- rle(is_dry)
  ends     <- cumsum(r$lengths)
  starts   <- ends - r$lengths + 1
  dry_runs <- which(r$values == TRUE)
  
  if (length(dry_runs) == 0)
    stop(label, ": no drought periods detected (no negative TSD months).")
  
  # ?????? Internal helper: attempt calibration on one dry run ?????????????????????????????????????????????????????????
  try_calibrate <- function(run_idx) {
    ds <- starts[run_idx]; de <- ends[run_idx]
    if ((de - ds + 1) < 6) return(NULL)   # need ≥6 points for regression
    d <- df_valid[ds:de, ] %>%
      mutate(drought_time = as.numeric(date - min(date)) / 30.44,  # months, real elapsed time
             local_cum_TSD = cumsum(TSD))
    fit <- lm(local_cum_TSD ~ drought_time, data = d)
    m_val  <- coef(fit)[2]; b_val <- coef(fit)[1]
    denom  <- m_val + b_val
    if (abs(denom) < 1e-6) return(NULL)
    p_val  <- 1 - (m_val / denom)
    q_val  <- -drought_class_c / denom
    list(drought_df = d, m = m_val, b = b_val, denom = denom,
         p = p_val, q = q_val,
         date_start = df_valid$date[ds], date_end = df_valid$date[de])
  }
  
  # ?????? Step 1: primary selection ??? deepest deficit endpoint ??????????????????????????????????????????????????????
  run_scores  <- sapply(dry_runs, function(ri) df_valid$cumulative_TSD[ends[ri]])
  primary_run <- dry_runs[which.min(run_scores)]
  cal         <- try_calibrate(primary_run)
  
  # Fallback A: longest run (original behaviour)
  if (is.null(cal) || nrow(cal$drought_df) < 6) {
    longest_run <- dry_runs[which.max(r$lengths[dry_runs])]
    cal_fb      <- try_calibrate(longest_run)
    if (!is.null(cal_fb) && nrow(cal_fb$drought_df) >= 3) {
      cal <- cal_fb
      message(label, ": deepest episode < 6 months; fell back to longest dry run.")
    }
  }
  
  # ?????? Step 2 (NEW): p-validity retry ?????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
  # If p ??? (0,1) the recursion diverges.  Try remaining dry runs in order of
  # decreasing deficit severity until a valid p is found.
  found <- !is.null(cal) && cal$p > 0 && cal$p < 1   # TRUE if initial cal already valid
  
  if (!found) {
    invalid_p <- if (!is.null(cal)) round(cal$p, 4) else "NA"
    message(label, ": initial calibration p = ", invalid_p,
            " is outside (0,1).  Searching alternative dry periods...")
    
    sorted_runs <- dry_runs[order(run_scores)]   # most-negative deficit first
    for (ri in sorted_runs) {
      cal_try <- try_calibrate(ri)
      if (!is.null(cal_try) && cal_try$p > 0 && cal_try$p < 1) {
        cal   <- cal_try
        found <- TRUE
        message(label, ": valid calibration found using episode ",
                format(cal$date_start), " ??? ", format(cal$date_end),
                "  (p = ", round(cal$p, 4), ")")
        break
      }
    }
    if (!found) {
      warning(label, ": no dry period yields p ??? (0,1).  ",
              "TSDI values will be numerically unstable.  ",
              "Consider adjusting DROUGHT_CLASS_C or inspecting data quality.")
    }
  }
  
  if (!found) {
    sorted_runs <- dry_runs[order(run_scores)]   # needed here too if not already defined
    best_ri  <- sorted_runs[which.min(abs(sapply(sorted_runs, function(ri) {
      ct <- try_calibrate(ri); if (is.null(ct)) Inf else abs(ct$p - 0.5)
    })))]
    best_cal <- try_calibrate(best_ri)
    stop(label, ": no dry period ??? 6 months yields p ??? (0,1). ",
         "Best available: p = ", if (!is.null(best_cal)) round(best_cal$p,4) else "NA",
         ". Reduce minimum episode length or inspect data.")
  }
  drought_df  <- cal$drought_df
  m <- cal$m; b <- cal$b; p <- cal$p; q <- cal$q
  
  cat("\n", label, " ??? dominant drought episode:\n",
      "   ", format(cal$date_start), " to ", format(cal$date_end),
      " (", nrow(drought_df), " months )\n")
  cat(label, " ??? drought monograph regression:\n",
      "   slope (m) =",     round(m, 4), "\n",
      "   intercept (b) =", round(b, 4), "\n")
  cat(label, " ??? derived coefficients:\n",
      "   p =", round(p, 6), "\n",
      "   q =", round(q, 6), "\n")
  
  if (p <= 0 || p >= 1)
    warning(label, ": p = ", round(p, 4), " is still outside (0,1).")
  
  # ?????? Recursive TSDI  (Eq. 2) ??????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
  n    <- nrow(df)
  TSDI <- rep(NA_real_, n)
  first_valid <- which(!is.na(df$TSD))[1]
  if (!is.na(first_valid))
    TSDI[first_valid] <- 0.02 * df$TSD[first_valid]
  
  # FIX (bug #2): the original index expression
  #   seq_len(n)[-seq_len(max(1L, first_valid - 1L))][-1L]
  # was meant to iterate i = (first_valid + 1):n, but when first_valid == 1
  # (i.e. the series starts with a valid TSD value -- true for Pipeline 1,
  # which begins Jan 1984 with full climatology), max(1L, first_valid-1L)
  # evaluates to 1, so seq_len(1) = 1, and the outer [-1] then ALSO strips the
  # first remaining element -- together dropping i = 2 entirely from the loop.
  # TSDI[2] was left NA, and when the loop resumed at i = 3, 'prev' silently
  # fell back to 0 instead of the true TSDI[2], injecting a wrong initial
  # condition one step into the series (the error then decays at rate p each
  # step, but it is a real bug, not just a style issue -- worse the closer p
  # is to 1). Replaced with an explicit, unambiguous sequence.
  if (!is.na(first_valid) && first_valid < n) {
    for (i in seq.int(first_valid + 1L, n)) {
      if (!is.na(df$TSD[i])) {
        prev    <- if (!is.na(TSDI[i - 1])) TSDI[i - 1] else 0
        TSDI[i] <- p * prev + q * df$TSD[i]
      } else {
        TSDI[i] <- TSDI[i - 1]
      }
    }
  }
  df[[label]] <- TSDI
  
  # ?????? PIT normalisation ??? N(0,1)  (Farahmand & AghaKouchak 2015) ?????????????????????????????????
  tsdi_valid <- df %>%
    filter(!is.na(.data[[label]])) %>%
    mutate(rank_tsdi      = rank(.data[[label]], ties.method = "average"),
           n_valid        = n(),
           u_blom         = (rank_tsdi - blom_a) / (n_valid + 1 - 2 * blom_a),
           norm_col_value = qnorm(u_blom)) %>%
    rename(!!paste0(label, "_norm") := norm_col_value) %>%
    dplyr::select(date, rank_tsdi, n_valid, u_blom, !!paste0(label, "_norm"))
  df <- left_join(df, tsdi_valid, by = "date")
  
  list(data          = df,
       p             = p, q = q, m = m, b = b,
       drought_start = cal$date_start,
       drought_end   = cal$date_end)
}

# ---- 5c. Diagnostic plots for a single pipeline ----------------------------
plot_pipeline <- function(res, label, basin_label,
                          storage_units_label = "Mm??") {
  df        <- res$data
  norm_col  <- paste0(label, "_norm")
  p_val     <- res$p
  q_val     <- res$q
  
  dband <- data.frame(xmin = res$drought_start, xmax = res$drought_end,
                      ymin = -Inf, ymax = Inf)
  
  p1 <- ggplot(df, aes(date, storage_m3)) +
    geom_rect(data = dband,
              aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
              fill = "lightblue", alpha = 0.3, inherit.aes = FALSE) +
    geom_line(linewidth = 0.8, colour = "#2c7bb6") +
    labs(title    = paste(basin_label, "??? Lake Storage"),
         subtitle = "Shaded band = dominant drought episode used for calibration",
         x = NULL, y = paste0("Storage (", storage_units_label, ")")) +
    theme_bw(base_size = 12)
  print(p1)
  
  p2 <- ggplot(df, aes(date, TSD)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_ribbon(aes(ymin = pmin(TSD, 0), ymax = 0),
                fill = "#d7191c", alpha = 0.3) +
    geom_ribbon(aes(ymin = 0, ymax = pmax(TSD, 0)),
                fill = "#1a9641", alpha = 0.3) +
    geom_line(linewidth = 0.6) +
    labs(title    = paste(basin_label, "??? Total Storage Deficit (TSD %)"),
         subtitle = "Red = deficit | Green = surplus",
         x = NULL, y = "TSD (%)") +
    theme_bw(base_size = 12)
  print(p2)
  
  thresh <- data.frame(
    yint  = c(-1, -2, -3, -4),
    label = c("Mild (-1)", "Moderate (-2)", "Severe (-3)", "Extreme (-4)")
  )
  p3 <- ggplot(df, aes(date, .data[[label]])) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_hline(data = thresh,
               aes(yintercept = yint, colour = label),
               linetype = "dotted", linewidth = 0.7) +
    geom_line(linewidth = 0.8) +
    scale_colour_manual(
      values = c("Mild (-1)"="#fee08b","Moderate (-2)"="#fc8d59",
                 "Severe (-3)"="#d73027","Extreme (-4)"="#7b0000"),
      name = "Palmer class") +
    labs(title    = paste(basin_label, "???", label, "(Palmer-style)"),
         subtitle = paste0("p = ", round(p_val, 4), "  |  q = ",
                           round(q_val, 4), "  |  C = ", DROUGHT_CLASS_C),
         x = NULL, y = label) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom")
  print(p3)
  
  df_std <- data.frame(date = df$date, value = df[[norm_col]])
  p4 <- plot_standard_basin_ts(df_std, index_name = label, color = "#2c7bb6")
  print(p4)
}

# ---- 5d. Summary console output -------------------------------------------
print_summary <- function(res, label, out_file) {
  df       <- res$data
  norm_col <- paste0(label, "_norm")
  cat("\n============================================================\n")
  cat(label, " ANALYSIS COMPLETE\n")
  cat("============================================================\n")
  cat("Time range : ", format(min(df$date)), " to ", format(max(df$date)), "\n")
  cat("Records    : ", nrow(df), "\n")
  cat("\nCalibration:\n")
  cat("  p              = ", round(res$p, 6), "\n")
  cat("  q              = ", round(res$q, 6), "\n")
  cat("  Drought period : ", format(res$drought_start),
      " to ", format(res$drought_end), "\n")
  cat("\nRaw ", label, " range    : ",
      round(min(df[[label]],    na.rm = TRUE), 4), " to ",
      round(max(df[[label]],    na.rm = TRUE), 4), "\n")
  cat("Standardised range : ",
      round(min(df[[norm_col]], na.rm = TRUE), 4), " to ",
      round(max(df[[norm_col]], na.rm = TRUE), 4), "\n")
  cat("Output CSV : ", out_file, "\n")
}

# ---- 5e. Lake regulation diagnostic  (NEW in v2) ---------------------------
# Reads the GeoLakes long file, filters to the Nechako watershed bounding box,
# assigns a regulation status to every lake, and produces:
#   ??? console table of per-lake statistics
#   ??? processed_data/Pipeline1_lake_regulation_check.csv
#   ??? Pipeline1_LakeRegulation_Diagnostics.pdf  (spatial map + time series)
#
# Regulation classification logic
# ???????????????????????????????????????????????????????????????????????????????????????????????????
# IMPOUNDED (direct)    : exact coordinate match to the known Nechako Reservoir
#                         GeoLakes entry (lat ??? 54.25??N, lon ??? 125.79??W).
# LIKELY IMPOUNDED      : centroid falls within the reservoir zone polygon
#                         (53.1???54.6??N, 127.8???125.0??W, elevation ??? 840 m).
#                         These are likely Ootsa / Whitesail / Eutsuk /
#                         Knewstubb / Natalkuz sub-basins that form the
#                         Nechako Reservoir system.
# DOWNSTREAM-REGULATED  : east of Kenney Dam on the Nechako mainstem corridor;
#                         not impounded but subject to regulated releases.
# NATURAL               : remainder of the watershed bounding box.
#
# NOTE: If the 'basin' column in the GeoLakes file contains a Nechako basin
#   label (e.g. "Nechako", "Fraser-Nechako"), the function will also print a
#   basin-filtered subset so you can cross-check the coordinate-based flags.
check_regulated_lakes <- function(
    lat_range  = NECHAKO_WS_LAT,
    lon_range  = NECHAKO_WS_LON,
    dam_lat    = KENNEY_DAM_LAT,
    dam_lon    = KENNEY_DAM_LON
) {
  cat("\n\n========================================================\n")
  cat("PIPELINE 1 ??? LAKE REGULATION DIAGNOSTIC\n")
  cat("  Checking which GeoLakes in the Nechako watershed are\n")
  cat("  affected by Kenney Dam (Nechako Reservoir, 1952)\n")
  cat("========================================================\n\n")
  
  if (!file.exists(GEOLAKES_FILE)) {
    warning("GEOLAKES_FILE not found: ", GEOLAKES_FILE,
            "\nSkipping regulation diagnostic.")
    return(invisible(NULL))
  }
  
  df_geo <- read_csv(GEOLAKES_FILE, show_col_types = FALSE)
  cat("GeoLakes file: ", nrow(df_geo), " rows |",
      n_distinct(df_geo$file_idx), " unique lakes (file_idx)\n\n")
  
  # ?????? Print unique basin / state labels to help manual cross-checking ?????????????????????
  if ("basin" %in% names(df_geo)) {
    cat("Unique 'basin' values (full file):\n")
    print(sort(unique(df_geo$basin)))
    cat("\n")
  }
  if ("state" %in% names(df_geo)) {
    cat("Unique 'state' values (full file):\n")
    print(sort(unique(df_geo$state)))
    cat("\n")
  }
  
  # ?????? Filter to watershed bounding box ???????????????????????????????????????????????????????????????????????????????????????????????????????????????
  df_ws <- df_geo %>%
    filter(between(.data[[GEOLAKES_LAT_COL]], lat_range[1], lat_range[2]),
           between(.data[[GEOLAKES_LON_COL]], lon_range[1], lon_range[2]))
  
  cat(sprintf("Lakes within bounding box [%.1f???%.1f??N, %.1f???%.1f??W]: %d rows, %d unique lakes\n\n",
              lat_range[1], lat_range[2], abs(lon_range[2]), abs(lon_range[1]),
              nrow(df_ws), n_distinct(df_ws$file_idx)))
  
  if (nrow(df_ws) == 0) {
    warning("No lakes found within the bounding box.  ",
            "Check NECHAKO_WS_LAT / NECHAKO_WS_LON settings.")
    return(invisible(NULL))
  }
  
  # ?????? Per-lake summary statistics ?????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
  lake_summary <- df_ws %>%
    group_by(file_idx, lake_name, lake_id,
             lat  = .data[[GEOLAKES_LAT_COL]],
             lon  = .data[[GEOLAKES_LON_COL]],
             lake_area_km2, elevation_m) %>%
    summarise(
      n_obs             = n(),
      n_valid_storage   = sum(!is.na(.data[[GEOLAKES_STORAGE_COL]])),
      pct_valid         = round(n_valid_storage / n_obs * 100, 1),
      date_start        = min(as.Date(date), na.rm = TRUE),
      date_end          = max(as.Date(date), na.rm = TRUE),
      mean_storage_Mm3  = mean(.data[[GEOLAKES_STORAGE_COL]], na.rm = TRUE),
      sd_storage_Mm3    = sd(  .data[[GEOLAKES_STORAGE_COL]], na.rm = TRUE),
      cv_pct            = round(
        sd(.data[[GEOLAKES_STORAGE_COL]], na.rm = TRUE) /
          mean(.data[[GEOLAKES_STORAGE_COL]], na.rm = TRUE) * 100, 1),
      range_storage_Mm3 = max(.data[[GEOLAKES_STORAGE_COL]], na.rm = TRUE) -
        min(.data[[GEOLAKES_STORAGE_COL]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_storage_Mm3))
  
  # FIX (bug #3): elevation_m is NA for every lake in this dataset (GeoLakes
  # does not populate it here). The `is.na(elevation_m) | elevation_m >=
  # elev_min_m` guard was silently treating every NA as a PASS, which means
  # the documented elevation >= 840 m criterion was never actually being
  # enforced -- in_reservoir_zone was really just a lat/lon box test. That is
  # not fixable by inferring elevation (we don't have it), so instead we make
  # the no-op explicit and loud rather than silent.
  n_no_elev <- sum(is.na(lake_summary$elevation_m))
  if (n_no_elev > 0) {
    cat(sprintf(
      "\n  NOTE: %d / %d lakes have no elevation_m value. The documented\n",
      n_no_elev, nrow(lake_summary)))
    cat(sprintf(
      "  elevation >= %g m reservoir-zone criterion CANNOT be applied to\n",
      RESERVOIR_ZONE$elev_min_m))
    cat("  these lakes and is being skipped (treated as pass-through) for\n")
    cat("  them -- in_reservoir_zone is effectively a lat/lon-only test for\n")
    cat("  these rows. Verify elevation_m is genuinely unavailable upstream\n")
    cat("  in the GeoLakes extraction if a real elevation filter is needed.\n\n")
  }
  
  # ?????? Regulation flag ?????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
  lake_summary <- lake_summary %>%
    mutate(
      # 1. Direct match: known Nechako Reservoir GeoLakes centroid
      is_nechako_reservoir = abs(lat - NECHAKO_LAT) <= NECHAKO_TOL * 3 &
        abs(lon - NECHAKO_LON) <= NECHAKO_TOL * 3,
      
      # 2. Reservoir zone: west of dam, within impoundment footprint
      # NOTE: elevation_m is missing for all/most lakes here (see console
      # diagnostic above) -- the elevation part of this test is a documented
      # pass-through (no-op) whenever elevation_m is NA, not an enforced
      # filter. See note above.
      in_reservoir_zone = between(lat, RESERVOIR_ZONE$lat[1], RESERVOIR_ZONE$lat[2]) &
        between(lon, RESERVOIR_ZONE$lon[1], RESERVOIR_ZONE$lon[2]) &
        (is.na(elevation_m) | elevation_m >= RESERVOIR_ZONE$elev_min_m),
      
      # 3. Downstream regulated: east of dam in Nechako River corridor
      downstream_of_dam = between(lat, 53.0, 54.5) &
        between(lon, -127.3, -122.5) &
        !(between(lat, RESERVOIR_ZONE$lat[1], RESERVOIR_ZONE$lat[2]) &
            between(lon, RESERVOIR_ZONE$lon[1], RESERVOIR_ZONE$lon[2])),
      
      regulation_status = case_when(
        is_nechako_reservoir             ~ "IMPOUNDED ??? Nechako Reservoir (direct match)",
        in_reservoir_zone & !is_nechako_reservoir ~ "LIKELY IMPOUNDED ??? within reservoir zone",
        downstream_of_dam                ~ "DOWNSTREAM ??? regulated releases from dam",
        TRUE                             ~ "NATURAL ??? likely unaffected"
      )
    )
  
  # ?????? Console table ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
  # FIX (bug #5): the previous print() relied on tibble's default terminal
  # width, which silently hid regulation_status off-screen ("i 1 more
  # variable: regulation_status <chr>") -- the one column this diagnostic
  # exists to report never actually showed up in the console. Fixed with
  # width = Inf plus a dedicated compact status table and a category-count
  # summary so the classification is unmissable.
  cat("?????? Per-lake summary (sorted by mean storage) ??????????????????????????????????????????????????????????????????\n")
  print(
    lake_summary %>%
      dplyr::select(file_idx, lake_name, lat, lon, lake_area_km2, elevation_m,
                    n_valid_storage, pct_valid,
                    mean_storage_Mm3, sd_storage_Mm3, cv_pct,
                    regulation_status) %>%
      mutate(across(c(lat, lon, lake_area_km2, elevation_m,
                      mean_storage_Mm3, sd_storage_Mm3), ~round(., 3))),
    n = Inf, width = Inf
  )
  
  cat("\n?????? Regulation status ??? compact view (guaranteed visible) ????????????????????????????????????????????????????\n")
  print(
    lake_summary %>%
      dplyr::select(file_idx, lake_name, lat, lon, regulation_status) %>%
      mutate(lat = round(lat, 3), lon = round(lon, 3)) %>%
      as.data.frame(),
    row.names = FALSE
  )
  
  cat("\n?????? Regulation status ??? category counts ???????????????????????????????????????????????????????????????????????????????????\n")
  print(
    lake_summary %>%
      count(regulation_status, name = "n_lakes") %>%
      arrange(desc(n_lakes)) %>%
      as.data.frame(),
    row.names = FALSE
  )
  
  # ?????? Optional basin-column cross-check ???????????????????????????????????????????????????????????????????????????????????????????????????????????????
  if ("basin" %in% names(df_ws)) {
    nechako_basins <- df_ws %>%
      filter(grepl("nechako|fraser.*nechako|nechako.*fraser",
                   basin, ignore.case = TRUE)) %>%
      dplyr::select(file_idx, lake_name, lat = .data[[GEOLAKES_LAT_COL]],
                    lon = .data[[GEOLAKES_LON_COL]], basin) %>%
      distinct()
    cat("\n?????? Lakes matched by 'basin' column (Nechako keyword) ???????????????????????????????????????\n")
    if (nrow(nechako_basins) > 0) print(nechako_basins, n = Inf) else
      cat("  No matches ??? check unique 'basin' values printed above.\n")
  }
  
  # ?????? Save CSV ??????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
  out_dir_proc <- file.path(OUT_DIR, "processed_data")
  dir.create(out_dir_proc, showWarnings = FALSE, recursive = TRUE)
  reg_csv <- file.path(out_dir_proc, "Pipeline1_lake_regulation_check.csv")
  write_csv(lake_summary, reg_csv)
  cat("\nRegulation table saved to:", reg_csv, "\n")
  
  # ?????? Plots ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
  reg_colours <- c(
    "IMPOUNDED ??? Nechako Reservoir (direct match)"  = "#d73027",
    "LIKELY IMPOUNDED ??? within reservoir zone"      = "#fc8d59",
    "DOWNSTREAM ??? regulated releases from dam"      = "#fee090",
    "NATURAL ??? likely unaffected"                   = "#4575b4"
  )
  
  pdf(file.path(OUT_DIR, "Pipeline1_LakeRegulation_Diagnostics.pdf"),
      width = 12, height = 8)
  
  # Plot A ??? spatial map
  p_map <- ggplot(lake_summary %>% filter(!is.na(mean_storage_Mm3)),
                  aes(x = lon, y = lat,
                      colour = regulation_status,
                      size   = log1p(pmax(mean_storage_Mm3, 0)))) +
    geom_point(alpha = 0.80) +
    # Kenney Dam marker
    geom_point(data = data.frame(lon = dam_lon, lat = dam_lat),
               aes(x = lon, y = lat), inherit.aes = FALSE,
               shape = 17, size = 5, colour = "black") +
    annotate("text", x = dam_lon + 0.05, y = dam_lat + 0.08,
             label = "Kenney Dam (1952)", hjust = 0, size = 3.5,
             fontface = "bold") +
    # Reservoir sub-basin labels
    geom_text(data = RESERVOIR_COMPONENTS,
              aes(x = lon, y = lat, label = name),
              inherit.aes = FALSE, size = 2.8, colour = "darkred",
              nudge_y = 0.06, hjust = 0.5) +
    scale_colour_manual(values = reg_colours) +
    scale_size_continuous(name = "log(mean storage + 1 Mm??)", range = c(2, 10)) +
    labs(
      title    = "GeoLakes in Nechako Watershed ??? Regulation Status",
      subtitle = paste0("Bounding box: ", lat_range[1], "???", lat_range[2],
                        "??N, ", abs(lon_range[2]), "???", abs(lon_range[1]),
                        "??W  |  ??? = Kenney Dam (~53.73??N, 127.22??W)"),
      x = "Longitude (??W)",  y = "Latitude (??N)",
      colour = "Regulation status"
    ) +
    coord_fixed(ratio = 1.4) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 9))
  print(p_map)
  
  # Plot B ??? storage time series faceted by lake
  lakes_with_data <- lake_summary %>%
    filter(n_valid_storage >= 5) %>% pull(file_idx)
  
  df_ts <- df_ws %>%
    filter(file_idx %in% lakes_with_data) %>%
    left_join(lake_summary %>%
                dplyr::select(file_idx, regulation_status),
              by = "file_idx") %>%
    mutate(
      lake_label = case_when(
        !is.na(lake_name) & lake_name != "" ~
          paste0(file_idx, ": ", lake_name),
        TRUE ~
          paste0("idx_", file_idx, "\n(",
                 round(.data[[GEOLAKES_LAT_COL]], 2), "??N, ",
                 round(.data[[GEOLAKES_LON_COL]], 2), "??W)")
      )
    )
  
  if (nrow(df_ts) > 0) {
    p_ts <- ggplot(df_ts,
                   aes(x = as.Date(date),
                       y = .data[[GEOLAKES_STORAGE_COL]],
                       colour = regulation_status)) +
      geom_line(linewidth = 0.5, na.rm = TRUE) +
      geom_point(size = 0.8, na.rm = TRUE) +
      facet_wrap(~ lake_label, scales = "free_y", ncol = 2) +
      scale_colour_manual(values = reg_colours) +
      scale_x_date(date_breaks = "10 years", date_labels = "%Y") +
      labs(title    = "Storage Time Series ??? GeoLakes in Nechako Watershed",
           subtitle = "Coloured by regulation status (free y-axis per lake)",
           x = NULL, y = "Storage (Mm??)", colour = NULL) +
      theme_bw(base_size = 9) +
      theme(legend.position  = "bottom",
            strip.text       = element_text(size = 7),
            axis.text.x      = element_text(angle = 45, hjust = 1))
    print(p_ts)
  }
  
  # Plot C ??? mean storage bar chart coloured by regulation status
  p_bar <- ggplot(
    lake_summary %>%
      filter(!is.na(mean_storage_Mm3)) %>%
      mutate(lake_label = paste0(
        file_idx, ": ",
        ifelse(is.na(lake_name) | lake_name == "", "(unnamed)", lake_name)
      )),
    aes(x = reorder(lake_label, mean_storage_Mm3),
        y = mean_storage_Mm3,
        fill = regulation_status)
  ) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = reg_colours) +
    labs(title = "Mean GeoLakes Storage ??? Nechako Watershed Lakes",
         x     = "Lake (file_idx: name)",
         y     = "Mean storage (Mm??)",
         fill  = "Regulation status") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 9))
  print(p_bar)
  
  dev.off()
  cat("Regulation diagnostic PDF saved to:",
      file.path(OUT_DIR, "Pipeline1_LakeRegulation_Diagnostics.pdf"), "\n")
  
  invisible(lake_summary)
}

# ============================================================================
# ?????? HELPER ??? Annual coverage assessment (Pipeline 1)
# ============================================================================
assess_annual_coverage <- function(df,
                                   min_coverage          = 0.70,
                                   require_consecutive   = 5L,
                                   expected_obs_per_year = NULL) {
  raw_counts <- df %>%
    group_by(year) %>%
    summarise(n_months_present = n(),
              n_valid          = sum(!is.na(storage_m3)),
              .groups          = "drop") %>%
    filter(year < year(Sys.Date()))
  
  if (is.null(expected_obs_per_year)) {
    expected_obs_per_year <- max(round(median(raw_counts$n_valid)), 1L)
    cat(sprintf("   expected_obs_per_year auto-detected: %d\n",
                expected_obs_per_year))
  } else {
    cat(sprintf("   expected_obs_per_year (user-supplied): %d\n",
                expected_obs_per_year))
  }
  
  cov_tbl <- raw_counts %>%
    mutate(coverage        = n_valid / expected_obs_per_year,
           meets_threshold = coverage >= min_coverage)
  
  cat("\n?????? Annual storage coverage (Pipeline 1) ??????????????????????????????????????????????????????????????????????????????\n")
  cat(sprintf("   Threshold: %.0f%%  |  Require %.0f consecutive qualifying years\n",
              min_coverage * 100, require_consecutive))
  print(
    cov_tbl %>%
      mutate(coverage_pct = sprintf("%4.0f%%", coverage * 100),
             flag         = if_else(meets_threshold, "  OK", "  SPARSE")) %>%
      dplyr::select(year, n_months_present, n_valid, coverage_pct, flag),
    n = Inf
  )
  
  years  <- cov_tbl$year
  meets  <- cov_tbl$meets_threshold
  first_stable <- NA_integer_
  for (i in seq_along(years)) {
    end_idx <- i + require_consecutive - 1L
    if (end_idx > length(years)) break
    if (all(meets[i:end_idx])) { first_stable <- years[i]; break }
  }
  if (is.na(first_stable)) {
    if (any(meets)) {
      first_stable <- years[which(meets)[1]]
      cat(sprintf("\n  NOTE: could not find %d consecutive qualifying years; ",
                  require_consecutive))
      cat(sprintf("falling back to first qualifying year: %d\n", first_stable))
    } else {
      stop("No years meet the coverage threshold.")
    }
  } else {
    cat(sprintf("\n  First year starting ???%d consecutive qualifying years: %d\n",
                require_consecutive, first_stable))
  }
  list(table = cov_tbl, first_stable_year = first_stable)
}

# ============================================================================
# ?????? HELPER ??? Coverage-jump diagnostic (run manually from console)
# ============================================================================
diagnose_coverage_jump <- function(split_year, input_file = P1_INPUT_FILE) {
  df_raw <- read_csv(input_file, show_col_types = FALSE) %>%
    mutate(date   = as.Date(date),
           period = if_else(year(date) < split_year, "pre", "post"))
  
  cat("\n========== Coverage-jump diagnostic ==========\n")
  cat("Split year:", split_year, "\n")
  
  dup_dates <- df_raw %>% count(date) %>% filter(n > 1)
  cat(sprintf("\n1. Duplicate 'date' rows: %d\n", nrow(dup_dates)))
  if (nrow(dup_dates) > 0) print(head(dup_dates, 10))
  
  if (all(c("n_lakes_contributing", "n_lakes_total") %in% names(df_raw))) {
    cat("\n2. Lake participation, pre vs post:\n")
    print(df_raw %>% group_by(period) %>%
            summarise(mean_n_contributing = mean(n_lakes_contributing, na.rm=TRUE),
                      mean_n_total        = mean(n_lakes_total, na.rm=TRUE),
                      .groups = "drop"))
  }
  if ("product" %in% names(df_raw)) {
    cat("\n3. 'product' field, pre vs post:\n")
    print(df_raw %>% count(period, product))
  }
  cat("\n4. total_storage_Mm3 magnitude, pre vs post:\n")
  print(df_raw %>% group_by(period) %>%
          summarise(mean_storage = mean(total_storage_Mm3, na.rm=TRUE),
                    sd_storage   = sd(total_storage_Mm3, na.rm=TRUE),
                    .groups = "drop"))
  
  invisible(df_raw)
}

# ============================================================================
# ?????? SECTION 6 ??????  PIPELINE 1  ???  Basin-wide GeoLakes TSDI
# ============================================================================
run_pipeline_1 <- function(min_annual_coverage    = 0.70,
                           require_consecutive   = 5L,
                           calibration_start_year = NULL,
                           quality_max            = 1L,
                           expected_obs_per_year  = NULL) {
  cat("\n\n========================================================\n")
  cat("PIPELINE 1 ??? Basin-wide TSDI  (GeoLakes aggregated)\n")
  qual_label <- if (!is.null(quality_max))
    sprintf("quality <= %d only", quality_max) else "all quality levels"
  cat(sprintf("  Input filter: %s\n", qual_label))
  cat("========================================================\n")
  
  df_raw <- read_csv(P1_INPUT_FILE, show_col_types = FALSE)
  cat("Columns found:", paste(names(df_raw), collapse = ", "), "\n")
  stopifnot("date" %in% names(df_raw), "total_storage_Mm3" %in% names(df_raw))
  
  if (!is.null(quality_max)) {
    if (!"quality" %in% names(df_raw)) {
      warning("quality_max set but no 'quality' column found; filter skipped.")
    } else {
      n_raw  <- nrow(df_raw)
      df_raw <- df_raw %>% filter(quality <= quality_max)
      n_kept <- nrow(df_raw)
      cat(sprintf("  Quality filter (???%d): kept %d / %d rows, dropped %d\n",
                  quality_max, n_kept, n_raw, n_raw - n_kept))
      if (n_kept == 0) stop("All rows removed by quality filter.")
    }
  }
  
  df_input <- df_raw %>%
    transmute(date       = as.Date(date),
              year       = year(as.Date(date)),
              month      = month(as.Date(date)),
              storage_m3 = total_storage_Mm3)
  
  cov <- assess_annual_coverage(df_input,
                                min_coverage          = min_annual_coverage,
                                require_consecutive   = require_consecutive,
                                expected_obs_per_year = expected_obs_per_year)
  
  start_yr <- if (!is.null(calibration_start_year)) {
    cat(sprintf("  Manually-set calibration start year: %d\n", calibration_start_year))
    calibration_start_year
  } else {
    cov$first_stable_year
  }
  
  df_calib <- df_input %>% filter(year >= start_yr) %>%
    dplyr::select(date, storage_m3)
  
  cat(sprintf("\n  Full record : %s to %s (%d months)\n",
              min(df_input$date), max(df_input$date), nrow(df_input)))
  cat(sprintf("  Calibration : %s to %s (%d months, %.0f%%)\n\n",
              min(df_calib$date), max(df_calib$date), nrow(df_calib),
              nrow(df_calib) / nrow(df_input) * 100))
  
  res <- compute_tsdi(df_in           = df_calib,
                      drought_class_c = DROUGHT_CLASS_C,
                      use_smoothing   = USE_SMOOTHING,
                      smooth_window   = SMOOTH_WINDOW,
                      blom_a          = BLOM_A,
                      label           = "TSDI")
  
  res$coverage_table    <- cov$table
  res$calibration_start <- start_yr
  
  out_file <- file.path(OUT_DIR, "Nechako_TSDI_Pipeline1_Output.csv")
  write_csv(res$data, out_file)
  
  dir.create(file.path(OUT_DIR, "processed_data"), showWarnings = FALSE, recursive = TRUE)
  write_csv(cov$table,
            file.path(OUT_DIR, "processed_data", "Pipeline1_annual_coverage.csv"))
  
  pdf(file.path(OUT_DIR, "Pipeline1_Diagnostics.pdf"), width = 10, height = 6)
  plot_pipeline(res, label = "TSDI",
                basin_label = "Nechako Basin (all lakes)",
                storage_units_label = "Mm??")
  dev.off()
  
  print_summary(res, label = "TSDI", out_file = out_file)
  invisible(res)
}

# ============================================================================
# ?????? SECTION 7 ??????  PIPELINE 2  ???  Nechako Reservoir NTSDI  (LOCAL CSV ONLY)
# ============================================================================
# CHANGE v2: GEOLAKES branch removed.  Pipeline 2 now reads only from the
# local authority CSV file.  NECHAKO_DATA_SOURCE is fixed to "LOCAL_DATA" in
# Section 3 above; the run_pipeline_2() function below no longer checks or
# branches on the data-source setting.
run_pipeline_2 <- function() {
  cat("\n\n========================================================\n")
  cat("PIPELINE 2 ??? Nechako Reservoir NTSDI\n")
  cat("  Data source: LOCAL_DATA (", LOCAL_DATA_TYPE, ")\n")
  cat("========================================================\n")
  
  cat("Reading local authority CSV: ", LOCAL_DATA_FILE, "\n")
  
  if (LOCAL_DATA_TYPE == "STORAGE") {
    df_local <- read_csv(LOCAL_DATA_FILE, show_col_types = FALSE)
    stopifnot(LOCAL_DATE_COL %in% names(df_local),
              LOCAL_VALUE_COL %in% names(df_local))
    df_input <- df_local %>%
      transmute(date       = as.Date(.data[[LOCAL_DATE_COL]],
                                     format = LOCAL_DATE_FORMAT),
                storage_m3 = .data[[LOCAL_VALUE_COL]])
    
  } else if (LOCAL_DATA_TYPE == "LEVEL") {
    df_local <- read_csv(LOCAL_DATA_FILE, show_col_types = FALSE)
    stopifnot(LOCAL_DATE_COL %in% names(df_local),
              LOCAL_VALUE_COL %in% names(df_local))
    df_input <- df_local %>%
      transmute(date        = as.Date(.data[[LOCAL_DATE_COL]],
                                      format = LOCAL_DATE_FORMAT),
                level_m_asl = .data[[LOCAL_VALUE_COL]]) %>%
      mutate(storage_m3 = level_to_storage(level_m_asl))
    
  } else if (LOCAL_DATA_TYPE == "LEVEL_DAILY") {
    cat("  Parsing Year + Day-of-Year format and converting to storage...\n")
    # FIX (bug #1): the file was previously assumed to be whitespace-delimited
    # with exactly 3 fields (Year, DOY, Level_m). In practice every row in the
    # live file triggered a readr parsing failure ("expected 3 columns, actual
    # 4 columns"), meaning the fixed 3-column spec below was WRONG and was
    # silently mis-parsing every row (readr drops/misaligns the unexpected
    # trailing field instead of erroring out).
    #
    # Fix: detect the actual field count from the file itself (via
    # count.fields()) instead of hard-coding it, print a diagnostic so the
    # extra column is visible/verifiable, and only then parse. If a 4th field
    # is present we assume it is a trailing quality/status flag (common in
    # BC Hydro / BCWIS daily level exports) and drop it after showing its
    # unique values -- but we no longer do this silently.
    # Read raw lines and strip leading/trailing whitespace first. The file
    # has trailing whitespace after Level_m on every row, which was making
    # read_table()'s positional column-guesser infer a phantom 4th (empty)
    # column even though count.fields() correctly saw 3. Trimming removes
    # the ambiguity so both diagnosis and parsing agree.
    lines_raw   <- readLines(LOCAL_DATA_FILE)
    lines_clean <- trimws(lines_raw)
    lines_clean <- lines_clean[nzchar(lines_clean)]   # drop blank lines
    
    n_fields_per_line <- utils::count.fields(textConnection(lines_clean), sep = "")
    field_tbl <- table(n_fields_per_line)
    cat("  Detected field counts across rows of", LOCAL_DATA_FILE, "(post-trim):\n")
    print(field_tbl)
    n_fields <- as.integer(names(field_tbl)[which.max(field_tbl)])
    cat("  Using n_fields =", n_fields,
        "(most common field count; used to build column spec)\n")
    
    if (any(n_fields_per_line != n_fields)) {
      n_bad <- sum(n_fields_per_line != n_fields)
      warning(n_bad, " row(s) in ", LOCAL_DATA_FILE,
              " do not match the dominant field count (", n_fields,
              "). These rows may parse incorrectly -- inspect them directly.")
    }
    
    if (n_fields == 3) {
      col_names_ld <- c("Year", "DOY", "Level_m")
      col_classes_ld <- c("integer", "integer", "numeric")
    } else if (n_fields == 4) {
      col_names_ld <- c("Year", "DOY", "Level_m", "Extra_field")
      col_classes_ld <- c("integer", "integer", "numeric", "character")
    } else {
      stop("LOCAL_DATA_FILE has an unexpected number of whitespace-delimited ",
           "fields (", n_fields, "). Expected 3 (Year, DOY, Level_m) or 4 ",
           "(Year, DOY, Level_m, extra flag/column). Inspect the file before ",
           "proceeding.")
    }
    
    df_raw_ld <- read.table(
      text        = lines_clean,
      header      = FALSE,
      col.names   = col_names_ld,
      colClasses  = col_classes_ld,
      strip.white = TRUE
    ) %>% as_tibble()
    
    # Verification gate: parsed row count must match cleaned line count.
    if (nrow(df_raw_ld) != length(lines_clean)) {
      stop("Parsed row count (", nrow(df_raw_ld), ") does not match cleaned ",
           "line count (", length(lines_clean), "). Parsing is misaligned -- ",
           "inspect the file before proceeding.")
    }
    
    if (n_fields == 4) {
      cat("  NOTE: file has a 4th whitespace-delimited field beyond",
          "Year/DOY/Level_m. Unique values of that field (first 20 shown):\n")
      print(utils::head(sort(unique(df_raw_ld$Extra_field)), 20))
      cat("  This 4th field is being DROPPED from the analysis below.\n",
          "  If it is a quality/status flag (e.g. estimated vs. measured),\n",
          "  consider filtering on it before proceeding.\n")
    }
    
    
    df_input <- df_raw_ld %>%
      filter(Year <= year(Sys.Date())) %>%
      mutate(date       = as.Date(paste(Year, DOY), format = "%Y %j"),
             storage_m3 = level_to_storage(Level_m)) %>%
      dplyr::select(date, storage_m3)
    
  } else {
    stop("LOCAL_DATA_TYPE must be 'STORAGE', 'LEVEL', or 'LEVEL_DAILY'.")
  }
  
  n_na <- sum(is.na(df_input$storage_m3))
  if (n_na > 0)
    cat("  WARNING:", n_na,
        "records produced NA storage (level outside hypsometric table).\n")
  
  n_before <- nrow(df_input)
  df_input <- df_input %>% filter(!is.na(storage_m3))
  n_after  <- nrow(df_input)
  if (n_before != n_after)
    cat("  Removed", n_before - n_after, "rows with NA storage.\n")
  
  cat(sprintf("  Valid storage records: %d  (%s to %s)\n",
              n_after, min(df_input$date), max(df_input$date)))
  df_input <- df_input %>% filter(date <= as_date("2026-01-01"))
  # Aggregate daily ??? monthly (TSDI algorithm expects monthly input)
  df_in <- df_input %>%
    mutate(year  = year(date), month = month(date)) %>%
    group_by(year, month) %>%
    summarise(
      storage_m3 = mean(storage_m3, na.rm = TRUE),
      date       = as.Date(paste(year[1], month[1], "01"), "%Y %m %d"),
      .groups    = "drop"
    ) %>%
    arrange(date) %>%
    dplyr::select(date, storage_m3)
  # ====== Run TSDI algorithm ==============================
  res <- compute_tsdi(df_in           = df_in, # <-- FIXED
                      drought_class_c = DROUGHT_CLASS_C,
                      use_smoothing   = USE_SMOOTHING,
                      smooth_window   = SMOOTH_WINDOW,
                      blom_a          = BLOM_A,
                      label           = "NTSDI")
  
  res$data <- res$data %>%
    mutate(data_source    = "LOCAL_DATA",
           local_data_type = LOCAL_DATA_TYPE)
  
  out_file <- file.path(OUT_DIR, "Nechako_NTSDI_Pipeline2_Output.csv")
  write_csv(res$data, out_file)
  
  # Storage units depend on the conversion path
  p2_units <- if (LOCAL_DATA_TYPE %in% c("LEVEL", "LEVEL_DAILY")) "m??" else "Mm??"
  
  pdf(file.path(OUT_DIR, "Pipeline2_Diagnostics.pdf"), width = 10, height = 6)
  plot_pipeline(
    res,
    label       = "NTSDI",
    basin_label = paste0("Nechako Reservoir [Local CSV: ", LOCAL_DATA_TYPE, "]"),
    storage_units_label = p2_units
  )
  dev.off()
  
  print_summary(res, label = "NTSDI", out_file = out_file)
  invisible(res)
}

# ============================================================================
# ?????? SECTION 8 ??????  COMPARISON PLOT  (only when PIPELINE_MODE == "BOTH")
# ============================================================================
plot_comparison <- function(res1, res2) {
  cat("\n\n== Generating comparison plot ==\n")
  comp <- inner_join(
    res1$data %>% dplyr::select(date, TSDI, TSDI_norm),
    res2$data %>% dplyr::select(date, NTSDI, NTSDI_norm),
    by = "date"
  )
  if (nrow(comp) == 0) {
    warning("No overlapping dates between pipelines.  Skipping comparison plot.")
    return(invisible(NULL))
  }
  cat("Overlapping records:", nrow(comp), "\n")
  
  comp_long_raw <- comp %>%
    pivot_longer(c(TSDI, NTSDI), names_to = "Index", values_to = "value") %>%
    mutate(Index = recode(Index,
                          TSDI  = "Pipeline 1: Basin-wide TSDI",
                          NTSDI = "Pipeline 2: Nechako Reservoir NTSDI"))
  
  pc1 <- ggplot(comp_long_raw, aes(date, value, colour = Index)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_line(linewidth = 0.8, alpha = 0.85) +
    scale_colour_manual(
      values = c("Pipeline 1: Basin-wide TSDI"         = "#2c7bb6",
                 "Pipeline 2: Nechako Reservoir NTSDI" = "#d7191c")) +
    labs(title    = "Comparison ??? Basin TSDI vs. Nechako Reservoir NTSDI",
         subtitle = paste0("Overlap: ", format(min(comp$date)), " ??? ",
                           format(max(comp$date)),
                           "  |  P2 data: LOCAL_DATA (", LOCAL_DATA_TYPE, ")"),
         x = NULL, y = "Index value", colour = NULL) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom")
  print(pc1)
  
  comp_long_norm <- comp %>%
    pivot_longer(c(TSDI_norm, NTSDI_norm),
                 names_to = "Index", values_to = "value") %>%
    mutate(Index = recode(Index,
                          TSDI_norm  = "Pipeline 1: Basin TSDI_norm",
                          NTSDI_norm = "Pipeline 2: Nechako Res. NTSDI_norm"))
  
  pc2 <- ggplot(comp_long_norm, aes(date, value, colour = Index)) +
    geom_hline(yintercept =  0,    linetype = "dashed",  colour = "grey50") +
    geom_hline(yintercept = -1.0,  linetype = "dotted",  colour = "#fc8d59") +
    geom_hline(yintercept = -1.5,  linetype = "dotted",  colour = "#d73027") +
    geom_line(linewidth = 0.8, alpha = 0.85) +
    scale_colour_manual(
      values = c("Pipeline 1: Basin TSDI_norm"          = "#2c7bb6",
                 "Pipeline 2: Nechako Res. NTSDI_norm"  = "#d7191c")) +
    labs(title    = "Comparison ??? Standardised TSDI_norm vs. NTSDI_norm",
         subtitle = "Dotted lines: moderate (???1.0) and severe (???1.5) thresholds",
         x = NULL, y = "Standardised index [N(0,1)]", colour = NULL) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom")
  print(pc2)
  
  corr_val <- round(cor(comp$TSDI_norm, comp$NTSDI_norm, use = "complete.obs"), 3)
  pc3 <- ggplot(comp, aes(TSDI_norm, NTSDI_norm)) +
    geom_point(alpha = 0.4, size = 1.5, colour = "#404040") +
    geom_smooth(method = "lm", se = TRUE, colour = "#2c7bb6", linewidth = 0.9) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4, size = 4,
             label = paste0("r = ", corr_val)) +
    labs(title    = "Scatter ??? Basin TSDI_norm vs. Nechako Reservoir NTSDI_norm",
         subtitle = "Dashed = 1:1 line; blue = OLS fit",
         x = "Basin TSDI_norm (Pipeline 1)",
         y = "Nechako Reservoir NTSDI_norm (Pipeline 2)") +
    theme_bw(base_size = 12)
  print(pc3)
  
  cat("Pearson r (TSDI_norm vs NTSDI_norm):", corr_val, "\n")
}

# ============================================================================
# ?????? SECTION 9 ??????  MAIN EXECUTION
# ============================================================================
res_p1 <- NULL
res_p2 <- NULL

if (PIPELINE_MODE %in% c("BASIN", "BOTH")) {
  res_p1 <- run_pipeline_1(quality_max = 1L)
  
  # ?????? NEW in v2: run lake regulation diagnostic after Pipeline 1 ???????????????????????????????????????
  # This checks which GeoLakes lakes in the Nechako watershed may be affected
  # by regulation from Kenney Dam (Nechako Reservoir, 1952).
  res_reg <- check_regulated_lakes()
  
  # Standalone basin-mean time series PNG (TSDI_norm), same onset/shading
  # convention and visual style as the companion SSI/SWEI/SSPI-1 scripts.
  plot_index_timeseries(
    dates_vec   = res_p1$data$date,
    values_vec  = res_p1$data$TSDI_norm,
    index_label = "TSDI",
    title_label = "TSDI (Basin-Wide)",
    out_dir     = OUT_DIR,
    onset       = -0.5
  )
}

if (PIPELINE_MODE %in% c("NECHAKO_RESERVOIR", "BOTH")) {
  res_p2 <- run_pipeline_2()
  
  # Standalone basin-mean time series PNG (NTSDI_norm)
  plot_index_timeseries(
    dates_vec   = res_p2$data$date,
    values_vec  = res_p2$data$NTSDI_norm,
    index_label = "NTSDI",
    title_label = "NTSDI (Nechako Reservoir)",
    out_dir     = OUT_DIR,
    onset       = -0.5
  )
}

if (PIPELINE_MODE == "BOTH" && !is.null(res_p1) && !is.null(res_p2)) {
  pdf(file.path(OUT_DIR, "Comparison_Plots.pdf"), width = 10, height = 6)
  plot_comparison(res_p1, res_p2)
  dev.off()
}

sink()
cat("\n\nDone. All outputs saved to:", normalizePath(OUT_DIR), "\n")
cat("Basin-mean time series plots saved to:", file.path(normalizePath(OUT_DIR), "figures"), "\n")