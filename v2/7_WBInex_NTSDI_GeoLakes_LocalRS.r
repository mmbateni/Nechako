# ============================================================================
# NECHAKO BASIN - TOTAL STORAGE DEFICIT INDEX  (Extended Dual-Pipeline)
# VERSION 3 - patched from 7_WBInex_NTSDI_GeoLakes_LocalRS.r
#
# CHANGES IN V3
#   1. compute_tsdi() now ALWAYS computes two indices:
#        (a) p/q recursion TSDI (Awange et al. 2016) + Blom normalisation
#        (b) TSD-only Blom normalisation (no memory term, no calibration)
#      and reports whichever one is valid:
#        index_reported / index_method columns =
#          "PQ_RECURSION"                  strict-gap calibration succeeded
#          "PQ_RECURSION_RELAXED_GAP_<d>d" calibration only worked after letting
#                                          dry runs span unobserved winter months
#          "TSD_BLOM"                      no valid p/q -> TSD-only fallback
#      A p/q result is accepted only if p is in (0,1), q is finite, the raw
#      recursion is finite, and cor(TSD, raw TSDI) >= PQ_MIN_COR.
#      Nothing is forced; both series are always written to the output CSV.
#   2. Pipeline 1: seasonal-sampling diagnostic (observed calendar months),
#      relaxed-gap calibration schedule, and a fallback instead of a hard failure.
#   3. Pipeline 2: NTSDI_START_YEAR drops the reservoir-filling years from
#      climatology/calibration/ranking; out-of-range levels are reported by year
#      and side (below/above the storage curve) and clamped + flagged by default.
#   4. Regulation diagnostic no longer depends on lake_name (all NA in this
#      dataset); coordinate-only screens are labelled as screens, and each lake
#      is compared with the reservoir storage-curve maximum.
#   5. Date-scale warnings fixed (annotate() rects use Date +/-Inf).
#   6. Comparison plot uses TSD (%) and the reported index, so it works whichever
#      method each pipeline ended up with.
#
# Based on:
#   Awange et al. (2016) Adv. Water Resources 94, 45-59
#   Yirdaw et al. (2008) J. Hydrology 356, 84-92
#   Farahmand & AghaKouchak (2015) Adv. Water Resources 76, 140-145
# ============================================================================
# --- SHARED FRAMEWORK UTILITIES ---
# Needed for the framework-wide drought-event definition (DROUGHT_ONSET,
# DROUGHT_END, EVENT_MIN_DURATION_BY_SCALE, identify_drought_events_ts()).
rm(list = ls())

# ---- Working directory (framework rule) -------------------------------------
# Run every script from the project root: set the working directory there
# BEFORE running (setwd(), or open the RStudio project). All paths below are
# derived from getwd(); no path is hard-coded.
WD_PATH    <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
UTILS_FILE <- file.path(WD_PATH, "DROUGHT_ANALYSIS_utils.R")
if (!file.exists(UTILS_FILE))
  stop("DROUGHT_ANALYSIS_utils.R not found in the working directory: ", WD_PATH,
       "\nSet the working directory to the project root before running this script.")

source(UTILS_FILE)

# TSDI / NTSDI are monthly storage indices with no accumulation window, so
# their events use the scale-1 minimum duration of EVENT_MIN_DURATION_BY_SCALE.
TSDI_EVENT_SCALE <- 1L

# --- OUTPUT DIRECTORY SETUP ---
OUT_DIR <- "ntsdi_results"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

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

# Date-class infinities for annotate("rect") on date axes (avoids the
# "A <numeric> value was passed to a Date scale" warnings). structure() is used
# instead of as.Date(-Inf) so it also works on R < 4.3.
DATE_NEG_INF <- structure(-Inf, class = "Date")
DATE_POS_INF <- structure( Inf, class = "Date")

# ============================================================================
# HELPER - standard time-series plot (used inside the multi-panel diagnostic PDF)
# ============================================================================
plot_standard_basin_ts <- function(df, index_name, color = "steelblue",
                                   title_suffix = NULL,
                                   event_scale = TSDI_EVENT_SCALE) {
  if (!all(c("date", "value") %in% names(df)))
    stop("df must have 'date' and 'value' columns.")
  df <- df[order(df$date), ]
  if (!any(is.finite(df$value))) return(invisible(NULL))
  
  # Framework-wide event definition (shared with Scripts 6, 11, 12, 14):
  # detected on the full monthly calendar, so missing months end an event and
  # runs are never joined across unobserved (e.g. winter) months.
  ev      <- identify_drought_events_ts(df$date, df$value, scale = event_scale)
  n_ev    <- nrow(ev)
  min_dur <- attr(ev, "min_duration")
  df$in_event <- attr(ev, "in_event")
  df <- df[is.finite(df$value), ]
  
  y_range <- max(df$value) - min(df$value)
  y_pad   <- y_range * 0.10
  y_lo    <- min(-3.5, min(df$value) - y_pad)
  y_hi    <- max( 3.5, max(df$value) + y_pad)
  
  band <- function(ymin, ymax, fill)
    ggplot2::annotate("rect", xmin = DATE_NEG_INF, xmax = DATE_POS_INF,
                      ymin = ymin, ymax = ymax, fill = fill, alpha = 0.15)
  drought_bands <- list(
    band(-Inf, -2.0, "#7B0025"), band(-2.0, -1.5, "#D73027"),
    band(-1.5, -1.0, "#F46D43"), band(-1.0,  1.0, "#2E7D32"),
    band( 1.0,  1.5, "#90EE90"), band( 1.5,  2.0, "#66BB6A"),
    band( 2.0,  Inf, "#4575B4")
  )
  
  ggplot2::ggplot(df, ggplot2::aes(x = date, y = value)) +
    drought_bands +
    ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = c(-1, 1), linetype = "dashed",
                        color = "gray50", linewidth = 0.3) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = ifelse(in_event, value, 0), ymax = 0),
      fill = "#c0392b", alpha = 0.35) +
    ggplot2::geom_line(color = color, linewidth = 0.6) +
    ggplot2::scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
    ggplot2::coord_cartesian(ylim = c(y_lo, y_hi)) +
    ggplot2::labs(
      title    = sprintf("%s Time Series%s", toupper(index_name),
                         if (!is.null(title_suffix)) paste0("  [", title_suffix, "]") else ""),
      subtitle = sprintf("%d drought events (onset < %.2f, end >= %.2f, min %d months)",
                         n_ev, DROUGHT_ONSET, DROUGHT_END, min_dur),
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
# HELPER - standalone basin time series PNG (shared across both pipelines)
# ============================================================================
# Same event definition as every other framework script (DROUGHT_ANALYSIS_utils.R):
# onset < DROUGHT_ONSET, termination >= DROUGHT_END, scale-dependent minimum
# duration. Shading covers event months only, so it matches the event count.
plot_index_timeseries <- function(dates_vec, values_vec, index_label,
                                  title_label, out_dir,
                                  method_label = NULL,
                                  event_scale = TSDI_EVENT_SCALE) {
  df <- data.frame(date = as.Date(dates_vec), value = values_vec)
  df <- df[order(df$date), ]
  if (!any(is.finite(df$value))) {
    cat(sprintf("  [PLOT] Skipping %s time series: no finite values\n", index_label))
    return(invisible(NULL))
  }
  
  # Detect on the full monthly calendar BEFORE dropping missing values, so an
  # unobserved month ends an event instead of silently joining two runs.
  ev       <- identify_drought_events_ts(df$date, df$value, scale = event_scale)
  n_events <- nrow(ev)
  min_dur  <- attr(ev, "min_duration")
  df$below <- attr(ev, "in_event")
  df <- df[is.finite(df$value), ]
  month_id <- as.integer(format(df$date, "%Y")) * 12L + as.integer(format(df$date, "%m"))
  df$grp   <- cumsum(c(1L, diff(as.integer(df$below)) != 0L | diff(month_id) != 1L))
  
  y_lo <- min(-2.2, floor(min(df$value) * 10) / 10)
  y_hi <- max( 2.2, ceiling(max(df$value) * 10) / 10)
  
  bnd <- function(ymin, ymax, fill, alpha)
    ggplot2::annotate("rect", xmin = DATE_NEG_INF, xmax = DATE_POS_INF,
                      ymin = ymin, ymax = ymax, fill = fill, alpha = alpha)
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = date, y = value)) +
    bnd( 1,   Inf, "#4472C4", 0.10) +
    bnd( 0,   1,   "#70AD47", 0.08) +
    bnd(-1,   0,   "#70AD47", 0.05) +
    bnd(-Inf, -1,  "#C00000", 0.08) +
    ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
    ggplot2::geom_hline(yintercept = c(-1, 1), color = "grey30",
                        linetype = "dashed", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = c(-1.5, 1.5), color = "grey55",
                        linetype = "dotted", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = c(-2, 2), color = "grey55",
                        linetype = "dotted", linewidth = 0.4) +
    ggplot2::geom_ribbon(data = subset(df, below),
                         ggplot2::aes(ymin = value, ymax = 0, group = grp),
                         fill = "#C0504D", alpha = 0.45) +
    ggplot2::geom_line(color = "#1F77B4", linewidth = 0.45) +
    ggplot2::scale_x_date(date_breaks = "5 years", date_labels = "%Y",
                          expand = ggplot2::expansion(mult = c(0.01, 0.015))) +
    ggplot2::scale_y_continuous(breaks = seq(-4, 4, 1)) +
    ggplot2::coord_cartesian(ylim = c(y_lo, y_hi)) +
    ggplot2::labs(
      title    = sprintf("%s Time Series", title_label),
      subtitle = sprintf("%d drought events  (onset < %.2f, end >= %.2f, min %d months)%s",
                         n_events, DROUGHT_ONSET, DROUGHT_END, min_dur,
                         if (!is.null(method_label)) paste0("  |  method: ", method_label) else ""),
      x = "Year", y = sprintf("%s Index Value", index_label)
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(color = "grey85"),
      panel.grid.major.y = ggplot2::element_blank(),
      plot.title         = ggplot2::element_text(face = "bold", size = 15),
      plot.subtitle      = ggplot2::element_text(color = "grey40", size = 11),
      axis.text.x        = ggplot2::element_text(angle = 30, hjust = 1),
      plot.margin        = ggplot2::margin(10, 16, 10, 10)
    )
  
  figures_dir <- file.path(out_dir, "figures")
  if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)
  fname <- file.path(figures_dir, sprintf("%s_basin_timeseries.png", tolower(index_label)))
  ggplot2::ggsave(fname, p, width = 11, height = 5.5, dpi = 200, bg = "white")
  cat(sprintf("  [PLOT] %s time series saved: %s  (%d events, method %s)\n",
              index_label, basename(fname), n_events,
              if (is.null(method_label)) "n/a" else method_label))
  invisible(p)
}

# ============================================================================
# -- SECTION 1 --  PIPELINE SELECTOR
# ============================================================================
VALID_PIPELINE_MODES <- c("BASIN", "NECHAKO_RESERVOIR", "BOTH")

# ---- To use in batch/scheduled jobs, set DEFAULT_PIPELINE_MODE below ----
DEFAULT_PIPELINE_MODE <- "BOTH"

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
  
  # Priority 1: command-line argument (e.g. Rscript script.r BOTH)
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1 && !is.na(normalize_pipeline_answer(args[1]))) {
    mode <- normalize_pipeline_answer(args[1])
    cat(sprintf("-> Pipeline set from command-line argument: %s\n", mode))
    return(mode)
  }
  
  # Priority 2: RStudio GUI dialog
  if (interactive() && requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    raw_input <- rstudioapi::showPrompt(
      title   = "Pipeline Selection",
      message = paste0(menu_text, "\nEnter BASIN / NECHAKO_RESERVOIR / BOTH (or 1/2/3):"),
      default = DEFAULT_PIPELINE_MODE
    )
    if (is.null(raw_input)) {
      cat(sprintf("-> Dialog cancelled: using default pipeline %s\n", DEFAULT_PIPELINE_MODE))
      return(DEFAULT_PIPELINE_MODE)
    }
    mode <- normalize_pipeline_answer(raw_input)
    while (is.na(mode)) {
      raw_input <- rstudioapi::showPrompt(
        title   = "Invalid input - try again",
        message = "Please enter BASIN, NECHAKO_RESERVOIR, BOTH, or 1/2/3:",
        default = DEFAULT_PIPELINE_MODE
      )
      if (is.null(raw_input)) {
        cat(sprintf("-> Dialog cancelled: using default pipeline %s\n", DEFAULT_PIPELINE_MODE))
        return(DEFAULT_PIPELINE_MODE)
      }
      mode <- normalize_pipeline_answer(raw_input)
    }
    cat(sprintf("-> Pipeline set from RStudio dialog: %s\n", mode))
    return(mode)
  }
  
  # Priority 3: plain readline() fallback
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
  
  # Priority 4: non-interactive fallback
  cat("\n", menu_text, sep = "")
  cat("Enter BASIN / NECHAKO_RESERVOIR / BOTH (or 1/2/3): ")
  ans  <- tryCatch(readLines(con = "stdin", n = 1), error = function(e) character(0))
  mode <- if (length(ans) >= 1) normalize_pipeline_answer(ans[1]) else NA_character_
  if (is.na(mode)) {
    cat(sprintf("-> No valid input received: using default pipeline %s\n", DEFAULT_PIPELINE_MODE))
    cat("  (Edit DEFAULT_PIPELINE_MODE at the top of the script to change this.)\n")
    mode <- DEFAULT_PIPELINE_MODE
  }
  mode
}

PIPELINE_MODE <- prompt_pipeline_mode()
cat("\nSelected PIPELINE_MODE:", PIPELINE_MODE, "\n\n")

# ============================================================================
# -- SECTION 2 --  PIPELINE 1  SETTINGS  (basin-wide GeoLakes)
# ============================================================================
P1_INPUT_FILE <- "Lakes/Glolakes/nechako_extracted/BEST_ESTIMATE_basin_total_Mm3.csv"
# Required columns: date, total_storage_Mm3

# Interpolation across short internal gaps only (monthly grid, in months).
# Observations in this record fall in an open-water season (about Apr-Sep), so
# the unobserved winter gap (about 6 months) is deliberately NOT interpolated.
P1_INTERP_MAXGAP <- 3L

# Dry-run gap schedule for p/q calibration (days between consecutive valid
# observations that still counts as the same dry run):
#   75  = strict   : runs cannot span the winter gap
#   250 = relaxed  : runs may span Sep -> following May (about 242 days).
#                    The relaxed result credits TSD to unobserved months and is
#                    labelled PQ_RECURSION_RELAXED_GAP_250d in the output.
P1_GAP_SCHEDULE_DAYS <- c(75, 250)

# ============================================================================
# -- SECTION 3 --  PIPELINE 2  SETTINGS  (Nechako Reservoir - LOCAL FILE ONLY)
# ============================================================================
NECHAKO_DATA_SOURCE <- "LOCAL_DATA"   # FIXED - do not change

# NOTE: despite the .csv/.txt naming, this file is whitespace-delimited
# (Year, DOY, Level_m - 3 fields, no header).
LOCAL_DATA_FILE    <- file.path("Lakes", "dailylevel.txt")

#   "STORAGE"      - file already contains storage in Mm^3
#   "LEVEL"        - file has a WSE column (m ASL); storage via hypsometric LUT
#   "LEVEL_DAILY"  - file has Year, Day-of-Year, Level_m (no header row)
LOCAL_DATA_TYPE   <- "LEVEL_DAILY"

LOCAL_VALUE_COL   <- "Level_m"
LOCAL_DATE_COL    <- "Date"
LOCAL_DATE_FORMAT <- "%Y-%m-%d"

# First year used for Pipeline 2 climatology, calibration and ranking.
# The record starts in 1955 while the reservoir was still filling after Kenney
# Dam closure (1952); the first 144 days are below the storage curve and the
# first calibrated "drought" (1955-06 to 1956-10) is a filling artefact.
# 1960 is a PLACEHOLDER - confirm the start of normal operations with the
# operator (BC Hydro / Rio Tinto) and edit here.
NTSDI_START_YEAR <- 1960L

# Levels outside the storage-curve range:
#   "CLAMP" - clamp to the curve's min/max level before conversion and flag the
#             affected days (n_clamped_days column). Below-range clamping makes
#             deficits slightly milder than reality, so check the printed table.
#   "DROP"  - set storage to NA (old behaviour; removes the extreme days).
HYPSO_OUT_OF_RANGE <- "CLAMP"

# P2 uses a single strict gap rule (continuous daily data aggregated monthly).
P2_GAP_SCHEDULE_DAYS <- c(75)

# GeoLakes long-format file (regulation diagnostic only)
GEOLAKES_FILE        <- "Lakes/Glolakes/nechako_extracted/LandsatPlusICESat2_long.csv"
GEOLAKES_DATE_COL    <- "date"
GEOLAKES_STORAGE_COL <- "storage_Mm3"
GEOLAKES_LAT_COL     <- "lat"
GEOLAKES_LON_COL     <- "lon"

# Previously proposed "Nechako Reservoir" GeoLakes centroid. UNVERIFIED: it sits
# near Burns Lake, not in the main Ootsa/Whitesail body. Used only to flag lakes.
UNVERIFIED_CENTROID_LAT <- 54.25
UNVERIFIED_CENTROID_LON <- -125.79
NECHAKO_TOL <- 0.01
RESERVOIR_NAME_PATTERN <- "nechako reservoir|ootsa|whitesail|eutsuk|knewstubb|natalkuz|tahtsa|tetachuck"

# Storage curve (stage-volume look-up table for LEVEL / LEVEL_DAILY)
HYPSO_TABLE <- NULL
load_hypso_table <- function(path = file.path("Lakes", "Storage curve.xlsx")) {
  if (!file.exists(path)) stop("Storage curve workbook not found: ", path)
  Storage_curve <- readxl::read_excel(path)
  col_elev <- grep("Elevation", names(Storage_curve), value = TRUE, ignore.case = TRUE)[1]
  col_stor <- grep("Storage", names(Storage_curve), value = TRUE, ignore.case = TRUE)[1]
  if (is.na(col_elev) || is.na(col_stor)) {
    stop("Could not identify Elevation/Storage columns in ", path, ". Found: ",
         paste(names(Storage_curve), collapse = ", "))
  }
  # Inherited assumptions: elevation is feet and storage is million cubic metres.
  elevation_raw <- as.numeric(Storage_curve[[col_elev]])
  storage_raw   <- as.numeric(Storage_curve[[col_stor]])
  med_elev    <- stats::median(elevation_raw, na.rm = TRUE)
  med_storage <- stats::median(storage_raw,   na.rm = TRUE)
  if (!is.finite(med_elev) || med_elev < 1000) {
    stop("Storage-curve elevation values do not look like feet (median = ",
         signif(med_elev, 5), "). Verify the workbook units before conversion.")
  }
  if (!is.finite(med_storage) || med_storage > 1e6) {
    stop("Storage-curve values do not look like million cubic metres (median = ",
         signif(med_storage, 5), "). Verify the workbook units before conversion.")
  }
  out <- data.frame(level_m = elevation_raw * 0.3048,
                    storage_m3 = storage_raw * 1e6) %>%
    filter(is.finite(level_m), is.finite(storage_m3)) %>% arrange(level_m)
  if (nrow(out) < 2L || anyDuplicated(out$level_m))
    stop("Storage curve must contain at least two unique, finite elevation values.")
  out
}

# ============================================================================
# -- SECTION 3.5 --  KENNEY DAM / REGULATION CHECK SETTINGS
# ============================================================================
KENNEY_DAM_LAT <- 53.566667
KENNEY_DAM_LON <- -124.95

# Nechako watershed coarse bounding box for the lake regulation scan
NECHAKO_WS_LAT  <- c(52.5, 55.5)
NECHAKO_WS_LON  <- c(-128.7, -122.0)

# Screening envelope for the Nechako Reservoir (NOT a hydrological polygon).
RESERVOIR_ZONE <- list(
  lat = c(53.1, 54.3),
  lon = c(-128.7, -124.9),
  elev_min_m = 850.4        # 2,790 ft
)

RESERVOIR_COMPONENTS <- data.frame(
  name = c("Ootsa Lake / Nechako Res.", "Whitesail Lake",
           "Eutsuk Lake", "Knewstubb Lake", "Natalkuz Lake", "Tahtsa Lake (west extent)"),
  lat  = c(53.60, 53.55, 53.35, 53.68, 53.78, 53.60),
  lon  = c(-126.00, -126.95, -127.25, -126.30, -126.55, -128.67)
)

# ============================================================================
# -- SECTION 4 --  SHARED ALGORITHM SETTINGS
# ============================================================================
DROUGHT_CLASS_C <- -3      # Palmer class: -1 mild | -2 moderate | -3 severe | -4 extreme
USE_SMOOTHING   <- FALSE
SMOOTH_WINDOW   <- 3       # months
BLOM_A          <- 0.44    # Blom plotting-position constant

P_VALID_RANGE   <- c(0.2, 1) # p must lie strictly inside this interval while was c(0, 1); 0.2 is a judgment call, not from the paperer
MIN_RUN_OBS     <- 6L      # minimum observations in a dry run used for calibration
PQ_MIN_COR      <- 0.5     # recursion rejected if cor(TSD, raw index) is below this

# ============================================================================
# -- SECTION 5 --  HELPER FUNCTIONS
# ============================================================================

# ---- 5a. Level -> storage (piecewise linear) --------------------------------
# clamp = TRUE : levels outside the curve are clamped to its min/max level
#                (caller reports how many days were affected).
# clamp = FALSE: levels outside the curve return NA with a warning.
level_to_storage <- function(level_vec, hypso = HYPSO_TABLE, clamp = FALSE) {
  stopifnot(
    is.data.frame(hypso),
    "level_m"    %in% names(hypso),
    "storage_m3" %in% names(hypso)
  )
  hypso <- hypso[order(hypso$level_m), ]
  lo <- min(hypso$level_m); hi <- max(hypso$level_m)
  if (clamp) {
    level_vec <- pmin(pmax(level_vec, lo), hi)
  } else {
    out_of_range <- level_vec < lo | level_vec > hi
    if (any(out_of_range, na.rm = TRUE)) {
      warning(sum(out_of_range, na.rm = TRUE),
              " level value(s) outside hypsometric table range [",
              lo, ", ", hi, " m].")
    }
  }
  approx(x = hypso$level_m, y = hypso$storage_m3,
         xout = level_vec, method = "linear", rule = 1)$y
}

# ---- 5b. Blom normalisation of a vector (NA preserved) -----------------------
blom_normalize <- function(x, a = BLOM_A) {
  out <- rep(NA_real_, length(x))
  ok  <- is.finite(x)
  n   <- sum(ok)
  if (n >= 3L)
    out[ok] <- qnorm((rank(x[ok], ties.method = "average") - a) / (n + 1 - 2 * a))
  out
}

# ---- 5c. Recursion  TSDI[i] = p^dt * TSDI[prev] + q * TSD[i]  ---------------
# dt = elapsed whole months since the previous valid observation (>= 1), so
# unobserved months (e.g. winter) decay the state instead of being skipped.
run_recursion <- function(df, p, q) {
  out   <- rep(NA_real_, nrow(df))
  valid <- which(is.finite(df$TSD))
  if (!length(valid)) return(out)
  state <- 0.02 * df$TSD[valid[1L]]
  out[valid[1L]] <- state
  last <- valid[1L]
  for (i in valid[-1L]) {
    dt <- max(1, round(as.numeric(df$date[i] - df$date[last]) / 30.4375))
    state  <- (p ^ dt) * state + q * df$TSD[i]
    out[i] <- state
    last   <- i
  }
  out
}

# ---- 5d. p/q calibration for ONE gap rule -----------------------------------
# Returns list(cal, diag, reason). cal is NULL when no eligible dry run yields
# p inside P_VALID_RANGE.
calibrate_pq <- function(df, gap_days, drought_class_c, min_run_obs, p_range) {
  valid_idx <- which(is.finite(df$TSD))
  if (length(valid_idx) < min_run_obs)
    return(list(cal = NULL, diag = NULL,
                reason = sprintf("only %d finite TSD observations", length(valid_idx))))
  
  valid_dates <- df$date[valid_idx]
  dry_flag    <- df$TSD[valid_idx] < 0
  gap_v       <- c(Inf, as.numeric(diff(valid_dates)))
  # A run continues only across adjacent valid negative observations separated
  # by <= gap_days. Non-dry observations and longer gaps split it.
  run_break      <- !dry_flag | gap_v > gap_days
  run_break[1L]  <- TRUE
  run_id         <- cumsum(run_break)
  dry_runs <- split(valid_idx[dry_flag], run_id[dry_flag])
  dry_runs <- Filter(function(idx) length(idx) > 0L, dry_runs)
  if (!length(dry_runs))
    return(list(cal = NULL, diag = NULL, reason = "no negative-TSD observations"))
  
  fit_run <- function(run_idx) {
    d <- df[run_idx, , drop = FALSE]
    tt <- as.numeric(d$date - min(d$date)) / 30.4375 + 1
    step_months <- c(1, diff(tt))
    d$drought_time  <- tt
    d$step_months   <- step_months
    d$local_cum_TSD <- cumsum(d$TSD * step_months)
    fit <- tryCatch(lm(local_cum_TSD ~ drought_time, data = d), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    cf <- coef(fit)
    if (length(cf) < 2L || any(!is.finite(cf))) return(NULL)
    m_val <- unname(cf[2L]); b_val <- unname(cf[1L])
    denom <- m_val + b_val
    if (!is.finite(denom) || abs(denom) < 1e-6) return(NULL)
    p_val <- 1 - (m_val / denom)
    q_val <- drought_class_c / denom
    if (!is.finite(p_val) || !is.finite(q_val)) return(NULL)
    list(drought_df = d, m = m_val, b = b_val, denom = denom,
         p = p_val, q = q_val,
         date_start = min(d$date), date_end = max(d$date))
  }
  
  # Diagnostic table for every dry run (same integration as the fit)
  diag_tbl <- dplyr::bind_rows(lapply(names(dry_runs), function(nm) {
    idx <- dry_runs[[nm]]
    out <- tibble(run = nm, n_obs = length(idx),
                  start = min(df$date[idx]), end = max(df$date[idx]),
                  sum_TSD = sum(df$TSD[idx]),
                  m = NA_real_, b = NA_real_, p = NA_real_)
    if (length(idx) >= 3L) {
      f <- fit_run(idx)
      if (!is.null(f)) { out$m <- f$m; out$b <- f$b; out$p <- f$p }
    }
    out
  })) %>%
    mutate(eligible = n_obs >= min_run_obs,
           valid_p  = is.finite(p) & p > p_range[1] & p < p_range[2]) %>%
    arrange(desc(n_obs))
  
  run_scores <- vapply(dry_runs, function(idx) sum(df$TSD[idx]), numeric(1))
  n_eligible <- sum(lengths(dry_runs) >= min_run_obs)
  
  cal <- NULL
  if (n_eligible > 0L) {
    for (ri in order(run_scores)) {                 # most severe episode first
      idx <- dry_runs[[ri]]
      if (length(idx) < min_run_obs) next
      cal_try <- fit_run(idx)
      if (!is.null(cal_try) && cal_try$p > p_range[1] && cal_try$p < p_range[2]) {
        cal <- cal_try
        break
      }
    }
  }
  reason <- if (!is.null(cal)) NA_character_ else if (n_eligible == 0L)
    sprintf("no dry run with >= %d observations (longest = %d)",
            min_run_obs, max(lengths(dry_runs)))
  else sprintf("%d eligible dry run(s), none gave p in (%g,%g)",
               n_eligible, p_range[1], p_range[2])
  list(cal = cal, diag = diag_tbl, reason = reason)
}

# ---- 5e. Core TSDI / NTSDI computation --------------------------------------
# Always computes TSD and the TSD-only Blom index (TSD_norm). Then tries to
# calibrate p/q (one attempt per entry of gap_schedule_days). If calibration
# succeeds AND the recursion passes the sanity checks, the p/q recursion + Blom
# index is reported; otherwise TSD_norm is reported. No coefficients are forced.
compute_tsdi <- function(df_in,
                         drought_class_c   = DROUGHT_CLASS_C,
                         use_smoothing     = USE_SMOOTHING,
                         smooth_window     = SMOOTH_WINDOW,
                         blom_a            = BLOM_A,
                         gap_schedule_days = c(75),
                         label             = "TSDI") {
  stopifnot("date" %in% names(df_in), "storage_m3" %in% names(df_in))
  
  df <- df_in %>%
    mutate(date  = as.Date(date),
           year  = year(date),
           month = month(date)) %>%
    arrange(date)
  
  if (use_smoothing) {
    df <- df %>%
      mutate(storage_m3 = zoo::rollmean(storage_m3, k = smooth_window,
                                        fill = NA, align = "center"))
    cat(label, ": applied", smooth_window, "-month rolling mean.\n")
  }
  
  # Monthly climatology (Awange et al. 2016, Eq. 1). Calendar months with no
  # observations (seasonal records) yield NA climatology and NA TSD.
  monthly_stats <- df %>%
    group_by(month) %>%
    summarise(
      clim_mean = suppressWarnings(mean(storage_m3, na.rm = TRUE)),
      clim_max  = suppressWarnings(max( storage_m3, na.rm = TRUE)),
      clim_min  = suppressWarnings(min( storage_m3, na.rm = TRUE)),
      .groups   = "drop"
    ) %>%
    mutate(
      clim_mean = ifelse(is.nan(clim_mean) | is.infinite(clim_mean), NA_real_, clim_mean),
      clim_max  = ifelse(is.infinite(clim_max), NA_real_, clim_max),
      clim_min  = ifelse(is.infinite(clim_min), NA_real_, clim_min)
    )
  df <- left_join(df, monthly_stats, by = "month")
  
  # TSD (%)  (guard against zero climatological range)
  df <- df %>%
    mutate(clim_range = clim_max - clim_min,
           TSD = ifelse(is.finite(storage_m3) & is.finite(clim_range) & clim_range > 0,
                        (storage_m3 - clim_mean) / clim_range * 100, NA_real_),
           cumulative_TSD = cumsum(replace_na(TSD, 0)))
  
  n_tsd <- sum(is.finite(df$TSD))
  if (n_tsd < 3L)
    stop(label, ": only ", n_tsd, " finite TSD observations; cannot compute any index.")
  
  # ---- (b) TSD-only Blom index: always available ----------------------------
  df$TSD_norm <- blom_normalize(df$TSD, blom_a)
  
  # ---- (a) p/q calibration over the gap schedule ----------------------------
  cal <- NULL; gap_used <- NA_real_; fail_notes <- character(0)
  for (g in gap_schedule_days) {
    cr <- calibrate_pq(df, gap_days = g, drought_class_c = drought_class_c,
                       min_run_obs = MIN_RUN_OBS, p_range = P_VALID_RANGE)
    if (!is.null(cr$diag) && nrow(cr$diag)) {
      cat("\n", label, ": dry-run calibration diagnostics, gap <= ", g,
          " d (up to 25 longest runs):\n", sep = "")
      print(utils::head(cr$diag, 25L), n = 25L, width = Inf)
      diag_file <- file.path(OUT_DIR,
                             sprintf("%s_calibration_run_diagnostics_gap%g.csv", label, g))
      tryCatch(readr::write_csv(cr$diag, diag_file), error = function(e)
        warning("Could not write calibration diagnostics: ", conditionMessage(e), call. = FALSE))
      cat("Calibration-run diagnostics saved to:", diag_file, "\n")
    }
    if (!is.null(cr$cal)) { cal <- cr$cal; gap_used <- g; break }
    note <- sprintf("gap <= %g d: %s", g, cr$reason)
    fail_notes <- c(fail_notes, note)
    message(label, ": p/q calibration unavailable (", note, ").")
  }
  
  # ---- sanity-check the recursion -------------------------------------------
  recursion_ok <- FALSE
  tsd_cor <- NA_real_
  TSDI_raw <- rep(NA_real_, nrow(df))
  if (!is.null(cal)) {
    TSDI_raw <- run_recursion(df, cal$p, cal$q)
    ok_vals  <- TSDI_raw[is.finite(df$TSD)]
    tsd_cor  <- suppressWarnings(stats::cor(df$TSD, TSDI_raw, use = "complete.obs"))
    if (!all(is.finite(ok_vals))) {
      fail_notes <- c(fail_notes, "recursion produced non-finite values")
    } else if (!is.finite(tsd_cor) || tsd_cor < PQ_MIN_COR) {
      fail_notes <- c(fail_notes,
                      sprintf("recursion rejected: cor(TSD, raw index) = %s < %g",
                              format(round(tsd_cor, 3)), PQ_MIN_COR))
    } else {
      recursion_ok <- TRUE
    }
    if (!recursion_ok) message(label, ": ", utils::tail(fail_notes, 1L))
  }
  
  df[[label]]                   <- if (recursion_ok) TSDI_raw else NA_real_
  df[[paste0(label, "_norm")]]  <- if (recursion_ok) blom_normalize(TSDI_raw, blom_a) else NA_real_
  
  # ---- choose the reported index --------------------------------------------
  if (recursion_ok) {
    method <- if (isTRUE(gap_used == gap_schedule_days[1L])) "PQ_RECURSION"
    else sprintf("PQ_RECURSION_RELAXED_GAP_%gd", gap_used)
    df$index_reported <- df[[paste0(label, "_norm")]]
    cat("\n", label, " - dominant drought episode:\n",
        "   ", format(cal$date_start), " to ", format(cal$date_end),
        " (", nrow(cal$drought_df), " observations )\n", sep = "")
    cat(label, " - drought monograph regression:  slope (m) =", round(cal$m, 4),
        " intercept (b) =", round(cal$b, 4), "\n")
    cat(label, " - derived coefficients:  p =", round(cal$p, 6),
        " q =", round(cal$q, 6), "\n")
    cat(label, " - cor(TSD, raw index) =", round(tsd_cor, 3), "\n")
    if (cal$p < 0.2)
      warning(label, ": p = ", round(cal$p, 3),
              " implies almost no storage memory; treat the recursion with caution.",
              call. = FALSE)
    if (!identical(method, "PQ_RECURSION"))
      warning(label, ": calibration only succeeded with the relaxed gap rule (",
              gap_used, " d); dry runs span unobserved months.", call. = FALSE)
  } else {
    method <- "TSD_BLOM"
    df$index_reported <- df$TSD_norm
    cat("\n", label, " - p/q recursion NOT available; reporting TSD-only Blom index.\n", sep = "")
    for (n in fail_notes) cat("   - ", n, "\n", sep = "")
  }
  df$index_method <- method
  
  list(data          = df,
       p             = if (recursion_ok) cal$p else NA_real_,
       q             = if (recursion_ok) cal$q else NA_real_,
       m             = if (recursion_ok) cal$m else NA_real_,
       b             = if (recursion_ok) cal$b else NA_real_,
       drought_start = if (recursion_ok) cal$date_start else as.Date(NA),
       drought_end   = if (recursion_ok) cal$date_end   else as.Date(NA),
       method        = method,
       gap_days_used = if (recursion_ok) gap_used else NA_real_,
       recursion_ok  = recursion_ok,
       tsd_cor       = tsd_cor,
       failure_notes = fail_notes)
}

# ---- 5f. Diagnostic plots for a single pipeline ----------------------------
plot_pipeline <- function(res, label, basin_label,
                          storage_units_label = "Mm^3") {
  df       <- res$data
  norm_col <- paste0(label, "_norm")
  
  p1 <- ggplot(df, aes(date, storage_m3))
  if (isTRUE(res$recursion_ok)) {
    dband <- data.frame(xmin = res$drought_start, xmax = res$drought_end,
                        ymin = -Inf, ymax = Inf)
    p1 <- p1 + geom_rect(data = dband,
                         aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                         fill = "lightblue", alpha = 0.3, inherit.aes = FALSE)
  }
  p1 <- p1 +
    geom_line(linewidth = 0.8, colour = "#2c7bb6", na.rm = TRUE) +
    geom_point(size = 0.7, colour = "#2c7bb6", na.rm = TRUE) +
    labs(title    = paste(basin_label, "- Lake Storage"),
         subtitle = if (isTRUE(res$recursion_ok))
           "Shaded band = dominant drought episode used for calibration"
         else "No valid p/q calibration episode (TSD-only fallback)",
         x = NULL, y = paste0("Storage (", storage_units_label, ")")) +
    theme_bw(base_size = 12)
  print(p1)
  
  p2 <- ggplot(df, aes(date, TSD)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_ribbon(aes(ymin = pmin(TSD, 0), ymax = 0), fill = "#d7191c", alpha = 0.3, na.rm = TRUE) +
    geom_ribbon(aes(ymin = 0, ymax = pmax(TSD, 0)), fill = "#1a9641", alpha = 0.3, na.rm = TRUE) +
    geom_line(linewidth = 0.6, na.rm = TRUE) +
    labs(title    = paste(basin_label, "- Total Storage Deficit (TSD %)"),
         subtitle = "Red = deficit | Green = surplus",
         x = NULL, y = "TSD (%)") +
    theme_bw(base_size = 12)
  print(p2)
  
  if (isTRUE(res$recursion_ok)) {
    thresh <- data.frame(
      yint  = c(-1, -2, -3, -4),
      label = c("Mild (-1)", "Moderate (-2)", "Severe (-3)", "Extreme (-4)")
    )
    p3 <- ggplot(df, aes(date, .data[[label]])) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_hline(data = thresh, aes(yintercept = yint, colour = label),
                 linetype = "dotted", linewidth = 0.7) +
      geom_line(linewidth = 0.8, na.rm = TRUE) +
      scale_colour_manual(
        values = c("Mild (-1)" = "#fee08b", "Moderate (-2)" = "#fc8d59",
                   "Severe (-3)" = "#d73027", "Extreme (-4)" = "#7b0000"),
        name = "Palmer class") +
      labs(title    = paste(basin_label, "-", label, "(Palmer-style, raw p/q recursion)"),
           subtitle = paste0("p = ", round(res$p, 4), "  |  q = ", round(res$q, 4),
                             "  |  C = ", DROUGHT_CLASS_C, "  |  ", res$method),
           x = NULL, y = label) +
      theme_bw(base_size = 12) + theme(legend.position = "bottom")
    print(p3)
    
    # Both standardised indices on one panel
    cmp <- data.frame(
      date   = rep(df$date, 2),
      value  = c(df[[norm_col]], df$TSD_norm),
      Method = rep(c("p/q recursion (Blom)", "TSD only (Blom)"), each = nrow(df)))
    r_cmp <- suppressWarnings(cor(df[[norm_col]], df$TSD_norm, use = "complete.obs"))
    p5 <- ggplot(cmp, aes(date, value, colour = Method)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_line(linewidth = 0.6, alpha = 0.85, na.rm = TRUE) +
      scale_colour_manual(values = c("p/q recursion (Blom)" = "#d7191c",
                                     "TSD only (Blom)"      = "#2c7bb6")) +
      labs(title = paste(basin_label, "- both standardised indices"),
           subtitle = paste0("Pearson r between them = ", round(r_cmp, 3)),
           x = NULL, y = "Standardised value [N(0,1)]", colour = NULL) +
      theme_bw(base_size = 12) + theme(legend.position = "bottom")
    print(p5)
  }
  
  df_std <- data.frame(date = df$date, value = df$index_reported)
  p4 <- plot_standard_basin_ts(df_std, index_name = label, color = "#2c7bb6",
                               title_suffix = res$method)
  if (!is.null(p4)) print(p4)
}

# ---- 5g. Summary console output --------------------------------------------
print_summary <- function(res, label, out_file) {
  df       <- res$data
  norm_col <- paste0(label, "_norm")
  rng <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) "n/a" else paste(round(min(x), 4), "to", round(max(x), 4))
  }
  cat("\n============================================================\n")
  cat(label, " ANALYSIS COMPLETE\n")
  cat("============================================================\n")
  cat("Time range : ", format(min(df$date)), " to ", format(max(df$date)), "\n")
  cat("Records    : ", nrow(df), "  (finite TSD: ", sum(is.finite(df$TSD)), ")\n", sep = "")
  cat("\nREPORTED METHOD : ", res$method, "\n", sep = "")
  if (isTRUE(res$recursion_ok)) {
    cat("\nCalibration (p/q recursion):\n")
    cat("  p              = ", round(res$p, 6), "\n")
    cat("  q              = ", round(res$q, 6), "\n")
    cat("  Drought period : ", format(res$drought_start), " to ", format(res$drought_end), "\n")
    cat("  Gap rule used  : <= ", res$gap_days_used, " days\n", sep = "")
    cat("  cor(TSD, raw)  : ", round(res$tsd_cor, 3), "\n", sep = "")
    cat("  Raw ", label, " range      : ", rng(df[[label]]), "\n", sep = "")
    cat("  ", label, "_norm range      : ", rng(df[[norm_col]]), "\n", sep = "")
  } else {
    cat("\np/q recursion: NOT AVAILABLE. Reason(s):\n")
    for (n in res$failure_notes) cat("  - ", n, "\n", sep = "")
  }
  cat("TSD_norm range (TSD-only Blom) : ", rng(df$TSD_norm), "\n", sep = "")
  cat("Reported index range           : ", rng(df$index_reported), "\n", sep = "")
  cat("Output CSV : ", out_file, "\n")
}

# ---- 5h. Lake regulation diagnostic -----------------------------------------
# Screening diagnostic only: NOT an authoritative basin-boundary assignment and
# NOT used to compute Pipeline 1. lake_name is all-NA in this dataset, so the
# status labels below describe COORDINATE screens, not confirmed regulation.
#   IMPOUNDED                : name contains "nechako reservoir"
#   LIKELY IMPOUNDED         : name matches a reservoir component AND is inside
#                              the screening envelope
#   CANDIDATE (unnamed)      : inside the screening envelope, no name available
#   DOWNSTREAM SCREEN        : east of Kenney Dam, outside the envelope
#   OUTSIDE SCREENING ZONES  : everything else
# Each lake's mean storage is also compared with the reservoir storage-curve
# maximum; a lake holding a tiny fraction of that cannot be the reservoir body.
check_regulated_lakes <- function(
    lat_range  = NECHAKO_WS_LAT,
    lon_range  = NECHAKO_WS_LON,
    dam_lat    = KENNEY_DAM_LAT,
    dam_lon    = KENNEY_DAM_LON
) {
  cat("\n\n========================================================\n")
  cat("PIPELINE 1 - LAKE REGULATION DIAGNOSTIC (coordinate screen)\n")
  cat("========================================================\n\n")
  
  if (!file.exists(GEOLAKES_FILE)) {
    warning("GEOLAKES_FILE not found: ", GEOLAKES_FILE, "\nSkipping regulation diagnostic.")
    return(invisible(NULL))
  }
  
  df_geo <- read_csv(GEOLAKES_FILE, show_col_types = FALSE)
  for (cc in c("lake_name", "lake_id", "lake_area_km2", "elevation_m"))
    if (!cc %in% names(df_geo)) df_geo[[cc]] <- NA
  cat("GeoLakes file: ", nrow(df_geo), " rows |",
      n_distinct(df_geo$file_idx), " unique lakes (file_idx)\n\n")
  
  if ("basin" %in% names(df_geo)) {
    cat("Unique 'basin' values (full file):\n"); print(sort(unique(df_geo$basin))); cat("\n")
    cat("NOTE: 'basin' is a parent-basin label (Fraser); it cannot identify the Nechako sub-basin.\n\n")
  }
  if ("state" %in% names(df_geo)) {
    cat("Unique 'state' values (full file):\n"); print(sort(unique(df_geo$state))); cat("\n")
  }
  
  df_ws <- df_geo %>%
    filter(between(.data[[GEOLAKES_LAT_COL]], lat_range[1], lat_range[2]),
           between(.data[[GEOLAKES_LON_COL]], lon_range[1], lon_range[2]))
  
  cat(sprintf("Lakes within bounding box [%.1f-%.1fN, %.1f-%.1fW]: %d rows, %d unique lakes\n\n",
              lat_range[1], lat_range[2], abs(lon_range[2]), abs(lon_range[1]),
              nrow(df_ws), n_distinct(df_ws$file_idx)))
  if (nrow(df_ws) == 0) {
    warning("No lakes found within the bounding box. Check NECHAKO_WS_LAT / NECHAKO_WS_LON.")
    return(invisible(NULL))
  }
  
  lake_summary <- df_ws %>%
    group_by(file_idx, lake_name, lake_id,
             lat = .data[[GEOLAKES_LAT_COL]],
             lon = .data[[GEOLAKES_LON_COL]],
             lake_area_km2, elevation_m) %>%
    summarise(
      n_obs             = n(),
      n_valid_storage   = sum(!is.na(.data[[GEOLAKES_STORAGE_COL]])),
      pct_valid         = round(n_valid_storage / n_obs * 100, 1),
      date_start        = min(as.Date(date), na.rm = TRUE),
      date_end          = max(as.Date(date), na.rm = TRUE),
      mean_storage_Mm3  = mean(.data[[GEOLAKES_STORAGE_COL]], na.rm = TRUE),
      sd_storage_Mm3    = sd(  .data[[GEOLAKES_STORAGE_COL]], na.rm = TRUE),
      cv_pct            = round(sd_storage_Mm3 / mean_storage_Mm3 * 100, 1),
      range_storage_Mm3 = max(.data[[GEOLAKES_STORAGE_COL]], na.rm = TRUE) -
        min(.data[[GEOLAKES_STORAGE_COL]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_storage_Mm3))
  
  n_no_name <- sum(is.na(lake_summary$lake_name) | !nzchar(as.character(lake_summary$lake_name)))
  if (n_no_name > 0)
    cat(sprintf("  NOTE: %d / %d lakes have no lake_name; name-based regulation matching is impossible for them.\n",
                n_no_name, nrow(lake_summary)))
  n_no_elev <- sum(is.na(lake_summary$elevation_m))
  if (n_no_elev > 0)
    cat(sprintf("  NOTE: %d / %d lakes have no elevation_m; the %g m elevation criterion is NOT applied.\n",
                n_no_elev, nrow(lake_summary), RESERVOIR_ZONE$elev_min_m))
  
  # Reservoir storage-curve maximum for scale comparison (optional)
  res_max_Mm3 <- tryCatch(max(load_hypso_table()$storage_m3) / 1e6,
                          error = function(e) NA_real_)
  if (is.finite(res_max_Mm3)) {
    cat(sprintf("  Reservoir storage-curve maximum: %.0f Mm3 | largest GeoLake mean: %.1f Mm3 (%.2f%%)\n\n",
                res_max_Mm3, max(lake_summary$mean_storage_Mm3, na.rm = TRUE),
                max(lake_summary$mean_storage_Mm3, na.rm = TRUE) / res_max_Mm3 * 100))
  } else {
    cat("  (Storage curve not available; skipping scale comparison.)\n\n")
  }
  
  lake_summary <- lake_summary %>%
    mutate(
      lake_name_key = stringr::str_to_lower(dplyr::coalesce(as.character(lake_name), "")),
      name_available = nzchar(lake_name_key),
      is_nechako_reservoir = name_available & stringr::str_detect(lake_name_key, "nechako reservoir"),
      named_reservoir_component = name_available &
        stringr::str_detect(lake_name_key, RESERVOIR_NAME_PATTERN),
      in_reservoir_envelope = between(lat, RESERVOIR_ZONE$lat[1], RESERVOIR_ZONE$lat[2]) &
        between(lon, RESERVOIR_ZONE$lon[1], RESERVOIR_ZONE$lon[2]),
      downstream_of_dam = between(lat, 53.0, 54.5) & lon > dam_lon & lon <= -122.5,
      matches_unverified_centroid =
        abs(lat - UNVERIFIED_CENTROID_LAT) <= NECHAKO_TOL &
        abs(lon - UNVERIFIED_CENTROID_LON) <= NECHAKO_TOL,
      mean_storage_pct_of_reservoir_max =
        if (is.finite(res_max_Mm3)) round(mean_storage_Mm3 / res_max_Mm3 * 100, 3) else NA_real_,
      regulation_status = case_when(
        is_nechako_reservoir ~ "IMPOUNDED - named Nechako Reservoir",
        named_reservoir_component & in_reservoir_envelope ~ "LIKELY IMPOUNDED - named reservoir component",
        in_reservoir_envelope ~ "CANDIDATE - inside reservoir envelope (unnamed; verify)",
        downstream_of_dam ~ "DOWNSTREAM SCREEN - verify river-corridor location",
        TRUE ~ "OUTSIDE SCREENING ZONES - not verified"
      )
    )
  
  cat("-- Per-lake summary (sorted by mean storage) ----------------------\n")
  print(
    lake_summary %>%
      dplyr::select(file_idx, lat, lon, n_valid_storage, pct_valid,
                    mean_storage_Mm3, sd_storage_Mm3, cv_pct,
                    mean_storage_pct_of_reservoir_max, matches_unverified_centroid,
                    regulation_status) %>%
      mutate(across(c(lat, lon, mean_storage_Mm3, sd_storage_Mm3), ~round(., 3))),
    n = Inf, width = Inf
  )
  
  cat("\n-- Regulation status - category counts -----------------------------\n")
  print(lake_summary %>% count(regulation_status, name = "n_lakes") %>%
          arrange(desc(n_lakes)) %>% as.data.frame(), row.names = FALSE)
  
  if (any(lake_summary$matches_unverified_centroid)) {
    cat("\n  NOTE: lake(s) ",
        paste(lake_summary$file_idx[lake_summary$matches_unverified_centroid], collapse = ", "),
        " match the previously proposed reservoir centroid. That point lies near Burns Lake,\n",
        "  not in the main Ootsa/Whitesail body; it is NOT treated as the reservoir.\n", sep = "")
  }
  
  out_dir_proc <- file.path(OUT_DIR, "processed_data")
  dir.create(out_dir_proc, showWarnings = FALSE, recursive = TRUE)
  reg_csv <- file.path(out_dir_proc, "Pipeline1_lake_regulation_check.csv")
  write_csv(lake_summary, reg_csv)
  cat("\nRegulation table saved to:", reg_csv, "\n")
  
  reg_colours <- c(
    "IMPOUNDED - named Nechako Reservoir"                    = "#d73027",
    "LIKELY IMPOUNDED - named reservoir component"           = "#fc8d59",
    "CANDIDATE - inside reservoir envelope (unnamed; verify)" = "#fdae61",
    "DOWNSTREAM SCREEN - verify river-corridor location"     = "#4575b4",
    "OUTSIDE SCREENING ZONES - not verified"                 = "#999999"
  )
  
  pdf(file.path(OUT_DIR, "Pipeline1_LakeRegulation_Diagnostics.pdf"), width = 12, height = 8)
  on.exit(if (grDevices::dev.cur() > 1L) grDevices::dev.off(), add = TRUE)
  
  p_map <- ggplot(lake_summary %>% filter(!is.na(mean_storage_Mm3)),
                  aes(x = lon, y = lat, colour = regulation_status,
                      size = log1p(pmax(mean_storage_Mm3, 0)))) +
    geom_point(alpha = 0.80) +
    geom_point(data = data.frame(lon = dam_lon, lat = dam_lat),
               aes(x = lon, y = lat), inherit.aes = FALSE,
               shape = 17, size = 5, colour = "black") +
    annotate("text", x = dam_lon + 0.05, y = dam_lat + 0.08,
             label = "Kenney Dam (1952)", hjust = 0, size = 3.5, fontface = "bold") +
    geom_text(data = RESERVOIR_COMPONENTS, aes(x = lon, y = lat, label = name),
              inherit.aes = FALSE, size = 2.8, colour = "darkred",
              nudge_y = 0.06, hjust = 0.5) +
    scale_colour_manual(values = reg_colours) +
    scale_size_continuous(name = "log(mean storage + 1 Mm^3)", range = c(2, 10)) +
    labs(title    = "GeoLakes in Nechako Watershed - coordinate-screen status",
         subtitle = paste0("Bounding box: ", lat_range[1], "-", lat_range[2], "N, ",
                           abs(lon_range[2]), "-", abs(lon_range[1]),
                           "W  |  black triangle = Kenney Dam (~53.57N, 124.95W)"),
         x = "Longitude (W)", y = "Latitude (N)", colour = "Screen status") +
    coord_fixed(ratio = 1.4) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom", legend.text = element_text(size = 9))
  print(p_map)
  
  lakes_with_data <- lake_summary %>% filter(n_valid_storage >= 5) %>% pull(file_idx)
  df_ts <- df_ws %>%
    filter(file_idx %in% lakes_with_data) %>%
    left_join(lake_summary %>% dplyr::select(file_idx, regulation_status), by = "file_idx") %>%
    mutate(lake_label = paste0("idx_", file_idx, "\n(",
                               round(.data[[GEOLAKES_LAT_COL]], 2), "N, ",
                               round(.data[[GEOLAKES_LON_COL]], 2), "W)"))
  if (nrow(df_ts) > 0) {
    p_ts <- ggplot(df_ts, aes(x = as.Date(date), y = .data[[GEOLAKES_STORAGE_COL]],
                              colour = regulation_status)) +
      geom_line(linewidth = 0.5, na.rm = TRUE) +
      geom_point(size = 0.8, na.rm = TRUE) +
      facet_wrap(~ lake_label, scales = "free_y", ncol = 2) +
      scale_colour_manual(values = reg_colours) +
      scale_x_date(date_breaks = "10 years", date_labels = "%Y") +
      labs(title = "Storage Time Series - GeoLakes in Nechako Watershed",
           subtitle = "Coloured by coordinate-screen status (free y-axis per lake)",
           x = NULL, y = "Storage (Mm^3)", colour = NULL) +
      theme_bw(base_size = 9) +
      theme(legend.position = "bottom", strip.text = element_text(size = 7),
            axis.text.x = element_text(angle = 45, hjust = 1))
    print(p_ts)
  }
  
  p_bar <- ggplot(
    lake_summary %>% filter(!is.na(mean_storage_Mm3)) %>%
      mutate(lake_label = paste0("idx_", file_idx)),
    aes(x = reorder(lake_label, mean_storage_Mm3), y = mean_storage_Mm3,
        fill = regulation_status)) +
    geom_col() + coord_flip() +
    scale_fill_manual(values = reg_colours) +
    labs(title = "Mean GeoLakes Storage - Nechako Watershed Lakes",
         x = "Lake (file_idx)", y = "Mean storage (Mm^3)", fill = "Screen status") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", legend.text = element_text(size = 9))
  print(p_bar)
  
  dev.off()
  cat("Regulation diagnostic PDF saved to:",
      file.path(OUT_DIR, "Pipeline1_LakeRegulation_Diagnostics.pdf"), "\n")
  invisible(lake_summary)
}

# ============================================================================
# -- HELPER - Annual coverage assessment (Pipeline 1)
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
    # Coverage is relative to the MEDIAN observed count per year, so this table
    # measures year-to-year consistency, not completeness against 12 months.
    expected_obs_per_year <- max(round(stats::median(raw_counts$n_valid)), 1L)
    cat(sprintf("   expected_obs_per_year auto-detected: %d (relative to median; not 12)\n",
                expected_obs_per_year))
  } else {
    cat(sprintf("   expected_obs_per_year (user-supplied): %d\n", expected_obs_per_year))
  }
  
  cov_tbl <- raw_counts %>%
    mutate(coverage        = n_valid / expected_obs_per_year,
           meets_threshold = coverage >= min_coverage)
  
  cat("\n-- Annual storage coverage (Pipeline 1) --------------------------\n")
  cat(sprintf("   Threshold: %.0f%%  |  Require %.0f consecutive qualifying years\n",
              min_coverage * 100, require_consecutive))
  print(
    cov_tbl %>%
      mutate(coverage_pct = sprintf("%4.0f%%", coverage * 100),
             flag         = if_else(meets_threshold, "  OK", "  SPARSE")) %>%
      dplyr::select(year, n_months_present, n_valid, coverage_pct, flag),
    n = Inf
  )
  
  years <- cov_tbl$year
  meets <- cov_tbl$meets_threshold
  first_stable <- NA_integer_
  for (i in seq_along(years)) {
    end_idx <- i + require_consecutive - 1L
    if (end_idx > length(years)) break
    if (all(meets[i:end_idx])) { first_stable <- years[i]; break }
  }
  if (is.na(first_stable)) {
    if (any(meets)) {
      first_stable <- years[which(meets)[1]]
      cat(sprintf("\n  NOTE: could not find %d consecutive qualifying years; ", require_consecutive))
      cat(sprintf("falling back to first qualifying year: %d\n", first_stable))
    } else {
      stop("No years meet the coverage threshold.")
    }
  } else {
    cat(sprintf("\n  First year starting %d consecutive qualifying years: %d\n",
                require_consecutive, first_stable))
  }
  list(table = cov_tbl, first_stable_year = first_stable)
}

# ============================================================================
# -- HELPER - Coverage-jump diagnostic (run manually from console)
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
            summarise(mean_n_contributing = mean(n_lakes_contributing, na.rm = TRUE),
                      mean_n_total        = mean(n_lakes_total, na.rm = TRUE),
                      .groups = "drop"))
  }
  if ("product" %in% names(df_raw)) {
    cat("\n3. 'product' field, pre vs post:\n")
    print(df_raw %>% count(period, product))
  }
  cat("\n4. total_storage_Mm3 magnitude, pre vs post:\n")
  print(df_raw %>% group_by(period) %>%
          summarise(mean_storage = mean(total_storage_Mm3, na.rm = TRUE),
                    sd_storage   = sd(total_storage_Mm3, na.rm = TRUE),
                    .groups = "drop"))
  invisible(df_raw)
}

# ============================================================================
# -- SECTION 6 --  PIPELINE 1  -  Basin-wide GeoLakes TSDI
# ============================================================================
run_pipeline_1 <- function(min_annual_coverage     = 0.70,
                           require_consecutive     = 5L,
                           calibration_start_year  = NULL,
                           quality_max             = 1L,
                           expected_obs_per_year   = NULL,
                           interp_maxgap           = P1_INTERP_MAXGAP,
                           gap_schedule_days       = P1_GAP_SCHEDULE_DAYS) {
  cat("\n\n========================================================\n")
  cat("PIPELINE 1 - Basin-wide TSDI  (GeoLakes aggregated)\n")
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
      cat(sprintf("  Quality filter (<=%d): kept %d / %d rows, dropped %d\n",
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
  } else cov$first_stable_year
  
  # Monthly calibration grid: normalise dates to calendar months, average any
  # duplicates, add absent months, and linearly interpolate ONLY short internal
  # gaps (maxgap months). Seasonal gaps (winter) stay NA by design.
  df_calib_observed <- df_input %>%
    filter(year >= start_yr) %>%
    mutate(date = lubridate::floor_date(as.Date(date), unit = "month")) %>%
    group_by(date) %>%
    summarise(storage_m3 = mean(storage_m3, na.rm = TRUE), .groups = "drop") %>%
    mutate(storage_m3 = ifelse(is.finite(storage_m3), storage_m3, NA_real_))
  
  if (nrow(df_calib_observed) < 2L || all(!is.finite(df_calib_observed$storage_m3))) {
    stop("Pipeline 1: fewer than two monthly dates with usable storage values.")
  }
  
  monthly_dates <- seq(from = min(df_calib_observed$date),
                       to   = max(df_calib_observed$date), by = "month")
  df_calib <- tibble(date = monthly_dates) %>%
    left_join(df_calib_observed, by = "date") %>%
    mutate(
      storage_observed = is.finite(storage_m3),
      storage_m3 = zoo::na.approx(storage_m3, x = as.numeric(date),
                                  maxgap = interp_maxgap, na.rm = FALSE),
      storage_interpolated = !storage_observed & is.finite(storage_m3)
    )
  
  n_interpolated <- sum(df_calib$storage_interpolated, na.rm = TRUE)
  n_unfilled     <- sum(!is.finite(df_calib$storage_m3))
  cat(sprintf("  Monthly calibration grid: %d months; interpolated %d month(s); %d month(s) remain missing.\n",
              nrow(df_calib), n_interpolated, n_unfilled))
  cat("  Note: interpolated months are estimates, not independent observations.\n")
  
  # Seasonal-sampling diagnostic
  obs_by_month <- df_calib %>% filter(storage_observed) %>%
    count(calendar_month = month(date), name = "n_observed")
  cat("\n  Observed months by calendar month:\n")
  print(as.data.frame(obs_by_month), row.names = FALSE)
  never_obs <- setdiff(1:12, obs_by_month$calendar_month)
  if (length(never_obs)) {
    cat(sprintf("  NOTE: calendar month(s) %s have NO observations: the record is seasonal\n",
                paste(never_obs, collapse = ", ")))
    cat(sprintf("  (ice/snow-limited retrievals). They cannot be interpolated (maxgap = %d) and stay NA.\n",
                interp_maxgap))
    cat("  Dry runs therefore cannot span winter under the strict gap rule; the relaxed\n")
    cat("  gap schedule entry (if any) is used as a labelled sensitivity fallback.\n")
  }
  
  cat(sprintf("\n  Full record : %s to %s (%d input rows)\n",
              min(df_input$date), max(df_input$date), nrow(df_input)))
  cat(sprintf("  Calibration : %s to %s (%d monthly steps)\n\n",
              min(df_calib$date), max(df_calib$date), nrow(df_calib)))
  
  res <- compute_tsdi(df_in             = df_calib,
                      drought_class_c   = DROUGHT_CLASS_C,
                      use_smoothing     = USE_SMOOTHING,
                      smooth_window     = SMOOTH_WINDOW,
                      blom_a            = BLOM_A,
                      gap_schedule_days = gap_schedule_days,
                      label             = "TSDI")
  
  res$coverage_table    <- cov$table
  res$calibration_start <- start_yr
  
  out_file <- file.path(OUT_DIR, "Nechako_TSDI_Pipeline1_Output.csv")
  write_csv(res$data, out_file)
  
  dir.create(file.path(OUT_DIR, "processed_data"), showWarnings = FALSE, recursive = TRUE)
  write_csv(cov$table, file.path(OUT_DIR, "processed_data", "Pipeline1_annual_coverage.csv"))
  
  pdf(file.path(OUT_DIR, "Pipeline1_Diagnostics.pdf"), width = 10, height = 6)
  plot_pipeline(res, label = "TSDI",
                basin_label = "Nechako Basin (all lakes)",
                storage_units_label = "Mm^3")
  dev.off()
  
  print_summary(res, label = "TSDI", out_file = out_file)
  invisible(res)
}

# ============================================================================
# -- SECTION 7 --  PIPELINE 2  -  Nechako Reservoir NTSDI  (LOCAL FILE ONLY)
# ============================================================================
run_pipeline_2 <- function(gap_schedule_days = P2_GAP_SCHEDULE_DAYS,
                           start_year        = NTSDI_START_YEAR) {
  if (LOCAL_DATA_TYPE %in% c("LEVEL", "LEVEL_DAILY")) {
    HYPSO_TABLE <<- load_hypso_table()
  }
  cat("\n\n========================================================\n")
  cat("PIPELINE 2 - Nechako Reservoir NTSDI\n")
  cat("  Data source: LOCAL_DATA (", LOCAL_DATA_TYPE, ")\n")
  cat("  Analysis start year (NTSDI_START_YEAR):", start_year, "\n")
  cat("========================================================\n")
  
  cat("Reading local authority file: ", LOCAL_DATA_FILE, "\n")
  
  clamp_flag <- identical(HYPSO_OUT_OF_RANGE, "CLAMP")
  
  if (LOCAL_DATA_TYPE == "STORAGE") {
    df_local <- read_csv(LOCAL_DATA_FILE, show_col_types = FALSE)
    stopifnot(LOCAL_DATE_COL %in% names(df_local), LOCAL_VALUE_COL %in% names(df_local))
    df_input <- df_local %>%
      transmute(date       = as.Date(.data[[LOCAL_DATE_COL]], format = LOCAL_DATE_FORMAT),
                storage_m3 = as.numeric(.data[[LOCAL_VALUE_COL]]) * 1e6,
                clamped    = FALSE)
    
  } else if (LOCAL_DATA_TYPE == "LEVEL") {
    df_local <- read_csv(LOCAL_DATA_FILE, show_col_types = FALSE)
    stopifnot(LOCAL_DATE_COL %in% names(df_local), LOCAL_VALUE_COL %in% names(df_local))
    rng <- range(HYPSO_TABLE$level_m)
    df_input <- df_local %>%
      transmute(date        = as.Date(.data[[LOCAL_DATE_COL]], format = LOCAL_DATE_FORMAT),
                level_m_asl = .data[[LOCAL_VALUE_COL]]) %>%
      mutate(clamped    = level_m_asl < rng[1] | level_m_asl > rng[2],
             storage_m3 = level_to_storage(level_m_asl, clamp = clamp_flag)) %>%
      dplyr::select(date, storage_m3, clamped)
    
  } else if (LOCAL_DATA_TYPE == "LEVEL_DAILY") {
    cat("  Parsing Year + Day-of-Year format and converting to storage...\n")
    # Read raw lines, trim whitespace (trailing blanks otherwise make the
    # column guesser infer a phantom 4th column), and detect the field count.
    lines_raw   <- readLines(LOCAL_DATA_FILE)
    lines_clean <- trimws(lines_raw)
    lines_clean <- lines_clean[nzchar(lines_clean)]
    
    n_fields_per_line <- utils::count.fields(textConnection(lines_clean), sep = "")
    field_tbl <- table(n_fields_per_line)
    cat("  Detected field counts across rows of", LOCAL_DATA_FILE, "(post-trim):\n")
    print(field_tbl)
    n_fields <- as.integer(names(field_tbl)[which.max(field_tbl)])
    cat("  Using n_fields =", n_fields, "(most common field count)\n")
    
    if (any(n_fields_per_line != n_fields)) {
      warning(sum(n_fields_per_line != n_fields), " row(s) in ", LOCAL_DATA_FILE,
              " do not match the dominant field count (", n_fields,
              "). These rows may parse incorrectly -- inspect them directly.")
    }
    
    if (n_fields == 3) {
      col_names_ld   <- c("Year", "DOY", "Level_m")
      col_classes_ld <- c("integer", "integer", "numeric")
    } else if (n_fields == 4) {
      col_names_ld   <- c("Year", "DOY", "Level_m", "Extra_field")
      col_classes_ld <- c("integer", "integer", "numeric", "character")
    } else {
      stop("LOCAL_DATA_FILE has an unexpected number of whitespace-delimited fields (",
           n_fields, "). Expected 3 or 4. Inspect the file before proceeding.")
    }
    
    df_raw_ld <- read.table(text = lines_clean, header = FALSE,
                            col.names = col_names_ld, colClasses = col_classes_ld,
                            strip.white = TRUE) %>% as_tibble()
    
    if (nrow(df_raw_ld) != length(lines_clean)) {
      stop("Parsed row count (", nrow(df_raw_ld), ") does not match cleaned line count (",
           length(lines_clean), "). Parsing is misaligned -- inspect the file.")
    }
    
    if (n_fields == 4) {
      cat("  NOTE: file has a 4th field beyond Year/DOY/Level_m (first 20 unique values):\n")
      print(utils::head(sort(unique(df_raw_ld$Extra_field)), 20))
      cat("  This 4th field is being DROPPED from the analysis below.\n")
    }
    
    # ---- out-of-range diagnostic (ALL years, before the start-year filter) ----
    rng <- range(HYPSO_TABLE$level_m)
    cat(sprintf("\n  Storage-curve level range: %.3f to %.3f m\n", rng[1], rng[2]))
    oor <- df_raw_ld %>%
      filter(Year <= year(Sys.Date())) %>%
      mutate(side = case_when(Level_m < rng[1] ~ "below_curve",
                              Level_m > rng[2] ~ "above_curve",
                              TRUE             ~ "in_range")) %>%
      filter(side != "in_range") %>%
      count(Year, side, name = "n_days")
    if (nrow(oor) > 0) {
      cat("  Out-of-range days by year and side:\n")
      print(oor, n = Inf)
      cat("  Totals by side:\n")
      print(oor %>% group_by(side) %>% summarise(n_days = sum(n_days), .groups = "drop") %>%
              as.data.frame(), row.names = FALSE)
      tryCatch(readr::write_csv(oor, file.path(OUT_DIR, "NTSDI_out_of_range_levels.csv")),
               error = function(e) NULL)
    } else cat("  No out-of-range level days.\n")
    cat(sprintf("  Out-of-range handling: %s\n", HYPSO_OUT_OF_RANGE))
    
    n_pre <- sum(df_raw_ld$Year < start_year)
    cat(sprintf("  Dropping %d daily rows before NTSDI_START_YEAR = %d (reservoir filling period).\n",
                n_pre, start_year))
    
    df_input <- df_raw_ld %>%
      filter(Year <= year(Sys.Date()), Year >= start_year) %>%
      mutate(date       = as.Date(paste(Year, DOY), format = "%Y %j"),
             clamped    = Level_m < rng[1] | Level_m > rng[2],
             storage_m3 = level_to_storage(Level_m, hypso = HYPSO_TABLE, clamp = clamp_flag)) %>%
      dplyr::select(date, storage_m3, clamped)
    
  } else {
    stop("LOCAL_DATA_TYPE must be 'STORAGE', 'LEVEL', or 'LEVEL_DAILY'.")
  }
  
  # For non-daily types the start-year filter is applied here
  df_input <- df_input %>% filter(!is.na(date), year(date) >= start_year)
  
  n_na <- sum(is.na(df_input$storage_m3))
  if (n_na > 0)
    cat("  WARNING:", n_na, "records produced NA storage.\n")
  n_before <- nrow(df_input)
  df_input <- df_input %>% filter(!is.na(storage_m3))
  if (n_before != nrow(df_input))
    cat("  Removed", n_before - nrow(df_input), "rows with NA storage.\n")
  cat(sprintf("  Valid storage records: %d  (%s to %s)\n",
              nrow(df_input), min(df_input$date), max(df_input$date)))
  if (sum(df_input$clamped, na.rm = TRUE) > 0)
    cat(sprintf("  Days clamped to curve limits and flagged: %d\n", sum(df_input$clamped, na.rm = TRUE)))
  
  # Monthly means: require >= 15 distinct valid days per month
  MIN_VALID_DAYS_PER_MONTH <- 15L
  df_monthly <- df_input %>%
    mutate(year = year(date), month = month(date)) %>%
    group_by(year, month) %>%
    summarise(
      n_valid_days   = n_distinct(date[!is.na(storage_m3)]),
      n_clamped_days = sum(clamped, na.rm = TRUE),
      storage_m3     = if (n_valid_days >= MIN_VALID_DAYS_PER_MONTH)
        mean(storage_m3, na.rm = TRUE) else NA_real_,
      date           = as.Date(paste(year[1], month[1], "01"), "%Y %m %d"),
      .groups = "drop"
    ) %>% arrange(date)
  
  if (nrow(df_monthly) == 0L)
    stop("No monthly observations remain after applying the valid-day threshold.")
  month_grid <- tibble(date = seq(min(df_monthly$date), max(df_monthly$date), by = "month"))
  df_in <- month_grid %>%
    left_join(df_monthly %>% dplyr::select(date, storage_m3, n_valid_days, n_clamped_days),
              by = "date") %>%
    arrange(date)
  
  # ====== Run TSDI algorithm ==============================
  res <- compute_tsdi(df_in             = df_in,
                      drought_class_c   = DROUGHT_CLASS_C,
                      use_smoothing     = USE_SMOOTHING,
                      smooth_window     = SMOOTH_WINDOW,
                      blom_a            = BLOM_A,
                      gap_schedule_days = gap_schedule_days,
                      label             = "NTSDI")
  
  if (isTRUE(res$recursion_ok)) {
    cat(sprintf("Correlation between TSD and raw NTSDI: %.3f\n", res$tsd_cor))
  }
  
  res$data <- res$data %>%
    mutate(data_source     = "LOCAL_DATA",
           local_data_type = LOCAL_DATA_TYPE,
           analysis_start_year = start_year)
  
  out_file <- file.path(OUT_DIR, "Nechako_NTSDI_Pipeline2_Output.csv")
  write_csv(res$data, out_file)
  
  pdf(file.path(OUT_DIR, "Pipeline2_Diagnostics.pdf"), width = 10, height = 6)
  plot_pipeline(
    res, label = "NTSDI",
    basin_label = paste0("Nechako Reservoir [Local file: ", LOCAL_DATA_TYPE,
                         ", from ", start_year, "]"),
    storage_units_label = "m^3"
  )
  dev.off()
  
  print_summary(res, label = "NTSDI", out_file = out_file)
  invisible(res)
}

# ============================================================================
# -- SECTION 8 --  COMPARISON PLOT  (only when PIPELINE_MODE == "BOTH")
# ============================================================================
# Uses TSD (%) and the REPORTED standardised index of each pipeline, so it works
# whichever method (p/q recursion or TSD-only Blom) each pipeline ended up with.
plot_comparison <- function(res1, res2) {
  cat("\n\n== Generating comparison plot ==\n")
  comp <- inner_join(
    res1$data %>%
      mutate(date_month = lubridate::floor_date(as.Date(date), "month")) %>%
      dplyr::select(date_month, TSD_p1 = TSD, idx_p1 = index_reported),
    res2$data %>%
      mutate(date_month = lubridate::floor_date(as.Date(date), "month")) %>%
      dplyr::select(date_month, TSD_p2 = TSD, idx_p2 = index_reported),
    by = "date_month"
  ) %>% dplyr::rename(date = date_month) %>%
    filter(is.finite(idx_p1), is.finite(idx_p2))
  
  if (nrow(comp) == 0) {
    warning("No overlapping dates with finite values in both pipelines. Skipping comparison plot.")
    return(invisible(NULL))
  }
  cat("Overlapping records:", nrow(comp), "\n")
  
  lab1 <- paste0("Pipeline 1: basin lakes [", res1$method, "]")
  lab2 <- paste0("Pipeline 2: Nechako Reservoir [", res2$method, "]")
  cols <- setNames(c("#2c7bb6", "#d7191c"), c(lab1, lab2))
  note <- "P1 = basin GeoLakes aggregate (reservoir inclusion unverified); P2 = regulated reservoir record"
  
  comp_tsd <- bind_rows(
    comp %>% transmute(date, value = TSD_p1, Index = lab1),
    comp %>% transmute(date, value = TSD_p2, Index = lab2))
  pc1 <- ggplot(comp_tsd, aes(date, value, colour = Index)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_line(linewidth = 0.7, alpha = 0.85) +
    scale_colour_manual(values = cols) +
    labs(title = "Comparison - Total Storage Deficit (TSD %)",
         subtitle = paste0("Overlap: ", format(min(comp$date)), " - ", format(max(comp$date)),
                           "\n", note),
         x = NULL, y = "TSD (%)", colour = NULL) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom")
  print(pc1)
  
  comp_idx <- bind_rows(
    comp %>% transmute(date, value = idx_p1, Index = lab1),
    comp %>% transmute(date, value = idx_p2, Index = lab2))
  pc2 <- ggplot(comp_idx, aes(date, value, colour = Index)) +
    geom_hline(yintercept =  0,   linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = -1.0, linetype = "dotted", colour = "#fc8d59") +
    geom_hline(yintercept = -1.5, linetype = "dotted", colour = "#d73027") +
    geom_line(linewidth = 0.7, alpha = 0.85) +
    scale_colour_manual(values = cols) +
    labs(title = "Comparison - reported standardised index",
         subtitle = "Dotted lines: moderate (-1.0) and severe (-1.5) thresholds",
         x = NULL, y = "Standardised index [N(0,1)]", colour = NULL) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom")
  print(pc2)
  
  r_p <- round(cor(comp$idx_p1, comp$idx_p2, use = "complete.obs"), 3)
  r_s <- round(cor(comp$idx_p1, comp$idx_p2, use = "complete.obs", method = "spearman"), 3)
  pc3 <- ggplot(comp, aes(idx_p1, idx_p2)) +
    geom_point(alpha = 0.4, size = 1.5, colour = "#404040") +
    geom_smooth(method = "lm", se = TRUE, colour = "#2c7bb6", linewidth = 0.9, formula = y ~ x) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4, size = 4,
             label = paste0("Pearson r = ", r_p, "\nSpearman rho = ", r_s)) +
    labs(title = "Scatter - reported index, Pipeline 1 vs Pipeline 2",
         subtitle = "Dashed = 1:1 line; blue = OLS fit",
         x = paste0("Pipeline 1 (", res1$method, ")"),
         y = paste0("Pipeline 2 (", res2$method, ")")) +
    theme_bw(base_size = 12)
  print(pc3)
  
  cat("Pearson r (reported indices):", r_p, " | Spearman rho:", r_s, "\n")
}

# ============================================================================
# -- SECTION 9 --  MAIN EXECUTION
# ============================================================================
run_main <- function() {
  sink(file.path(OUT_DIR, "console_log.txt"), append = TRUE, split = TRUE)
  on.exit({
    while (grDevices::dev.cur() > 1L) grDevices::dev.off()
    while (sink.number() > 0L) sink()
  }, add = TRUE)
  res_p1 <- NULL
  res_p2 <- NULL
  
  if (PIPELINE_MODE %in% c("BASIN", "BOTH")) {
    # A hard failure here (missing file, <3 finite TSD values, ...) must not stop
    # the independent reservoir pipeline. A p/q calibration failure is NOT a
    # hard failure any more: compute_tsdi() falls back to the TSD-only index.
    res_p1 <- tryCatch(
      run_pipeline_1(quality_max = 1L),
      error = function(e) {
        msg <- conditionMessage(e)
        warning("Pipeline 1 (basin-wide TSDI) was skipped: ", msg, call. = FALSE)
        cat("\nPIPELINE 1 STATUS: FAILED / SKIPPED\nReason: ", msg, "\n", sep = "")
        NULL
      }
    )
    
    res_reg <- tryCatch(check_regulated_lakes(), error = function(e) {
      warning("Lake-regulation diagnostic failed and was skipped: ", conditionMessage(e))
      NULL
    })
    
    if (!is.null(res_p1)) {
      plot_index_timeseries(
        dates_vec    = res_p1$data$date,
        values_vec   = res_p1$data$index_reported,
        index_label  = "TSDI",
        title_label  = "TSDI (Basin-Wide)",
        out_dir      = OUT_DIR,
        method_label = res_p1$method
      )
    } else {
      cat("Pipeline 1 outputs and comparison requiring TSDI will not be created.\n")
    }
  }
  
  if (PIPELINE_MODE %in% c("NECHAKO_RESERVOIR", "BOTH")) {
    res_p2 <- tryCatch(
      run_pipeline_2(),
      error = function(e) {
        msg <- conditionMessage(e)
        warning("Pipeline 2 (reservoir NTSDI) failed: ", msg, call. = FALSE)
        cat("\nPIPELINE 2 STATUS: FAILED\nReason: ", msg, "\n", sep = "")
        NULL
      }
    )
    if (!is.null(res_p2)) {
      plot_index_timeseries(
        dates_vec    = res_p2$data$date,
        values_vec   = res_p2$data$index_reported,
        index_label  = "NTSDI",
        title_label  = "NTSDI (Nechako Reservoir)",
        out_dir      = OUT_DIR,
        method_label = res_p2$method
      )
    }
  }
  
  if (PIPELINE_MODE == "BOTH" && !is.null(res_p1) && !is.null(res_p2)) {
    pdf(file.path(OUT_DIR, "Comparison_Plots.pdf"), width = 10, height = 6)
    plot_comparison(res_p1, res_p2)
    dev.off()
  }
  
  # ---- pipeline status + which method was reported ---------------------------
  status_of <- function(res, requested)
    if (!is.null(res)) "SUCCESS" else if (requested) "FAILED / SKIPPED" else "NOT REQUESTED"
  pipeline_status <- data.frame(
    Pipeline = c("Pipeline 1: basin-wide TSDI", "Pipeline 2: reservoir NTSDI"),
    Status   = c(status_of(res_p1, PIPELINE_MODE %in% c("BASIN", "BOTH")),
                 status_of(res_p2, PIPELINE_MODE %in% c("NECHAKO_RESERVOIR", "BOTH"))),
    ReportedMethod = c(if (!is.null(res_p1)) res_p1$method else NA_character_,
                       if (!is.null(res_p2)) res_p2$method else NA_character_),
    PQ_Recursion_Valid = c(if (!is.null(res_p1)) res_p1$recursion_ok else NA,
                           if (!is.null(res_p2)) res_p2$recursion_ok else NA),
    p = c(if (!is.null(res_p1)) res_p1$p else NA_real_,
          if (!is.null(res_p2)) res_p2$p else NA_real_),
    q = c(if (!is.null(res_p1)) res_p1$q else NA_real_,
          if (!is.null(res_p2)) res_p2$q else NA_real_),
    GapDaysUsed = c(if (!is.null(res_p1)) res_p1$gap_days_used else NA_real_,
                    if (!is.null(res_p2)) res_p2$gap_days_used else NA_real_),
    Cor_TSD_RawIndex = c(if (!is.null(res_p1)) res_p1$tsd_cor else NA_real_,
                         if (!is.null(res_p2)) res_p2$tsd_cor else NA_real_),
    Notes = c(if (!is.null(res_p1)) paste(res_p1$failure_notes, collapse = " | ") else "",
              if (!is.null(res_p2)) paste(res_p2$failure_notes, collapse = " | ") else ""),
    stringsAsFactors = FALSE
  )
  utils::write.csv(pipeline_status, file.path(OUT_DIR, "pipeline_execution_status.csv"),
                   row.names = FALSE, na = "")
  cat("\n")
  print(pipeline_status[, c("Pipeline", "Status", "ReportedMethod", "PQ_Recursion_Valid", "p", "q")],
        row.names = FALSE)
  
  spatial_support_manifest <- data.frame(
    Pipeline = c("Pipeline 1: TSDI", "Pipeline 2: NTSDI"),
    InputFile = c(P1_INPUT_FILE, LOCAL_DATA_FILE),
    ReportedMethod = c(if (!is.null(res_p1)) res_p1$method else NA_character_,
                       if (!is.null(res_p2)) res_p2$method else NA_character_),
    SpatialSupport = c(
      "Basin-wide aggregate supplied in the input file (date, total_storage_Mm3); open-water-season observations only",
      paste0("Single Nechako Reservoir local record (", LOCAL_DATA_TYPE,
             ", from ", NTSDI_START_YEAR, "); not a basin-wide average")
    ),
    SpatialAggregationPerformedHere = c("NO; upstream input aggregation is inherited",
                                        "NO; one reservoir time series"),
    Interpretation = c(
      "Basin-wide lake-storage deficit index on the pre-aggregated input series; winter months are unobserved (NA)",
      "Reservoir-specific storage/level deficit index; interpret as regulated reservoir conditions"
    ),
    RasterBased = c("NO", "NO"),
    BoundaryCaveat = c(
      "No spatial reassignment of source lakes to basin/subbasin polygons; reservoir inclusion in the aggregate is unverified",
      "No station/lake averaging or subbasin assignment is performed"
    ),
    stringsAsFactors = FALSE
  )
  write.csv(spatial_support_manifest, file.path(OUT_DIR, "spatial_support_metadata.csv"),
            row.names = FALSE, na = "")
  cat("Spatial-support metadata saved to:",
      file.path(normalizePath(OUT_DIR), "spatial_support_metadata.csv"), "\n")
  
  cat("\n\nDone. All outputs saved to:", normalizePath(OUT_DIR), "\n")
  cat("Index time series plots saved to:", file.path(normalizePath(OUT_DIR), "figures"), "\n")
}

tryCatch(
  run_main(),
  error = function(e) {
    message("NTSDI script failed: ", conditionMessage(e))
    stop(e)
  },
  finally = {
    while (grDevices::dev.cur() > 1L) grDevices::dev.off()
    while (sink.number() > 0L) sink()
  }
)