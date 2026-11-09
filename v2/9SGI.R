#=============================================================================
# 9SGI.R  -  SGI (Standardised Groundwater level Index), Nechako River Basin
#
# Method : Bloomfield & Marchant (2013), HESS 17, 4769-4787,
#          doi:10.5194/hess-17-4769-2013. Non-parametric normal-scores
#          transform of monthly groundwater level, applied per well and per
#          calendar month; basin value = median SGI across qualifying wells.
#
# Data   : 1. BC PGOWN monthly medians (primary, Environmental Reporting BC)
#          2. NRCan GIN/NGWD           [stub]
#          3. ECCC WSC (tidyhydat)     [stub]
#          4. BC Real-Time Water Data  [stub]
#          5. BC GWELLS API            [stub]
#
# FRAMEWORK CONTRACT (same rules as Scripts 1, 3, 4, 5, 6, 7, 8):
#   - Working directory = project root; every path is derived from getwd();
#     shared helpers come from DROUGHT_ANALYSIS_utils.R.
#   - Basin = dissolved union of the 7 authoritative subbasins
#     (nechako_subbasin_union(load_nechako_subbasins())), identical to every
#     other script. The old Spatial/nechakoBound_dissolve.kmz boundary is no
#     longer used, so wells are screened against the same polygon as the
#     ERA5-Land basin mask.
#   - Analysis window 1950-01 to 2025-12 (as Scripts 4 and 16); calibration
#     reference = the WHOLE record of each well inside that window (as SWEI,
#     SSPI, SSI-SM, SSMI). No sub-period climatology.
#   - Plotting position: Gringorten (r - 0.44) / (n + 0.12), the empirical
#     rule used by SPI/SRI/SWEI/SSPI/SSI-SM (replaces the Hazen (2r-1)/(2n)
#     position of the previous version). Values clamped to +/-4.75, as
#     SPI/SRI/SSI-SM.
#   - Minimum sample: SGI_MIN_YEARS_PER_MONTH years of the same calendar month
#     per well (10, same as MIN_YEARS_PER_MONTH in 6StreamFlow_Index_WSCLocal.R,
#     the other point-support index). Calendar months below this are NA.
#   - Drought events: identify_drought_events_ts() on the full monthly
#     calendar (onset < DROUGHT_ONSET, end >= DROUGHT_END, scale-1 minimum
#     duration, intensity = mean(DROUGHT_END - x), severity = I x D). Missing
#     months terminate an event; gaps are never bridged.
#   - Classes: McKee et al. (1993) W3..D3 classes, the same bins written to
#     the SPI / SSI-SM summary files (see sgi_class()).
#   - Spatial support: SGI is a WELL-NETWORK MEDIAN, not an area-weighted
#     basin mean. 14basin_timeseries.R (Part 6) therefore treats it like the
#     station SSI and TSDI/NTSDI series, not like the gridded basin means.
#
# OUTPUTS (sgi_results/):
#   sgi_nechako_multisource_basin.csv  monthly basin series, long format
#                                      (date, month_label, sgi_median, ...);
#                                      read by 10MSDI.r and 14basin_timeseries.R
#   sgi_01_basin_median_by_month.csv   Year x Jan..Dec (Script 4/8 wide layout)
#   sgi_well_timeseries.csv            per-well GWL and SGI, long format
#   sgi_selected_wells.csv             qualifying wells, coordinates, record stats
#   sgi_drought_events.csv             framework drought events (basin series)
#   sgi_summary.txt                    run summary / QA
#   figures/sgi_basin_timeseries.png   basin series, framework event shading
#   figures/sgi_basin_extent_check.png basin polygon + candidate wells
#   figures/wells/sgi_well_<id>.png    one panel per qualifying well
#   cache/                             downloaded PGOWN tables
#
# EXECUTION ORDER: 1, 2b, 3, 4, 5, 6, 7, 8, 9 (this script) -> 10 -> 11 ... 16
#=============================================================================

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
utils_load_packages(c("httr", "dplyr", "tidyr", "readr", "purrr", "lubridate",
                      "ggplot2", "sf", "terra"))

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a

# SGI is a monthly state index with no accumulation window, so its events use
# the scale-1 minimum duration of EVENT_MIN_DURATION_BY_SCALE (as SSI, TSDI,
# NTSDI and SSPI-1).
SGI_EVENT_SCALE <- 1L

# ---- Paths -------------------------------------------------------------------
OUT_DIR      <- file.path(WD_PATH, "sgi_results")
FIG_DIR      <- file.path(OUT_DIR, "figures")
WELL_FIG_DIR <- file.path(FIG_DIR, "wells")
CACHE_DIR    <- file.path(OUT_DIR, "cache")
for (d in c(OUT_DIR, FIG_DIR, WELL_FIG_DIR, CACHE_DIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Date-class infinities for annotate("rect") on date axes (as Script 7).
DATE_NEG_INF <- structure(-Inf, class = "Date")
DATE_POS_INF <- structure( Inf, class = "Date")

#---- 1. Configuration -------------------------------------------------------
cfg <- list(
  # === Source toggles (PGOWN only enabled by default) ===
  src_pgown       = TRUE,   # BC PGOWN monthly medians (working)
  src_gin_ngwd    = FALSE,  # NRCan GIN/NGWD - requires manual cache
  src_wsc         = FALSE,  # Water Survey of Canada (tidyhydat)
  src_bc_realtime = FALSE,  # BC Real-Time (ArcGIS)
  src_gwells      = FALSE,  # BC GWELLS API
  
  # === PGOWN (primary) ===
  url_pgown_ts   = "https://catalogue.data.gov.bc.ca/dataset/a74f1b97-17f7-499b-84e7-6455e169e425/resource/35162de7-803e-42f2-93d4-563330bb7dc9/download/monthlywell_ts.csv",
  url_pgown_attr = "https://catalogue.data.gov.bc.ca/dataset/a74f1b97-17f7-499b-84e7-6455e169e425/resource/a8933793-eadb-4a9c-992c-da4f6ac8ca51/download/gw_well_table.csv",
  pgown_dataset_id = "a74f1b97-17f7-499b-84e7-6455e169e425",
  
  # === GIN/NGWD (Federal) ===
  gin_cache      = file.path(CACHE_DIR, "gin_nechako_subset.csv"),
  
  # === WSC (tidyhydat) ===
  wsc_radius_km  = 50,
  wsc_daily_agg  = "median",
  
  # === BC Real-Time (ArcGIS) ===
  bc_rt_layer    = "https://maps.gov.bc.ca/arcserver/rest/services/pub/Water/BC_Water_Levels/MapServer/1",
  
  # === GWELLS API ===
  gwells_api_base = "https://apps.nrs.gov.bc.ca/gwells/api/v2",
  gwells_min_obs  = 5,
  
  # === Framework analysis window / reference (whole record) ===
  analysis_start = as.Date("1950-01-01"),
  analysis_end   = as.Date("2025-12-01"),
  
  # === Standardisation ===
  min_years_per_month = 10L,    # years of the same calendar month per well
  z_clip              = 4.75,   # same clamp as SPI / SRI / SSI-SM
  
  # === Well quality filters ===
  min_years    = 10,            # record span (years) for PGOWN wells
  max_gap_frac = 0.50,          # max fraction of missing months in the span
  min_wells_per_month = 1L,     # basin median reported only with >= this many wells
  
  # === Basin extent sanity check (coarse; the 7-subbasin union is authoritative)
  # Nechako watershed is ~47,000 km^2; the Stuart-Takla headwaters push the
  # northern extent to ~56.2 N.
  basin_area_km2_expected_range = c(20000, 70000),
  basin_bbox_lat_range = c(52.0, 57.0),
  basin_bbox_lon_range = c(-128.5, -121.0),
  
  # === Sign convention ===
  # TRUE  -> PGOWN med_gwl is depth-to-water below ground (larger = drier), so
  #          it is negated: positive SGI = wetter (Bloomfield & Marchant 2013).
  # FALSE -> med_gwl is already an elevation/head (masl); no sign flip.
  # Automated checks (magnitude heuristic, CKAN metadata wording, known 2023
  # drought) warn loudly on disagreement, but confirm against the PGOWN
  # field definitions before trusting the classes.
  invert_gw_sign = TRUE,
  
  # === Output ===
  target_months = NULL,          # optional character vector of "YYYY-MM" to print
  basin_file    = "sgi_nechako_multisource_basin.csv"
)

stopifnot(cfg$min_years_per_month >= 5L, cfg$z_clip > 0,
          cfg$analysis_start < cfg$analysis_end)

#---- 2. Helpers -------------------------------------------------------------
fetch_csv <- function(url, cache_path, desc = "") {
  if (!is.null(cache_path) && file.exists(cache_path)) {
    cat(sprintf("  [cache] %s\n", cache_path))
    return(readr::read_csv(cache_path, show_col_types = FALSE))
  }
  cat(sprintf("  [download] %s ...\n", desc))
  dest <- if (!is.null(cache_path)) cache_path else tempfile(fileext = ".csv")
  resp <- httr::GET(url, httr::write_disk(dest, overwrite = TRUE), httr::timeout(120))
  httr::stop_for_status(resp, task = paste("download", desc))
  readr::read_csv(dest, show_col_types = FALSE)
}

detect_col <- function(df, patterns, label) {
  m <- grep(paste(patterns, collapse = "|"), names(df), ignore.case = TRUE, value = TRUE)
  if (!length(m)) stop(sprintf("Cannot find '%s' column. Tried: %s", label, paste(patterns, collapse = ", ")))
  if (length(m) > 1) cat(sprintf("  '%s': multiple matches (%s); using '%s'\n", label, paste(m, collapse = ", "), m[1]))
  m[1]
}

month_start <- function(d) as.Date(format(as.Date(d), "%Y-%m-01"))

# Gringorten normal scores for one well x calendar month sample (whole record).
# Zero variance -> 0 (neutral), as the SPI/SRI zero-variance rule.
sgi_normal_scores <- function(x, min_n = cfg$min_years_per_month, z_clip = cfg$z_clip) {
  z  <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  n  <- sum(ok)
  if (n < min_n) return(z)
  v <- stats::var(x[ok])
  if (!is.finite(v) || v < .Machine$double.eps) { z[ok] <- 0; return(z) }
  r <- rank(x[ok], ties.method = "average")
  p <- (r - 0.44) / (n + 0.12)
  z[ok] <- pmax(-z_clip, pmin(z_clip, stats::qnorm(p)))
  z
}

# McKee et al. (1993) classes, identical bins to the SPI / SSI-SM summaries.
sgi_class <- function(x) {
  dplyr::case_when(
    !is.finite(x) ~ NA_character_,
    x >=  2.0     ~ "W3 Extremely wet",
    x >=  1.5     ~ "W2 Severely wet",
    x >=  1.0     ~ "W1 Moderately wet",
    x >  -1.0     ~ "N0 Near normal",
    x >  -1.5     ~ "D1 Moderately dry",
    x >  -2.0     ~ "D2 Severely dry",
    TRUE          ~ "D3 Extremely dry"
  )
}

# Full monthly calendar so identify_drought_events_ts() sees missing months.
complete_monthly <- function(df, from = min(df$date), to = max(df$date)) {
  cal <- data.frame(date = seq(month_start(from), month_start(to), by = "month"))
  dplyr::left_join(cal, df, by = "date")
}

# Framework event table, one row per event (columns as Script 7 SWSI events).
sgi_event_table <- function(dates, values, n_wells) {
  ev <- identify_drought_events_ts(dates, values, scale = SGI_EVENT_SCALE)
  in_ev <- attr(ev, "in_event")
  if (is.null(in_ev) || length(in_ev) != length(values)) in_ev <- rep(FALSE, length(values))
  in_ev[is.na(in_ev)] <- FALSE
  r <- rle(in_ev); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1L
  runs <- which(r$values)
  if (length(runs) != nrow(ev))
    warning(sprintf("SGI: %d event runs vs %d events from identify_drought_events_ts().",
                    length(runs), nrow(ev)), call. = FALSE)
  tab <- dplyr::bind_rows(lapply(seq_along(runs), function(k) {
    idx <- starts[runs[k]]:ends[runs[k]]; v <- values[idx]
    inten <- mean(DROUGHT_END - v, na.rm = TRUE)
    data.frame(index = "SGI", event_id = k,
               start_date = dates[min(idx)], end_date = dates[max(idx)],
               duration_months = length(idx),
               sgi_min = min(v, na.rm = TRUE), sgi_mean = mean(v, na.rm = TRUE),
               intensity = inten, severity = inten * length(idx),
               peak_class = sgi_class(min(v, na.rm = TRUE)),
               mean_n_wells = mean(n_wells[idx], na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
  list(events = tab, in_event = in_ev, min_duration = attr(ev, "min_duration"))
}

# Standalone time-series PNG - same layout, event rule and shading as
# plot_index_timeseries() in 7WBInex_NTSDI_GeoLakes_LocalRS.r.
plot_index_timeseries <- function(dates_vec, values_vec, index_label, title_label,
                                  file, subtitle_extra = NULL,
                                  clamp_floor = NULL, clamp_label = NULL,
                                  event_scale = SGI_EVENT_SCALE,
                                  width = 11, height = 5.5) {
  df <- data.frame(date = as.Date(dates_vec), value = values_vec)
  df <- df[order(df$date), ]
  if (!any(is.finite(df$value))) {
    cat(sprintf("  [PLOT] Skipping %s time series: no finite values\n", index_label))
    return(invisible(NULL))
  }
  # Detect on the full monthly calendar BEFORE dropping missing values.
  ev       <- identify_drought_events_ts(df$date, df$value, scale = event_scale)
  n_events <- nrow(ev)
  min_dur  <- attr(ev, "min_duration")
  df$below <- attr(ev, "in_event")
  df <- df[is.finite(df$value), ]
  month_id <- as.integer(format(df$date, "%Y")) * 12L + as.integer(format(df$date, "%m"))
  df$grp   <- cumsum(c(1L, diff(as.integer(df$below)) != 0L | diff(month_id) != 1L))
  
  y_lo <- min(-2.2, floor(min(df$value) * 10) / 10,
              if (!is.null(clamp_floor)) floor(clamp_floor * 10) / 10 - 0.1 else Inf)
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
                        linetype = "dotted", linewidth = 0.4)
  
  # Plotting-position floor of the shortest well/month sample used (see Section 7).
  if (!is.null(clamp_floor))
    p <- p +
    ggplot2::geom_hline(yintercept = clamp_floor, color = "grey35",
                        linetype = "dashed", linewidth = 0.45) +
    ggplot2::annotate("text", x = min(df$date), y = clamp_floor, hjust = 0, vjust = -0.6,
                      label = clamp_label %||% sprintf("plotting-position floor (SGI = %.3f)", clamp_floor),
                      color = "grey35", size = 3.2, fontface = "italic")
  
  if (any(df$below))
    p <- p + ggplot2::geom_ribbon(data = subset(df, below),
                                  ggplot2::aes(ymin = value, ymax = 0, group = grp),
                                  fill = "#C0504D", alpha = 0.45)
  
  p <- p +
    ggplot2::geom_line(color = "#1F77B4", linewidth = 0.45) +
    ggplot2::scale_x_date(date_breaks = "5 years", date_labels = "%Y",
                          expand = ggplot2::expansion(mult = c(0.01, 0.015))) +
    ggplot2::scale_y_continuous(breaks = seq(-5, 5, 1)) +
    ggplot2::coord_cartesian(ylim = c(y_lo, y_hi)) +
    ggplot2::labs(
      title    = sprintf("%s Time Series  |  NECHAKO BASIN", title_label),
      subtitle = sprintf("%d drought events  (onset < %.2f, end >= %.2f, min %d months)%s",
                         n_events, DROUGHT_ONSET, DROUGHT_END, min_dur,
                         if (!is.null(subtitle_extra)) paste0("  |  ", subtitle_extra) else ""),
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
  ggplot2::ggsave(file, p, width = width, height = height, dpi = 200, bg = "white")
  cat(sprintf("  [PLOT] %s saved: %s  (%d events)\n", index_label, basename(file), n_events))
  invisible(p)
}

#---- 3. Basin boundary (framework: union of the 7 subbasins) ----------------
cat("\n============================================================\n")
cat("SGI CALCULATION (Bloomfield & Marchant 2013)\n")
cat("============================================================\n\n")

# load_nechako_subbasins() returns the framework's own "nechako_subbasins"
# object (not a SpatVector), so only the UNION is passed to terra directly, as
# in every other script. The individual polygons are used for well labels and
# the extent figure only, extracted tolerantly; if that fails, SGI is computed
# unchanged and only the subbasin label is skipped.
subbasins_to_vect <- function(x) {
  to_v <- function(e) {
    if (inherits(e, "SpatVector")) return(e)
    if (inherits(e, c("sf", "sfc"))) return(terra::vect(e))
    NULL
  }
  v <- to_v(x)
  if (!is.null(v)) return(v)
  if (is.list(x)) {
    parts <- lapply(x, function(e) tryCatch(to_v(e), error = function(err) NULL))
    keep  <- !vapply(parts, is.null, logical(1))
    if (!any(keep)) return(NULL)
    parts <- parts[keep]
    nm    <- names(parts)
    parts <- lapply(seq_along(parts), function(i) {
      pv <- parts[[i]]
      if (length(parts) > 1L && !is.null(nm) && nzchar(nm[i]) && !"region_name" %in% names(pv))
        pv$region_name <- nm[i]
      pv
    })
    return(if (length(parts) == 1L) parts[[1]] else do.call(rbind, parts))
  }
  NULL
}

subbasins_obj <- load_nechako_subbasins()
basin    <- nechako_subbasin_union(subbasins_obj)
basin_ll <- terra::project(basin, "EPSG:4326")
subb_ll  <- tryCatch({
  v <- subbasins_to_vect(subbasins_obj)
  if (is.null(v)) NULL else terra::project(v, "EPSG:4326")
}, error = function(e) NULL)
cat("Basin boundary loaded (union of 7 subbasins)\n")
if (is.null(subb_ll))
  cat("  NOTE: individual subbasin polygons could not be extracted from load_nechako_subbasins();\n",
      "       wells are screened against the basin union and the subbasin label is left NA.\n")

#---- 3b. Basin extent sanity check ------------------------------------------
cat("\n=== Basin extent sanity check ===\n")
basin_area_km2 <- sum(terra::expanse(basin_ll, unit = "km"))
bb <- as.vector(terra::ext(basin_ll))   # xmin, xmax, ymin, ymax
cat(sprintf("  Basin polygon area: %.0f km^2 (expected range: %.0f - %.0f km^2)\n",
            basin_area_km2, cfg$basin_area_km2_expected_range[1], cfg$basin_area_km2_expected_range[2]))
cat(sprintf("  Basin bbox: lat %.3f to %.3f | lon %.3f to %.3f\n", bb[3], bb[4], bb[1], bb[2]))
extent_ok <- basin_area_km2 >= cfg$basin_area_km2_expected_range[1] &&
  basin_area_km2 <= cfg$basin_area_km2_expected_range[2] &&
  bb[3] >= cfg$basin_bbox_lat_range[1] && bb[4] <= cfg$basin_bbox_lat_range[2] &&
  bb[1] >= cfg$basin_bbox_lon_range[1] && bb[2] <= cfg$basin_bbox_lon_range[2]
if (!extent_ok) {
  cat(strrep("!", 70), "\n")
  cat("  WARNING: basin polygon area/extent falls OUTSIDE the expected Nechako range.\n")
  cat("  Check the subbasin layer before trusting the well count.\n")
  cat(strrep("!", 70), "\n")
} else {
  cat("  Basin area and bbox fall within the expected Nechako range.\n")
}

# Subbasin attribute used to label wells (first character column, tolerant).
subb_name_col <- if (is.null(subb_ll)) NA_character_ else {
  df_sb <- as.data.frame(subb_ll)
  cand  <- intersect(c("region_name", "region_id", "name", "Name", "NAME"), names(df_sb))
  if (length(cand)) cand[1] else names(df_sb)[vapply(df_sb, is.character, logical(1))][1]
}

#---- 4. Source loaders -------------------------------------------------------
# Every loader returns NULL or a data frame with columns:
#   well_num, well_name, lat, lon, date (month start), gw (sign-corrected), source_id
# so the merge, spatial filter and well table below work for any source.

check_sign_convention <- function(raw_vals, col_val) {
  cat(sprintf("\n  --- Sign-convention check (invert_gw_sign = %s) ---\n", cfg$invert_gw_sign))
  raw_vals <- suppressWarnings(as.numeric(raw_vals))
  cat(sprintf("  Sample raw '%s' values: %s\n", col_val,
              paste(round(utils::head(raw_vals, 10), 3), collapse = ", ")))
  
  # Heuristic 1: magnitude. Depth-to-water in BC PGOWN is typically 0-30 m;
  # elevation/head in the BC interior is several hundred metres.
  med_abs <- stats::median(abs(raw_vals), na.rm = TRUE)
  says_depth <- med_abs < 50
  cat(sprintf("  Median |%s|: %.2f -> heuristic suggests: %s\n", col_val, med_abs,
              if (says_depth) "depth-below-ground (should invert)" else "elevation/head (should NOT invert)"))
  if (says_depth != cfg$invert_gw_sign) {
    cat(strrep("!", 70), "\n")
    cat("  WARNING: magnitude heuristic DISAGREES with cfg$invert_gw_sign.\n")
    cat("  Every drought/wet class downstream may be inverted.\n")
    cat(strrep("!", 70), "\n")
  }
  
  # Heuristic 2: wording of the PGOWN CKAN metadata (best effort).
  dict_hit <- tryCatch({
    meta_url <- sprintf("https://catalogue.data.gov.bc.ca/api/3/action/package_show?id=%s",
                        cfg$pgown_dataset_id)
    resp <- httr::GET(meta_url, httr::timeout(20))
    if (httr::status_code(resp) == 200) {
      meta <- httr::content(resp, as = "parsed", type = "application/json")
      txt <- tolower(paste(meta$result$notes %||% "",
                           paste(sapply(meta$result$resources, function(r) r$description %||% ""),
                                 collapse = " ")))
      depth_kw <- c("depth to water", "depth below ground", "below ground surface", "btoc", "depth to gw")
      elev_kw  <- c("elevation", "metres above sea level", "masl", "above sea level", "hydraulic head")
      list(checked = TRUE,
           found_depth = depth_kw[vapply(depth_kw, grepl, logical(1), x = txt, fixed = TRUE)],
           found_elev  = elev_kw[vapply(elev_kw, grepl, logical(1), x = txt, fixed = TRUE)])
    } else list(checked = FALSE, status = httr::status_code(resp))
  }, error = function(e) list(checked = FALSE, error = conditionMessage(e)))
  
  if (isTRUE(dict_hit$checked)) {
    if (length(dict_hit$found_depth) || length(dict_hit$found_elev)) {
      cat(sprintf("  CKAN metadata mentions: depth terms [%s] | elevation terms [%s]\n",
                  paste(dict_hit$found_depth, collapse = ", "), paste(dict_hit$found_elev, collapse = ", ")))
      if (length(dict_hit$found_depth) && !length(dict_hit$found_elev) && !cfg$invert_gw_sign)
        cat("  WARNING: metadata suggests depth-below-ground, but invert_gw_sign = FALSE.\n")
      if (length(dict_hit$found_elev) && !length(dict_hit$found_depth) && cfg$invert_gw_sign)
        cat("  WARNING: metadata suggests elevation/head, but invert_gw_sign = TRUE.\n")
    } else cat("  CKAN metadata fetched but contains no explicit depth/elevation wording.\n")
  } else {
    cat(sprintf("  Could not fetch/parse CKAN metadata (%s).\n",
                dict_hit$error %||% paste("HTTP", dict_hit$status)))
  }
  cat("  ACTION: confirm med_gwl's convention against the PGOWN field definitions\n")
  cat("  before trusting the SGI classes.\n")
  cat(strrep("-", 70), "\n")
  invisible(says_depth)
}

load_pgown <- function() {
  if (!cfg$src_pgown) return(NULL)
  cat("Loading PGOWN (Environmental Reporting BC) ...\n")
  
  attrs_raw <- fetch_csv(cfg$url_pgown_attr, file.path(CACHE_DIR, "pgown_attrs.csv"), "well attributes")
  names(attrs_raw) <- tolower(names(attrs_raw))
  idx_cols <- grep("^\\.\\.\\.\\d+$", names(attrs_raw), value = TRUE)
  if (length(idx_cols)) attrs_raw <- dplyr::select(attrs_raw, -dplyr::all_of(idx_cols))
  
  col_well <- detect_col(attrs_raw, c("^well_num$", "well_id", "obs.*well"), "well ID")
  col_lat  <- detect_col(attrs_raw, c("^lat$", "^latitude$"), "latitude")
  col_lon  <- detect_col(attrs_raw, c("^lon$", "^long$", "^longitude$"), "longitude")
  col_name <- tryCatch(detect_col(attrs_raw, c("^well_name$", "name"), "well name"),
                       error = function(e) NULL)
  
  attrs <- data.frame(
    well_num  = as.character(attrs_raw[[col_well]]),
    lat       = suppressWarnings(as.numeric(attrs_raw[[col_lat]])),
    lon       = suppressWarnings(as.numeric(attrs_raw[[col_lon]])),
    well_name = if (!is.null(col_name)) as.character(attrs_raw[[col_name]]) else NA_character_,
    stringsAsFactors = FALSE
  ) |>
    dplyr::filter(is.finite(lat), is.finite(lon)) |>
    dplyr::distinct(well_num, .keep_all = TRUE)
  
  ts_raw <- fetch_csv(cfg$url_pgown_ts, file.path(CACHE_DIR, "pgown_ts.csv"), "monthly medians")
  names(ts_raw) <- tolower(names(ts_raw))
  col_wn  <- detect_col(ts_raw, c("^well_num$", "well_id"), "well ID")
  col_dt  <- detect_col(ts_raw, c("^date$", "datetime"), "date")
  col_val <- detect_col(ts_raw, c("med_gwl", "gwl", "value", "level"), "GWL")
  
  check_sign_convention(ts_raw[[col_val]], col_val)
  
  gw_raw <- suppressWarnings(as.numeric(ts_raw[[col_val]]))
  ts <- data.frame(
    well_num = as.character(ts_raw[[col_wn]]),
    date     = month_start(ts_raw[[col_dt]]),
    gw       = if (cfg$invert_gw_sign) -gw_raw else gw_raw,
    source_id = "PGOWN",
    stringsAsFactors = FALSE
  ) |>
    dplyr::filter(!is.na(date), is.finite(gw))
  if (!nrow(ts)) { cat("  No valid PGOWN records after cleaning.\n"); return(NULL) }
  
  out <- dplyr::inner_join(ts, attrs, by = "well_num")
  if (!nrow(out)) { cat("  No PGOWN wells have coordinates in the attribute table.\n"); return(NULL) }
  cat(sprintf("  PGOWN: %d wells with coordinates, %d records (before spatial filter)\n",
              dplyr::n_distinct(out$well_num), nrow(out)))
  out
}

# GIN/NGWD: requires a manual download of the Nechako subset to cfg$gin_cache.
load_gin <- function() {
  if (!cfg$src_gin_ngwd) return(NULL)
  cat("Loading GIN/NGWD (NRCan) ... [STUB - requires manual cache]\n")
  if (!file.exists(cfg$gin_cache)) {
    cat(sprintf("  [manual] Download GIN subset for the Nechako bbox to %s\n", cfg$gin_cache))
    return(NULL)
  }
  NULL   # implementation pending cache availability
}
load_wsc <- function() {
  if (!cfg$src_wsc) return(NULL)
  cat("Loading WSC groundwater (tidyhydat) ... [STUB]\n")
  NULL
}
load_bc_realtime <- function() {
  if (!cfg$src_bc_realtime) return(NULL)
  cat("Loading BC Real-Time (ArcGIS) ... [STUB]\n")
  NULL
}
load_gwells <- function() {
  if (!cfg$src_gwells) return(NULL)
  cat("Loading GWELLS API (BC) ... [STUB]\n")
  NULL
}

#---- 5. Load, aggregate to well-month, spatial filter -----------------------
cat("\n=== Loading groundwater time series ===\n")
sources <- list(PGOWN = load_pgown(), GIN = load_gin(), WSC = load_wsc(),
                BC_RT = load_bc_realtime(), GWELLS = load_gwells())
sources_active <- names(purrr::compact(sources))
ts_all <- dplyr::bind_rows(purrr::compact(sources))
if (!nrow(ts_all)) stop("No groundwater data loaded from any source.")

# One value per well, month and source (median of any sub-monthly records),
# then keep the highest-priority source for each well-month.
source_priority <- c("PGOWN", "GIN", "WSC", "BC_RT", "GWELLS")
ts_dedup <- ts_all |>
  dplyr::group_by(well_num, date, source_id) |>
  dplyr::summarise(gw = stats::median(gw, na.rm = TRUE),
                   well_name = dplyr::first(well_name), lat = dplyr::first(lat),
                   lon = dplyr::first(lon), .groups = "drop") |>
  dplyr::mutate(src_rank = match(source_id, source_priority)) |>
  dplyr::group_by(well_num, date) |>
  dplyr::slice_min(src_rank, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::select(-src_rank) |>
  dplyr::filter(date >= cfg$analysis_start, date <= cfg$analysis_end)
cat(sprintf("Combined: %d wells, %d well-months inside %s to %s\n",
            dplyr::n_distinct(ts_dedup$well_num), nrow(ts_dedup),
            format(cfg$analysis_start, "%Y-%m"), format(cfg$analysis_end, "%Y-%m")))
if (!nrow(ts_dedup)) stop("No groundwater records inside the framework analysis window.")

# Spatial filter against the 7-subbasin union (same polygon as the ERA5-Land mask).
wells_xy <- dplyr::distinct(ts_dedup, well_num, well_name, lat, lon, source_id)
wells_v  <- terra::vect(wells_xy, geom = c("lon", "lat"), crs = "EPSG:4326", keepgeom = TRUE)
in_basin <- as.vector(terra::is.related(wells_v, basin_ll, "intersects"))
wells_xy$in_basin <- in_basin
wells_xy$subbasin <- NA_character_
if (!is.na(subb_name_col) && any(in_basin)) {
  hit <- tryCatch(terra::extract(subb_ll, wells_v[in_basin]), error = function(e) NULL)
  if (!is.null(hit) && subb_name_col %in% names(hit)) {
    hit <- hit[!duplicated(hit$id.y), , drop = FALSE]
    wells_xy$subbasin[which(in_basin)[hit$id.y]] <- as.character(hit[[subb_name_col]])
  }
}
cat(sprintf("  Distinct candidate wells: %d | inside basin: %d\n", nrow(wells_xy), sum(in_basin)))

# Diagnostic: basin polygon + every candidate well.
p_extent <- ggplot2::ggplot() +
  {if (!is.null(subb_ll))
    ggplot2::geom_sf(data = sf::st_as_sf(subb_ll), fill = "grey85", colour = "grey55", alpha = 0.5)} +
  ggplot2::geom_sf(data = sf::st_as_sf(basin_ll), fill = NA, colour = "grey20", linewidth = 0.6) +
  ggplot2::geom_point(data = wells_xy, ggplot2::aes(lon, lat,
                                                    colour = ifelse(in_basin, "Inside basin", "Outside basin")), size = 1.6) +
  ggplot2::scale_colour_manual(values = c("Inside basin" = "#4575b4", "Outside basin" = "#d73027")) +
  # PGOWN is province-wide: show the basin plus a margin, so nearby wells that
  # fall just outside the polygon are visible without zooming out to all of BC.
  ggplot2::coord_sf(xlim = bb[1:2] + c(-1.5, 1.5), ylim = bb[3:4] + c(-1.0, 1.0)) +
  ggplot2::labs(title = "Basin extent check: candidate wells vs. 7-subbasin union",
                subtitle = sprintf("%d of %d candidate wells fall inside the basin", sum(in_basin), nrow(wells_xy)),
                colour = NULL, x = NULL, y = NULL) +
  ggplot2::theme_bw(base_size = 10)
ggplot2::ggsave(file.path(FIG_DIR, "sgi_basin_extent_check.png"), p_extent,
                width = 8, height = 7, dpi = 150, bg = "white")
cat("  Diagnostic plot saved: figures/sgi_basin_extent_check.png\n")

if (!any(in_basin)) stop("No groundwater wells fall inside the 7-subbasin union.")
ts_basin <- dplyr::filter(ts_dedup, well_num %in% wells_xy$well_num[in_basin])

#---- 6. Well quality filters -------------------------------------------------
month_index <- function(d) as.integer(format(d, "%Y")) * 12L + as.integer(format(d, "%m"))
well_stats <- ts_basin |>
  dplyr::group_by(well_num) |>
  dplyr::summarise(
    n_months     = dplyr::n(),
    date_min     = min(date), date_max = max(date),
    record_years = (month_index(max(date)) - month_index(min(date)) + 1L) / 12,
    n_expected   = month_index(max(date)) - month_index(min(date)) + 1L,
    gap_frac     = 1 - n_months / n_expected,
    sources_used = paste(sort(unique(source_id)), collapse = ","),
    .groups = "drop"
  )
is_pgown <- grepl("PGOWN", well_stats$sources_used)
well_stats$pass <- well_stats$record_years >= ifelse(is_pgown, cfg$min_years, pmax(5, cfg$min_years * 0.5)) &
  well_stats$gap_frac <= cfg$max_gap_frac * ifelse(is_pgown, 1, 1.5)
wells_ok <- dplyr::filter(well_stats, pass)
cat(sprintf("Wells after quality filters: %d of %d\n", nrow(wells_ok), nrow(well_stats)))
if (!nrow(wells_ok)) stop("No wells pass the filters. Relax cfg$min_years or cfg$max_gap_frac.")

ts_ok <- dplyr::filter(ts_basin, well_num %in% wells_ok$well_num)

#---- 7. Compute SGI per well x calendar month (whole record) ----------------
cat("\nComputing SGI (Gringorten normal scores per well x calendar month) ...\n")
sgi_long <- ts_ok |>
  dplyr::mutate(cal_month = as.integer(format(date, "%m"))) |>
  dplyr::group_by(well_num, cal_month) |>
  dplyr::mutate(n_cal = sum(is.finite(gw)), sgi = sgi_normal_scores(gw)) |>
  dplyr::ungroup() |>
  dplyr::arrange(well_num, date)

n_sgi <- sum(is.finite(sgi_long$sgi))
cat(sprintf("  Finite SGI values: %d of %d well-months (%.1f%%)\n",
            n_sgi, nrow(sgi_long), 100 * n_sgi / nrow(sgi_long)))
if (!n_sgi) stop("No well/calendar month has >= ", cfg$min_years_per_month,
                 " years. Lower cfg$min_years_per_month (framework minimum 5).")

# Plotting-position floors. With Gringorten the lowest SGI a sample of n years
# can take is qnorm(0.56 / (n + 0.12)). The SHORTEST well/month sample (n_min)
# saturates first (least negative floor): below this line the shortest wells
# can no longer get drier, so the basin median is carried by the longer
# records. The LONGEST sample (n_max) sets the lowest SGI attainable anywhere.
n_used <- sgi_long$n_cal[is.finite(sgi_long$sgi)]
n_min  <- min(n_used); n_max <- max(n_used)
sgi_clamp_floor <- max(-cfg$z_clip, stats::qnorm(0.56 / (n_min + 0.12)))
sgi_floor_nmax  <- max(-cfg$z_clip, stats::qnorm(0.56 / (n_max + 0.12)))
sgi_clamp_label <- sprintf("shortest-record floor (n=%d \u2192 SGI=%.3f)", n_min, sgi_clamp_floor)
cat(sprintf("  Plotting-position floors: shortest sample n=%d -> %.3f | longest n=%d -> %.3f\n",
            n_min, sgi_clamp_floor, n_max, sgi_floor_nmax))

#---- 8. Basin aggregate: median across wells --------------------------------
sgi_basin <- sgi_long |>
  dplyr::filter(is.finite(sgi)) |>
  dplyr::group_by(date) |>
  dplyr::summarise(
    sgi_median      = stats::median(sgi),
    sgi_mean        = mean(sgi),
    sgi_sd          = if (dplyr::n() > 1) stats::sd(sgi) else NA_real_,
    n_wells         = dplyr::n(),
    sources_contrib = paste(sort(unique(source_id)), collapse = ","),
    .groups = "drop"
  ) |>
  dplyr::mutate(sgi_median = ifelse(n_wells >= cfg$min_wells_per_month, sgi_median, NA_real_))

sgi_basin <- complete_monthly(sgi_basin) |>
  dplyr::mutate(month_label = format(date, "%Y-%m"),
                n_wells     = ifelse(is.na(n_wells), 0L, n_wells))

ev_res <- sgi_event_table(sgi_basin$date, sgi_basin$sgi_median, sgi_basin$n_wells)
sgi_basin <- sgi_basin |>
  dplyr::mutate(drought_cat      = sgi_class(sgi_median),
                in_drought_event = ev_res$in_event) |>
  dplyr::select(date, month_label, sgi_median, sgi_mean, sgi_sd, n_wells,
                sources_contrib, drought_cat, in_drought_event)

#---- 8b. Sign-convention cross-check against a known drought year -----------
# Summer-fall 2023 was a severe, widely documented BC drought. With the sign
# set correctly the basin SGI should be negative in Aug-Oct 2023.
win <- dplyr::filter(sgi_basin, month_label %in% c("2023-08", "2023-09", "2023-10"))
mean_2023 <- if (nrow(win)) mean(win$sgi_median, na.rm = TRUE) else NA_real_
if (is.finite(mean_2023)) {
  cat(sprintf("\n  Known-drought cross-check: mean SGI Aug-Oct 2023 = %.2f\n", mean_2023))
  if (mean_2023 > 0.5) {
    cat(strrep("!", 70), "\n")
    cat("  WARNING: Aug-Oct 2023 is a documented severe BC drought but SGI is strongly\n")
    cat("  POSITIVE here - independent evidence that cfg$invert_gw_sign may be backwards.\n")
    cat(strrep("!", 70), "\n")
  } else if (mean_2023 < 0) {
    cat("  Consistent with the 2023 drought (negative SGI) - supports the sign convention.\n")
  }
} else cat("\n  Known-drought cross-check skipped: no Aug-Oct 2023 basin SGI.\n")

#---- 9. Outputs --------------------------------------------------------------
basin_csv <- file.path(OUT_DIR, cfg$basin_file)
readr::write_csv(sgi_basin, basin_csv, na = "NA")

wide <- data.frame(Year = sort(unique(as.integer(format(sgi_basin$date, "%Y")))))
for (m in 1:12) {
  s <- sgi_basin[as.integer(format(sgi_basin$date, "%m")) == m, ]
  wide[[month.abb[m]]] <- s$sgi_median[match(wide$Year, as.integer(format(s$date, "%Y")))]
}
utils::write.csv(wide, file.path(OUT_DIR, "sgi_01_basin_median_by_month.csv"),
                 row.names = FALSE, na = "")

readr::write_csv(sgi_long |> dplyr::select(well_num, date, gw, sgi, n_cal, source_id),
                 file.path(OUT_DIR, "sgi_well_timeseries.csv"), na = "NA")

selected_wells <- wells_ok |>
  dplyr::left_join(dplyr::distinct(wells_xy, well_num, .keep_all = TRUE) |>
                     dplyr::select(well_num, well_name, lat, lon, subbasin), by = "well_num") |>
  dplyr::select(well_num, well_name, lat, lon, subbasin, sources_used, date_min, date_max,
                record_years, n_months, gap_frac) |>
  dplyr::arrange(well_num)
readr::write_csv(selected_wells, file.path(OUT_DIR, "sgi_selected_wells.csv"), na = "")
readr::write_csv(ev_res$events, file.path(OUT_DIR, "sgi_drought_events.csv"), na = "")

cat("\n=== Qualifying wells ===\n")
print(as.data.frame(selected_wells[, c("well_num", "well_name", "subbasin", "record_years", "gap_frac")]),
      row.names = FALSE)

cat("\n========================================================\n")
cat("  SGI - Nechako River Basin (multi-source)\n")
cat(sprintf("  Sources active: %s\n", paste(sources_active, collapse = ", ")))
fin <- sgi_basin[is.finite(sgi_basin$sgi_median), ]
cat(sprintf("  Wells: %d | Basin months with SGI: %d | %s to %s\n",
            nrow(selected_wells), nrow(fin), min(fin$month_label), max(fin$month_label)))
cat(sprintf("  Drought events: %d (onset < %.2f, end >= %.2f, min %d months)\n",
            nrow(ev_res$events), DROUGHT_ONSET, DROUGHT_END, ev_res$min_duration))
cat("========================================================\n\n")
show <- if (!is.null(cfg$target_months)) dplyr::filter(sgi_basin, month_label %in% cfg$target_months) else utils::tail(fin, 24)
print(as.data.frame(show[, c("month_label", "sgi_median", "n_wells", "drought_cat", "in_drought_event")]),
      row.names = FALSE)

#---- 10. Figures -------------------------------------------------------------
plot_index_timeseries(
  sgi_basin$date, sgi_basin$sgi_median, "SGI", "SGI (Groundwater, Well-Network Median)",
  file.path(FIG_DIR, "sgi_basin_timeseries.png"),
  subtitle_extra = sprintf("%d wells | %s", nrow(selected_wells), paste(sources_active, collapse = "+")),
  clamp_floor = sgi_clamp_floor, clamp_label = sgi_clamp_label
)

cat("\n=== Saving individual well SGI plots ===\n")
for (w in unique(sgi_long$well_num)) {
  dw <- complete_monthly(dplyr::filter(sgi_long, well_num == w) |> dplyr::select(date, sgi))
  if (!any(is.finite(dw$sgi))) next
  wn <- selected_wells$well_name[selected_wells$well_num == w][1]
  plot_index_timeseries(
    dw$date, dw$sgi, "SGI",
    sprintf("SGI - Well %s%s", w, if (!is.na(wn) && nzchar(wn)) paste0(" (", wn, ")") else ""),
    file.path(WELL_FIG_DIR, sprintf("sgi_well_%s.png", gsub("[^A-Za-z0-9_-]+", "_", w))),
    subtitle_extra = paste("source:", paste(unique(sgi_long$source_id[sgi_long$well_num == w]), collapse = ",")),
    width = 10, height = 4.5
  )
}

#---- 11. Summary file --------------------------------------------------------
summary_file <- file.path(OUT_DIR, "sgi_summary.txt")
z <- sgi_basin$sgi_median
lines <- c(
  "============================================================",
  "SGI SUMMARY - Nechako River Basin",
  "============================================================",
  sprintf("Run date            : %s", format(Sys.time(), "%Y-%m-%d %H:%M")),
  "Method              : Bloomfield & Marchant (2013) normal scores, per well x calendar month",
  "Plotting position   : Gringorten (r - 0.44)/(n + 0.12); zero variance -> 0",
  sprintf("Clamp               : +/-%.2f", cfg$z_clip),
  sprintf("Reference           : whole record inside %s to %s",
          format(cfg$analysis_start, "%Y-%m"), format(cfg$analysis_end, "%Y-%m")),
  sprintf("Min years/cal month : %d", cfg$min_years_per_month),
  sprintf("Well filters        : record >= %g yr (PGOWN), gap fraction <= %.2f",
          cfg$min_years, cfg$max_gap_frac),
  sprintf("Basin boundary      : union of 7 subbasins (%.0f km^2)", basin_area_km2),
  sprintf("Sign convention     : invert_gw_sign = %s", cfg$invert_gw_sign),
  sprintf("Sources active      : %s", paste(sources_active, collapse = ", ")),
  sprintf("Wells               : %d candidate | %d in basin | %d qualifying",
          nrow(wells_xy), sum(in_basin), nrow(selected_wells)),
  sprintf("Basin aggregation   : median across wells (min %d wells/month)", cfg$min_wells_per_month),
  sprintf("Basin months        : %d with SGI, %s to %s", nrow(fin), min(fin$month_label), max(fin$month_label)),
  sprintf("Wells per month     : median %.0f, min %d, max %d",
          stats::median(fin$n_wells), min(fin$n_wells), max(fin$n_wells)),
  sprintf("Plotting floors     : shortest sample n=%d -> %.3f | longest n=%d -> %.3f",
          n_min, sgi_clamp_floor, n_max, sgi_floor_nmax),
  sprintf("Known drought 2023  : mean Aug-Oct SGI = %s",
          if (is.finite(mean_2023)) sprintf("%.2f", mean_2023) else "n/a"),
  "",
  sprintf("Drought events      : %d (identify_drought_events_ts, scale %d; onset < %.2f, end >= %.2f, min %d months)",
          nrow(ev_res$events), SGI_EVENT_SCALE, DROUGHT_ONSET, DROUGHT_END, ev_res$min_duration),
  "",
  "Class frequency (basin median, McKee et al. 1993 bins):",
  sprintf("  W3 Extremely wet  (>= 2.0)       : %.1f%%", 100 * mean(z >= 2, na.rm = TRUE)),
  sprintf("  W2 Severely wet   [1.5, 2)       : %.1f%%", 100 * mean(z >= 1.5 & z < 2, na.rm = TRUE)),
  sprintf("  W1 Moderately wet [1, 1.5)       : %.1f%%", 100 * mean(z >= 1 & z < 1.5, na.rm = TRUE)),
  sprintf("  N0 Near normal    (-1, 1)        : %.1f%%", 100 * mean(z > -1 & z < 1, na.rm = TRUE)),
  sprintf("  D1 Moderately dry (-1.5, -1]     : %.1f%%", 100 * mean(z <= -1 & z > -1.5, na.rm = TRUE)),
  sprintf("  D2 Severely dry   (-2, -1.5]     : %.1f%%", 100 * mean(z <= -1.5 & z > -2, na.rm = TRUE)),
  sprintf("  D3 Extremely dry  (<= -2.0)      : %.1f%%", 100 * mean(z <= -2, na.rm = TRUE)),
  "",
  "NOTE: the basin value is a well-network MEDIAN (point support), not an",
  "area-weighted ERA5-Land basin mean; the number of contributing wells changes",
  "through time (see n_wells)."
)
writeLines(lines, summary_file)

cat("\n============================================================\n")
cat("SGI CALCULATION COMPLETE\n")
cat("============================================================\n")
cat(sprintf("Output directory: %s\n", normalizePath(OUT_DIR)))
cat(sprintf("Basin series    : %s\n", basin_csv))
cat(sprintf("Summary         : %s\n", summary_file))