####################################################################################
# DROUGHT_ANALYSIS_utils.R  ·  SHARED UTILITIES FOR SPI, SPEI, AND SWEI
####################################################################################
## A. PROJECT CONFIGURATION ───────────────────────────────────────────
WD_PATH <- Sys.getenv("NECHAKO_WD", "D:/Nechako_Drought/Nechako/")
BASIN_SHP       <- file.path(WD_PATH, "Spatial/nechakoBound_dissolve.kmz")

## Private helper: unzip a KMZ file into a self-cleaning temp directory and
## call FUN(kml_path), returning FUN's result.  The temp directory is deleted
## by on.exit() when this function returns — so FUN must finish reading the
## KML before .read_kmz() exits (which it always does, since FUN is called
## synchronously inside the function body).
.read_kmz <- function(kmz_path, FUN) {
  tmp <- tempfile()
  dir.create(tmp, showWarnings = FALSE)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  utils::unzip(kmz_path, exdir = tmp)
  kml <- list.files(tmp, pattern = "\\.kml$", full.names = TRUE, recursive = TRUE)
  if (!length(kml))
    stop(".read_kmz: no .kml file found inside ", basename(kmz_path))
  FUN(kml[1])
}

## Helper: load basin as a single dissolved SpatVector (terra).
## terra::vect() cannot read KMZ directly on Windows (missing GDAL KMZ driver).
## .read_kmz() handles the unzip/cleanup; terra::vect() reads the inner KML.
load_basin_vect <- function(path = BASIN_SHP) {
  if (!file.exists(path)) stop("Basin file not found: ", path)
  v <- if (tolower(tools::file_ext(path)) == "kmz")
    .read_kmz(path, terra::vect)
  else
    terra::vect(path)
  if (nrow(v) > 1L) v <- terra::aggregate(v)
  v
}

## Index-specific directories
SPI_SEAS_DIR     <- file.path(WD_PATH, "spi_results_seasonal/")
SPEI_SEAS_DIR    <- file.path(WD_PATH, "spei_results_seasonal/")
SPEI_THW_SEAS_DIR <- file.path(WD_PATH, "spei_results_seasonal_thw/")   # Thornthwaite PET
SWEI_SEAS_DIR    <- file.path(WD_PATH, "swei_results_seasonal/")

## Output directories
TREND_DIR        <- file.path(WD_PATH,  "temporal_drought/")
BASIN_PLOT_DIR   <- file.path(WD_PATH,  "basin_averaged_plots/")
POINT_PLOT_DIR   <- file.path(WD_PATH,  "point_timeseries_plots/")
GIF_DIR <- file.path(TREND_DIR, "drought_gifs/")
CACHE_DIR <- file.path(WD_PATH, "temporal_drought/cache/")
TELE_DIR <- file.path(WD_PATH, "teleconnections")

## Figure output constants (used by save_figure() and all w* scripts)
FIG_DPI         <- 300L   # PNG resolution
FIG_WIDTH_WIDE  <- 14     # wide panel figures (inches)
FIG_WIDTH_STD   <- 10     # standard single-panel figures (inches)
FIG_HEIGHT_STD  <-  8     # standard figure height (inches)

## CRS and thresholds
EQUAL_AREA_CRS  <- "EPSG:3005"
DROUGHT_ONSET      <- -1.0
DROUGHT_END        <-  0.0
SEVERE_THRESHOLD   <- -1.3
EXTREME_THRESHOLD  <- -1.6
DROUGHT_THRESHOLD  <- -1.0

## Timescales
SPI_SCALES          <- c(1, 2, 3, 6, 12)
SPEI_SCALES         <- c(1, 2, 3, 6, 12)
SWEI_SCALE          <- 3  # Single timescale only
INDICES             <- c("spi", "spei", "spei_thw", "swei")
TIMESCALES_STANDARD <- SPI_SCALES          # used by w4 for time series plots
TIMESCALES_SPATIAL  <- SPI_SCALES          # used by w4 for spatial figures
MONTH_NAMES    <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

## B. BASIN-AVERAGED CSV LOADER ─────────────────────────────────────
#' Load a pre-computed basin-averaged drought index CSV and return a tidy
#' long-format data.frame(date, value).
#'
#' AUTHORITATIVE SOURCE: these CSVs are produced by w1_trend_test.R via
#' save_basin_avg_from_pixels() (see Section B2 below).  The methodology is:
#'   (1) compute the drought index at every basin pixel from per-pixel inputs,
#'   (2) take the area-weighted spatial mean of those pixel-level index values.
#' This is the standard hydrometeorological approach: area-weighted averaging
#' of standardised index values that were calibrated on the same spatial basis.
#' Run w1_trend_test.R before w1_basin_timeseries.R.
#'
#' File naming convention: {index}_{scale:02d}_basin_averaged_by_month.csv
#' File structure:
#'   Year | Jan | Feb | Mar | … | Dec
#'   1950 | -1.88 |  NA | …
#'
#' @param index_type  "spi", "spei", or "swei"
#' @param scale       Integer timescale (e.g. 1, 3, 6, 12)
#' @param data_dir    Directory containing the CSV; if NULL, the index-specific
#'                    seasonal directory from section A is used automatically.
#' @return data.frame with columns \code{date} (Date) and \code{value} (numeric),
#'         sorted by date with NA rows removed.  Returns NULL with a warning if
#'         the file is not found (usually means w3 has not been run yet).
load_basin_avg_csv <- function(index_type, scale, data_dir = NULL) {
  dir_map <- list(spi      = SPI_SEAS_DIR,
                  spei     = SPEI_SEAS_DIR,
                  spei_thw = SPEI_THW_SEAS_DIR,
                  swei     = SWEI_SEAS_DIR)
  dir     <- if (!is.null(data_dir)) data_dir else {
    d <- dir_map[[tolower(index_type)]]
    if (is.null(d)) stop("Unknown index_type: ", index_type)
    d
  }
  f <- file.path(dir, sprintf("%s_%02d_basin_averaged_by_month.csv",
                              tolower(index_type), as.integer(scale)))
  if (!file.exists(f)) {
    warning(sprintf("Basin-averaged CSV not found: %s", f))
    return(NULL)
  }
  df <- as.data.frame(data.table::fread(f))  # faster than read.csv for wide CSVs
  if (!"Year" %in% names(df))
    stop("Expected 'Year' column in ", basename(f))
  mon_cols <- intersect(MONTH_NAMES, names(df))
  if (!length(mon_cols))
    stop("No month name columns found in ", basename(f),
         ".  Expected: ", paste(MONTH_NAMES, collapse = ", "))
  long <- do.call(rbind, lapply(seq_along(mon_cols), function(mi) {
    data.frame(
      date  = as.Date(paste(df$Year, mi, "01", sep = "-")),
      value = as.numeric(df[[mon_cols[mi]]]),
      stringsAsFactors = FALSE)
  }))
  long <- long[order(long$date), ]
  long <- long[!is.na(long$value) & is.finite(long$value), ]
  rownames(long) <- NULL
  long
}

## B2. BASIN-AVERAGED CSV WRITER ────────────────────────────────────
#' Compute an area-weighted basin-average time series from a pixel × month
#' matrix and write the result as the standard basin_averaged_by_month.csv
#' expected by load_basin_avg_csv().
#'
#' Called by w1_trend_test.R::process_index() after the pixel time-series
#' matrix has been assembled.  The area weights are the equal-area cell sizes
#' extracted from the raster template (already in EPSG:3005).
#'
#' Methodology:
#'   basin_avg[t] = Σ( index_value[pixel, t] × cell_area[pixel] )
#'                / Σ( cell_area[pixel] )          (sum over non-NA pixels)
#'
#' @param ts_matrix   Numeric matrix (n_valid_pixels × n_months).  Rows
#'                    correspond to valid basin pixels; columns to calendar
#'                    months in chronological order.
#' @param area_weights Numeric vector of length n_valid_pixels; cell areas in
#'                    m² in the equal-area CRS (from terra::cellSize on the
#'                    BC Albers projected raster template).
#' @param dates       Date vector of length n_months giving the first day of
#'                    each month represented by a column of ts_matrix.
#' @param index_type  "spi", "spei", "swei", or any other index key.
#' @param scale       Integer timescale (e.g. 1, 3, 6, 12).
#' @param out_dir     Directory to write the CSV.  Should be the index-specific
#'                    seasonal results directory (SPI_SEAS_DIR etc.) so that
#'                    load_basin_avg_csv() can find it without specifying data_dir.
#' @return Invisible data.frame (Year × month-name columns) as written to disk.
save_basin_avg_from_pixels <- function(ts_matrix, area_weights, dates,
                                       index_type, scale, out_dir) {
  ## ── Input validation ───────────────────────────────────────────────
  if (!is.matrix(ts_matrix))
    stop("ts_matrix must be a numeric matrix (pixels × months)")
  if (nrow(ts_matrix) != length(area_weights))
    stop(sprintf(
      "ts_matrix has %d rows but area_weights has %d elements",
      nrow(ts_matrix), length(area_weights)))
  if (ncol(ts_matrix) != length(dates))
    stop(sprintf(
      "ts_matrix has %d columns but dates has %d elements",
      ncol(ts_matrix), length(dates)))
  
  ## ── Guard: replace zero / negative / NA areas with column median ───
  bad_w <- is.na(area_weights) | !is.finite(area_weights) | area_weights <= 0
  if (any(bad_w)) {
    med_w <- median(area_weights[!bad_w], na.rm = TRUE)
    area_weights[bad_w] <- if (is.finite(med_w) && med_w > 0) med_w else 1
    warning(sprintf(
      "save_basin_avg_from_pixels: %d bad area weight(s) replaced with %.0f m²",
      sum(bad_w), area_weights[which(bad_w)[1]]))
  }
  
  ## ── Area-weighted mean for every time step ─────────────────────────
  basin_avg <- vapply(seq_len(ncol(ts_matrix)), function(t) {
    v  <- ts_matrix[, t]
    ok <- !is.na(v) & is.finite(v)
    if (!any(ok)) return(NA_real_)
    sum(v[ok] * area_weights[ok]) / sum(area_weights[ok])
  }, numeric(1L))
  
  ## ── Reshape to Year × Month wide format ────────────────────────────
  yrs <- as.integer(format(dates, "%Y"))
  mos <- as.integer(format(dates, "%m"))
  all_years <- sort(unique(yrs))
  
  df <- data.frame(Year = all_years, stringsAsFactors = FALSE)
  for (mi in seq_along(MONTH_NAMES)) {
    hits    <- which(mos == mi)
    row_idx <- match(yrs[hits], all_years)        # vectorised: all hits at once
    valid   <- !is.na(row_idx)
    col_vals <- rep(NA_real_, length(all_years))
    col_vals[row_idx[valid]] <- basin_avg[hits[valid]]
    df[[MONTH_NAMES[mi]]] <- col_vals
  }
  
  ## ── Write CSV ───────────────────────────────────────────────────────
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_file <- file.path(
    out_dir,
    sprintf("%s_%02d_basin_averaged_by_month.csv",
            tolower(index_type), as.integer(scale)))
  utils::write.csv(df, out_file, row.names = FALSE)
  cat(sprintf(
    "  ✓ Basin-avg CSV written (area-weighted, %d pixels, %d months): %s\n",
    nrow(ts_matrix), ncol(ts_matrix), basename(out_file)))
  invisible(df)
}

## B3. PACKAGE BOOTSTRAP ──────────────────────────────────────────────
utils_load_packages <- function(pkgs) {
  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing)) {
    message("Installing missing packages: ", paste(missing, collapse = ", "))
    install.packages(missing, repos = "https://cloud.r-project.org/")
  }
  invisible(lapply(pkgs, library, character.only = TRUE, quietly = TRUE))
  cat("✓ Packages loaded:", paste(pkgs, collapse = ", "), "\n")
}

####################################################################################
# C. NETCDF FILE FINDERS ────────────────────────────────────────────
####################################################################################
find_seasonal_nc_files <- function(data_dir, index_type, scale) {
  idx <- tolower(index_type)
  pattern <- sprintf("^%s_%02d_month\\d{2}_[A-Za-z]+\\.nc$", idx, scale)
  files <- list.files(data_dir, pattern = pattern,
                      full.names = TRUE, ignore.case = TRUE, recursive = TRUE)
  if (!length(files)) {
    warning(sprintf("No seasonal files found for %s-%d in: %s",
                    toupper(index_type), scale, data_dir))
    return(character(0))
  }
  mon_nums <- as.integer(regmatches(basename(files),
                                    regexpr("(?<=month)\\d{2}", basename(files), perl = TRUE)))
  files[order(mon_nums)]
}

find_swei_seasonal_files <- function(data_dir, index_type, scale) {  # ← CHANGED
  pattern <- sprintf("^swei_%02d_month\\d{2}_[A-Za-z]+\\.nc$", scale)
  files <- list.files(data_dir, pattern = pattern,
                      full.names = TRUE, ignore.case = TRUE, recursive = TRUE)
  if (!length(files)) {
    warning(sprintf("No SWEI seasonal files found for scale %d in: %s", scale, data_dir))
    return(character(0))
  }
  mon_nums <- as.integer(regmatches(basename(files),
                                    regexpr("(?<=month)\\d{2}", basename(files), perl = TRUE)))
  files[order(mon_nums)]
}

## D. NETCDF DATE PARSING ────────────────────────────────────────────
extract_dates_from_nc <- function(nc_file, n_layers = NULL) {
  r_temp <- terra::rast(nc_file)
  n_layers <- if (is.null(n_layers)) terra::nlyr(r_temp) else n_layers
  dates  <- tryCatch({
    d  <- terra::time(r_temp)
    if (!is.null(d) && length(d) > 0) {
      d_conv  <- tryCatch(as.Date(d), error = function(e) as.Date(d, origin = "1970-01-01"))
      if (all(!is.na(d_conv)) &&
          min(d_conv) >= as.Date("1900-01-01") &&
          max(d_conv) <= as.Date("2030-01-01")) {
        return(.pad_dates(d_conv, n_layers))
      }
    }
    NULL
  }, error = function(e) NULL)
  if (is.null(dates)) {
    cat("  [dates] Using fallback monthly sequence (Jan 1950)\n")
    d <- seq(as.Date("1950-01-01"), by = "month", length.out = n_layers)
    if (tail(d, 1) > as.Date("2025-12-31")) {
      end <- as.Date("2025-12-01")
      start <- end - lubridate::months(n_layers - 1)
      d <- seq(start, by = "month", length.out = n_layers)
    }
    d
  } else dates
}

.pad_dates <- function(d, n) {
  if (length(d) == n) return(d)
  if (length(d) > n) return(d[1:n])
  extra <- seq(tail(d, 1) + lubridate::months(1), by = "month", length.out = n - length(d))
  c(d, extra)
}

## E. BASIN SPATIAL HELPERS ──────────────────────────────────────────
precompute_basin_geometry  <- function(basin_shp_path, raster_file_path, crs = EQUAL_AREA_CRS) {
  basin_v  <- terra::vect(basin_shp_path)
  r_tmpl  <- terra::rast(raster_file_path)[[1]]
  r_proj  <- terra::project(r_tmpl, crs)
  b_proj  <- terra::project(basin_v, crs)
  total_area  <- sum(terra::expanse(b_proj, unit = "m"))
  cell_areas  <- terra::cellSize(r_proj, unit = "m")
  areas_masked  <- terra::mask(cell_areas, b_proj)
  basin_rast  <- terra::rasterize(b_proj, r_proj, cover = TRUE)
  area_v  <- as.vector(terra::values(areas_masked, na.rm = FALSE))
  cover_v  <- as.vector(terra::values(basin_rast, na.rm = FALSE))
  eff_areas  <- area_v * cover_v
  valid_idx  <- which(!is.na(cover_v) & cover_v > 0)
  cat(sprintf("  Basin area: %.1f km²  |  Valid cells: %d\n",
              total_area / 1e6, length(valid_idx)))
  list(total_basin_area = total_area,
       effective_areas = eff_areas,
       valid_cell_idx = valid_idx,
       equal_area_crs = crs,
       basin_shp_path = basin_shp_path,
       raster_file_path = raster_file_path)
}

compute_layer_statistics  <- function(layer_idx, raster_file, basin_geom, threshold = DROUGHT_THRESHOLD) {
  r  <- terra::rast(raster_file)[[layer_idx]]
  r  <- terra::project(r, basin_geom$equal_area_crs)
  v  <- as.vector(terra::values(r, na.rm = FALSE))
  vc  <- v[basin_geom$valid_cell_idx]
  ac  <- basin_geom$effective_areas[basin_geom$valid_cell_idx]
  ok  <- !is.na(vc) & is.finite(vc)
  if (!any(ok)) return(list(mean_value = NA_real_, valid_fraction = 0,
                            drought_fraction = NA_real_, drought_area = 0,
                            n_cells_valid = 0, n_cells_drought = 0))
  valid_area  <- sum(ac[ok])
  mean_val  <- sum(vc[ok] * ac[ok]) / valid_area
  valid_frac  <- valid_area / basin_geom$total_basin_area
  drought_ok  <- ok & (vc <= threshold)
  drought_area  <- sum(ac[drought_ok])
  list(mean_value = mean_val,
       valid_fraction = valid_frac,
       drought_fraction = drought_area / basin_geom$total_basin_area,
       drought_area = drought_area,
       n_cells_valid = sum(ok),
       n_cells_drought = sum(drought_ok))
}

## F. DROUGHT EVENT DETECTION ────────────────────────────────────────
detect_drought_events <- function(df,
                                  onset_threshold       = DROUGHT_ONSET,
                                  termination_threshold = DROUGHT_END,
                                  min_duration          = 1,
                                  severe_thr            = SEVERE_THRESHOLD,
                                  extreme_thr           = EXTREME_THRESHOLD) {
  stopifnot(all(c("date", "value") %in% names(df)))
  df    <- df[order(df$date), ]
  vals  <- df$value
  dates <- df$date
  
  ## Canonical empty frame — returned when there are no events
  empty_events <- data.frame(
    event_id        = integer(),
    start_date      = as.Date(character()),
    end_date        = as.Date(character()),
    duration_months = integer(),
    min_value       = numeric(),
    severity        = character(),
    stringsAsFactors = FALSE)
  
  classify_severity <- function(min_v) {
    if      (min_v <= extreme_thr) "Extreme"
    else if (min_v <= severe_thr)  "Severe"
    else                           "Moderate"
  }
  
  ## Pre-allocated list accumulator — O(1) append, avoids the O(n²) rbind()
  ## that the previous <<- closure caused by copying the growing data.frame
  ## on every event.  do.call(rbind, event_list) is called once at the end.
  event_list <- list()
  
  close_event <- function(s, e) {
    dur <- e - s + 1L
    if (dur < min_duration) return()
    mv <- min(vals[s:e], na.rm = TRUE)
    event_list[[length(event_list) + 1L]] <<- data.frame(
      event_id        = length(event_list) + 1L,
      start_date      = dates[s],
      end_date        = dates[e],
      duration_months = dur,
      min_value       = mv,
      severity        = classify_severity(mv),
      stringsAsFactors = FALSE)
  }
  
  in_drought <- FALSE
  start_idx  <- NA_integer_
  for (i in seq_along(vals)) {
    if (!in_drought && !is.na(vals[i]) && vals[i] < onset_threshold) {
      in_drought <- TRUE
      start_idx  <- i
    } else if (in_drought && !is.na(vals[i]) && vals[i] >= termination_threshold) {
      close_event(start_idx, i - 1L)
      in_drought <- FALSE
    }
  }
  if (in_drought) close_event(start_idx, length(vals))   # drought still open at end
  
  if (!length(event_list)) return(empty_events)
  do.call(rbind, event_list)
}

## F2. SPATIAL HELPERS (used by w4_trends_visualization.R) ──────────
#' Load basin shapefile and reproject to equal-area CRS.
#' KMZ files are handled via the shared .read_kmz() helper.
load_basin <- function(shp_path = BASIN_SHP, crs = EQUAL_AREA_CRS) {
  if (!file.exists(shp_path)) stop("Basin file not found: ", shp_path)
  b <- if (tolower(tools::file_ext(shp_path)) == "kmz")
    .read_kmz(shp_path, function(kml) sf::st_read(kml, quiet = TRUE))
  else
    sf::st_read(shp_path, quiet = TRUE)
  if (nrow(b) > 1L) b <- sf::st_as_sf(sf::st_union(b))
  if (sf::st_crs(b)$input != crs) b <- sf::st_transform(b, crs)
  cat(sprintf("✓ Basin loaded (%d polygon(s) → dissolved, CRS: %s)\n", nrow(b), crs))
  b
}

#' Clip a data.table of lon/lat points to basin using sf::st_intersects
clip_to_basin  <- function(dt, basin_sf, crs = EQUAL_AREA_CRS) {
  basin_v    <- terra::vect(basin_sf)
  basin_ext  <- terra::ext(basin_v)
  dt_pre     <- dt[lon >= basin_ext[1] & lon <= basin_ext[2] &
                     lat >= basin_ext[3] & lat <= basin_ext[4]]
  if (nrow(dt_pre) == 0) { warning("No points within basin extent"); return(dt_pre) }
  pts_sf  <- sf::st_as_sf(terra::vect(dt_pre, geom = c("lon", "lat"), crs = crs))
  idx     <- sapply(sf::st_intersects(pts_sf, basin_sf, sparse = TRUE), length) > 0
  dt_pre[idx, ]
}

#' Create a raster template aligned to basin extent
create_raster_template  <- function(data_dt, basin_sf) {
  data_dt  <- data_dt[!is.na(lon) & !is.na(lat)]
  if (nrow(data_dt) == 0) return(NULL)
  u_lons  <- unique(data_dt$lon); u_lats  <- unique(data_dt$lat)
  if (length(u_lons) >= 2 && length(u_lats) >= 2) {
    dx   <- median(diff(sort(u_lons)), na.rm = TRUE)
    dy   <- median(diff(sort(u_lats)), na.rm = TRUE)
    res  <- min(dx, dy, na.rm = TRUE)
  } else {
    ext  <- terra::ext(terra::vect(basin_sf))
    res  <- min(ext[2]-ext[1], ext[4]-ext[3]) / 50
  }
  if (!is.finite(res) || res <= 0) res  <- 5000
  terra::rast(ext = terra::ext(terra::vect(basin_sf)),
              resolution = res, crs = EQUAL_AREA_CRS)
}

#' Rasterise a column from a data.table onto a template and mask to basin
create_raster_from_points <- function(dt, template, value_col, basin_sf) {
  if (is.null(template)) return(NULL)
  empty <- function() {
    r <- terra::rast(template); terra::values(r) <- NA
    terra::mask(r, terra::vect(basin_sf))
  }
  if (nrow(dt) == 0 || !value_col %in% names(dt)) return(empty())
  dt <- dt[!is.na(lon) & !is.na(lat)]
  if (nrow(dt) == 0 || all(is.na(dt[[value_col]]))) return(empty())
  pts <- terra::vect(dt, geom = c("lon","lat"), crs = terra::crs(template))
  r   <- terra::rasterize(pts, template, field = value_col, touches = TRUE)
  terra::mask(r, terra::vect(basin_sf))
}

#' 3×3 focal mean smoother (passes NA through gracefully)
smooth_raster <- function(r) {
  if (is.null(r) || all(is.na(terra::values(r)))) return(r)
  terra::focal(r, w = matrix(1, 3, 3), fun = mean, na.rm = TRUE)
}

#' Plot a raster with basin overlay; handles all-NA gracefully.
#' Relies on basin_sf_global being set in the calling script.
plot_raster_clean <- function(r, main, zlim = NULL, col, breaks = NULL,
                              legend = TRUE, categorical = FALSE, legend_title = "") {
  if (is.null(r) || all(is.na(terra::values(r)))) {
    plot.new(); text(0.5, 0.5, "No valid data", cex = 1.2, col = "red")
    return(invisible())
  }
  if (!categorical) r <- smooth_raster(r)
  plg <- list(cex = 0.9)
  if (nchar(legend_title)) plg$title <- legend_title
  if (!is.null(breaks)) {
    terra::plot(r, main = main, cex.main = 1.0, col = col, breaks = breaks,
                axes = FALSE, box = FALSE, legend = legend, plg = plg, colNA = NA)
  } else {
    terra::plot(r, main = main, cex.main = 1.0, col = col, zlim = zlim,
                axes = FALSE, box = FALSE, legend = legend, plg = plg, colNA = NA)
  }
  plot(sf::st_geometry(basin_sf_global),
       add = TRUE, col = NA, border = "black", lwd = 2.5)
}

## G. GGPLOT HELPERS ─────────────────────────────────────────────────
## Helper: make one annotation_custom rect band (avoids Date-scale numeric warning)
.band_rect <- function(ymin, ymax, fill, alpha = 0.10) {
  ggplot2::annotation_custom(
    grob = grid::rectGrob(
      gp = grid::gpar(fill = scales::alpha(fill, alpha), col = NA)),
    xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax)
}

drought_band_layers  <- function() {
  list(
    .band_rect(-Inf,  -2,     "#8B0000"),   # Exceptional drought
    .band_rect(-2,    -1.5,   "#FF0000"),   # Extreme drought
    .band_rect(-1.5,  -1,     "#FF8C00"),   # Severe drought
    .band_rect(-1,    -0.5,   "#FFA500"),   # Moderate drought
    .band_rect( 0.5,   1,     "#90EE90"),   # Abnormally moist
    .band_rect( 1,     1.5,   "#00FF00"),   # Moderately moist
    .band_rect( 1.5,   2,     "#008000"),   # Very moist
    .band_rect( 2,     Inf,   "#006400"),   # Exceptionally moist
    ggplot2::geom_hline(yintercept = c(-2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2),
                        linetype = "dotted", color = "gray50", linewidth = 0.3)
  )
}

shared_ts_theme <- function(base_size = 12) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "gray30"),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom"
    )
}

index_colours <- c(spi = "#e41a1c", spei = "#377eb8", swei = "#1E90FF")

calc_dynamic_ylim <- function(v, base = c(-4, 4), pad_frac = 0.10) {
  v <- v[is.finite(v)]
  if (!length(v)) return(base)
  r <- range(v)
  span <- diff(r); if (!is.finite(span) || span == 0) span <- 1
  c(min(base[1], r[1] - pad_frac * span),
    max(base[2], r[2] + pad_frac * span))
}

safe_pdf <- function(filename, width = 12, height = 8) {
  while (grDevices::dev.cur() > 1) try(grDevices::dev.off(), silent = TRUE)
  tryCatch({
    grDevices::pdf(filename, width = width, height = height, family = "Helvetica")
    TRUE
  }, error = function(e) {
    cat("ERROR: Cannot open PDF '", filename, "': ", e$message, "\n")
    FALSE
  })
}

#' Unified figure-saving helper — replaces the scattered mix of pdf()/dev.off()
#' and ggsave() calls across the w* scripts.
#'
#' Two modes:
#'   ggplot2 mode  — pass a ggplot object to \code{plot_obj}; saves both
#'                   <stem>.pdf and <stem>.png automatically.
#'   base-R mode   — leave \code{plot_obj = NULL}; opens a pdf() device so
#'                   the caller can draw with base graphics and close the
#'                   device with dev.off() as usual.  No .png is written in
#'                   this mode (base-R PNG requires a separate png() call).
#'
#' @param plot_obj  ggplot2 object, or NULL for base-R pdf output.
#' @param stem      File path WITHOUT extension; .pdf and .png are appended.
#' @param width     Figure width in inches  (default FIG_WIDTH_WIDE).
#' @param height    Figure height in inches (default FIG_HEIGHT_STD).
#' @param dpi       PNG resolution in ppi   (default FIG_DPI).
#' @return Invisible NULL.
save_figure <- function(plot_obj = NULL,
                        stem,
                        width  = FIG_WIDTH_WIDE,
                        height = FIG_HEIGHT_STD,
                        dpi    = FIG_DPI) {
  pdf_path <- paste0(stem, ".pdf")
  png_path <- paste0(stem, ".png")
  
  if (is.null(plot_obj)) {
    ## Base-R mode: open a pdf device; caller draws, then calls dev.off()
    safe_pdf(pdf_path, width = width, height = height)
    return(invisible(NULL))
  }
  
  ## ggplot2 mode: save both formats in one call
  tryCatch(
    ggplot2::ggsave(pdf_path, plot_obj,
                    width = width, height = height, units = "in", device = "pdf"),
    error = function(e) cat(sprintf("  ⚠ PDF save failed (%s): %s\n",
                                    basename(pdf_path), e$message))
  )
  tryCatch(
    ggplot2::ggsave(png_path, plot_obj,
                    width = width, height = height, units = "in",
                    dpi = dpi, device = "png"),
    error = function(e) cat(sprintf("  ⚠ PNG save failed (%s): %s\n",
                                    basename(png_path), e$message))
  )
  cat(sprintf("  ✓ Figure saved: %s (.pdf + .png)\n", basename(stem)))
  invisible(NULL)
}

## H. EXCEL SUMMARY EXPORT ───────────────────────────────────────────
export_summary_excel <- function(timescales, index_types, ts_loader_fn, output_file) {
  wb <- openxlsx::createWorkbook()
  hdr_style <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#D3D3D3")
  stats_df <- data.frame(Timescale = character(), Index = character(),
                         Mean = numeric(), Median = numeric(),
                         StdDev = numeric(), Min = numeric(), Max = numeric(),
                         Drought_Events_2020_2025 = integer(),
                         stringsAsFactors = FALSE)
  for (idx in index_types) {
    for (sc in timescales) {
      tryCatch({
        df  <- ts_loader_fn(idx, sc)
        df_rec  <- df[df$date >= as.Date("2020-01-01") & df$date <= as.Date("2025-12-31"), ]
        ev_count  <- nrow(detect_drought_events(df_rec))
        stats_df  <- rbind(stats_df, data.frame(
          Timescale = sprintf("%02d", sc),
          Index = toupper(idx),
          Mean = mean(df$value, na.rm = TRUE),
          Median = median(df$value, na.rm = TRUE),
          StdDev = sd(df$value, na.rm = TRUE),
          Min = min(df$value, na.rm = TRUE),
          Max = max(df$value, na.rm = TRUE),
          Drought_Events_2020_2025 = ev_count,
          stringsAsFactors = FALSE))
      }, error = function(e)
        cat(sprintf("  ⚠ Stats for %s-%02d: %s\n", toupper(idx), sc, e$message)))
    }
  }
  openxlsx::addWorksheet(wb, "Summary_Statistics")
  openxlsx::writeData(wb, "Summary_Statistics", stats_df)
  openxlsx::addStyle(wb, "Summary_Statistics", hdr_style,
                     rows = 1, cols = seq_len(ncol(stats_df)))
  openxlsx::saveWorkbook(wb, output_file, overwrite = TRUE)
  cat(sprintf("✓ Excel summary saved: %s\n", basename(output_file)))
  invisible(wb)
}
## I. W9 SHARED-STATE LOADER ─────────────────────────────────────────
#' Load and validate the shared-state RDS written by w9_atmospheric_diagnostics.R.
#'
#' w10c calls readRDS() on this file and unpacks every key into its local
#' environment.  If w9 is ever modified (keys renamed, added, or removed) the
#' downstream scripts previously failed with obscure NULL-reference errors deep
#' inside their analysis code.  This wrapper validates the schema up front and
#' gives a clear, actionable error message instead.
#'
#' Required keys reflect the full set documented in w9 Section 5.  If w9 adds
#' a new key that w10c depends on, add it to \code{required_keys} here.
#'
#' @param path  Full path to w9_shared_state.rds (produced by w9).
#' @return The named list exactly as saved by w9, after schema validation.
read_w9_state <- function(path) {
  if (!file.exists(path))
    stop("w9 shared-state RDS not found:\n  ", path,
         "\n  Run w9_atmospheric_diagnostics.R first.")
  
  st <- readRDS(path)
  
  required_keys <- c(
    ## Data objects
    "base_date_index", "z500_ridge_ts", "slp_nwbc_ts", "sst_nepac_ts",
    "map_layers",
    ## Anomaly NetCDF paths
    "anomaly_nc_z500", "anomaly_nc_slp", "anomaly_nc_sst",
    ## Path scalars
    "WD_PATH", "DATA_DIR", "OUT_DIR", "SPI_SEAS_DIR", "SPEI_SEAS_DIR",
    ## Analysis window
    "START_YEAR", "END_YEAR", "CLIM_START", "CLIM_END",
    "RECENT_START", "RECENT_END",
    "DROUGHT_FOCUS_START", "DROUGHT_FOCUS_END",
    ## Spatial boxes
    "RIDGE_LON_MIN", "RIDGE_LON_MAX", "RIDGE_LAT_MIN", "RIDGE_LAT_MAX",
    "SST_MEAN_LON_MIN", "SST_MEAN_LON_MAX", "SST_MEAN_LAT_MIN", "SST_MEAN_LAT_MAX",
    ## EOF domain (w10c)
    "EOF_LON_MIN", "EOF_LON_MAX", "EOF_LAT_MIN", "EOF_LAT_MAX",
    ## Plot settings
    "FIGURE_DPI", "FIGURE_WIDTH_WIDE", "FIGURE_WIDTH_STD", "MONTH_LABELS"
  )
  
  missing_keys <- setdiff(required_keys, names(st))
  if (length(missing_keys))
    stop("w9 shared-state is missing ", length(missing_keys), " required key(s):\n",
         "  ", paste(missing_keys, collapse = ", "), "\n",
         "  The RDS schema may have changed since w9 last ran.\n",
         "  Re-run w9_atmospheric_diagnostics.R to regenerate the RDS.")
  
  cat(sprintf("✓ w9 shared state loaded (%d keys): %s\n",
              length(st), basename(path)))
  st
}

#source("utils_teleconnection_addon.R")
cat("✓ DROUGHT_ANALYSIS_utils.R loaded\n")