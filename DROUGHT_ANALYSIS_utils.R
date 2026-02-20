####################################################################################
# DROUGHT_ANALYSIS_utils.R  ·  SHARED UTILITIES FOR SPI, SPEI, AND SWEI
####################################################################################
## A. PROJECT CONFIGURATION ───────────────────────────────────────────
WD_PATH <- Sys.getenv("NECHAKO_WD", "D:/Nechako_Drought/Nechako/")
BASIN_SHP       <- file.path(WD_PATH, "Spatial/nechakoBound_dissolve.shp")

## Index-specific directories
SPI_SEAS_DIR    <- file.path(WD_PATH, "spi_results_seasonal/")
SPEI_SEAS_DIR   <- file.path(WD_PATH, "spei_results_seasonal/")
SWEI_SEAS_DIR   <- file.path(WD_PATH, "swei_results_seasonal/")

## Output directories
TREND_DIR        <- file.path(WD_PATH,  "temporal_drought/")
BASIN_PLOT_DIR   <- file.path(WD_PATH,  "basin_averaged_plots/")
POINT_PLOT_DIR   <- file.path(WD_PATH,  "point_timeseries_plots/")
GIF_DIR <- file.path(TREND_DIR, "drought_gifs/")
CACHE_DIR        <- file.path(WD_PATH,  "cache/")

## CRS and thresholds
EQUAL_AREA_CRS  <- "EPSG:3005"
DROUGHT_ONSET      <- -1.0
DROUGHT_END        <-  0.0
SEVERE_THRESHOLD   <- -1.3
EXTREME_THRESHOLD  <- -1.6
DROUGHT_THRESHOLD  <- -1.0

## Timescales
SPI_SCALES          <- c(1, 3, 6, 12)
SPEI_SCALES         <- c(1, 3, 6, 12)
SWEI_SCALE          <- 3  # Single timescale only
INDICES             <- c("spi", "spei", "swei")
TIMESCALES_STANDARD <- SPI_SCALES          # used by w4 for time series plots
TIMESCALES_SPATIAL  <- SPI_SCALES          # used by w4 for spatial figures
MONTH_NAMES    <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

## B. PACKAGE BOOTSTRAP ──────────────────────────────────────────────
utils_load_packages <- function(pkgs) {
  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing)) {
    message("Installing missing packages: ", paste(missing, collapse = ", "))
    install.packages(missing, repos = "https://cloud.r-project.org/")
  }
  invisible(lapply(pkgs, library, character.only = TRUE, quietly = TRUE))
  cat("✓ Packages loaded:", paste(pkgs, collapse = ", "), "\n")
}

## C. NETCDF FILE FINDERS ────────────────────────────────────────────
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

find_swei_seasonal_files <- function(data_dir, scale) {
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
    cat("  [dates] Using fallback monthly sequence (Jan 1954)\n")
    d  <- seq(as.Date("1954-01-01"), by = "month", length.out = n_layers)
    if (tail(d, 1) > as.Date("2024-12-31")) {
      end  <- as.Date("2024-12-01")
      start  <- end - lubridate::months(n_layers - 1)
      d  <- seq(start, by = "month", length.out = n_layers)
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
  total_area  <- terra::expanse(b_proj, unit = "m")
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
detect_drought_events  <- function(df,
                                   onset_threshold = DROUGHT_ONSET,
                                   termination_threshold = DROUGHT_END,
                                   min_duration = 1,
                                   severe_thr = SEVERE_THRESHOLD,
                                   extreme_thr = EXTREME_THRESHOLD) {
  stopifnot(all(c("date", "value") %in% names(df)))
  df  <- df[order(df$date), ]
  vals  <- df$value
  dates  <- df$date
  events  <- data.frame(event_id = integer(),
                        start_date = as.Date(character()),
                        end_date = as.Date(character()),
                        duration_months = integer(),
                        min_value = numeric(),
                        severity = character(),
                        stringsAsFactors = FALSE)
  in_drought  <- FALSE
  start_idx  <- NA
  classify_severity  <- function(min_v) {
    if (min_v <= extreme_thr) "Extreme"
    else if (min_v <= severe_thr) "Severe"
    else "Moderate"
  }
  add_event  <- function(s, e) {
    dur  <- e - s + 1
    if (dur < min_duration) return()
    mv  <- min(vals[s:e], na.rm = TRUE)
    events  <<- rbind(events, data.frame(
      event_id = nrow(events) + 1,
      start_date = dates[s],
      end_date = dates[e],
      duration_months = dur,
      min_value = mv,
      severity = classify_severity(mv),
      stringsAsFactors = FALSE))
  }
  for (i in seq_along(vals)) {
    if (!in_drought && !is.na(vals[i]) && vals[i] < onset_threshold) {
      in_drought  <- TRUE; start_idx  <- i
    } else if (in_drought && !is.na(vals[i]) && vals[i] >= termination_threshold) {
      add_event(start_idx, i-1); in_drought  <- FALSE
    }
  }
  if (in_drought) add_event(start_idx, length(vals))
  events
}

## F2. SPATIAL HELPERS (used by w4_trends_visualization.R) ──────────
#' Load basin shapefile and reproject to equal-area CRS
load_basin <- function(shp_path = BASIN_SHP, crs = EQUAL_AREA_CRS) {
  if (!file.exists(shp_path)) stop("Basin shapefile not found: ", shp_path)
  b <- sf::st_read(shp_path, quiet = TRUE)
  if (sf::st_crs(b)$input != crs) b <- sf::st_transform(b, crs)
  cat(sprintf("✓ Basin loaded (%d polygon(s), CRS: %s)\n", nrow(b), crs))
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

cat("✓ DROUGHT_ANALYSIS_utils.R loaded\n")