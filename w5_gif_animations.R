####################################################################################
# w5_gif_animations.R  ·  ANIMATED GIFS FOR ALL INDICES
####################################################################################

source("DROUGHT_ANALYSIS_utils.R")
utils_load_packages(c("magick", "scales", "lubridate", "patchwork",
                      "ggplot2", "terra"))

if (!dir.exists(WD_PATH)) stop("Working directory not found: ", WD_PATH)
setwd(WD_PATH)
dir.create(GIF_DIR, showWarnings = FALSE, recursive = TRUE)

GIF_START_YEAR <- 1981
GIF_FPS <- 4
DROUGHT_COLOURS <- c("#8B0000", "#FF0000", "#FFA500", "#FFFF00",
                     "#ADFF2F", "#00FA9A", "#1E90FF")

build_stack <- function(data_dir, index_type, scale, find_fn) {
  files <- find_fn(data_dir, index_type, scale)
  if (!length(files)) return(NULL)
  stacks <- lapply(files, terra::rast)
  dates <- do.call(c, lapply(seq_along(stacks), function(i)
    extract_dates_from_nc(files[i], terra::nlyr(stacks[[i]]))))
  combined <- do.call(c, stacks)
  list(stack = combined, dates = dates)
}

build_frame <- function(raster, basin_shp, label, scale, date_lbl) {
  df <- as.data.frame(raster, xy = TRUE, na.rm = FALSE)
  colnames(df)[1:3] <- c("x", "y", "value")
  
  # Extract basin geometry - ensure unique column names
  basin_geom_raw <- as.data.frame(terra::geom(basin_shp))
  basin_geom <- data.frame(
    x = basin_geom_raw$x,
    y = basin_geom_raw$y,
    group = 1
  )
  
  ggplot2::ggplot() +
    ggplot2::geom_raster(data = df, ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_path(data = basin_geom, ggplot2::aes(x = x, y = y, group = group),
                       color = "darkgreen", linewidth = 1.2) +
    ggplot2::scale_fill_gradientn(colors = DROUGHT_COLOURS,
                                  limits = c(-4, 4), name = label) +
    ggplot2::labs(title = sprintf("%s-%d", label, scale), subtitle = date_lbl) +
    shared_ts_theme(11) +
    ggplot2::coord_equal()
}

render_gif <- function(index_type, scale, data_dir, find_fn, basin_shp) {
  label <- toupper(index_type)
  seasonal <- build_stack(data_dir, index_type, scale, find_fn)
  if (is.null(seasonal)) return()
  
  cat(sprintf("\n%s-%d: %d frames\n", label, scale, terra::nlyr(seasonal$stack)))
  frames <- vector("list", terra::nlyr(seasonal$stack))
  
  for (i in seq_along(frames)) {
    if (i %% 24 == 1) cat(sprintf("  Frame %d\n", i))
    p <- build_frame(seasonal$stack[[i]], basin_shp, label, scale,
                     format(seasonal$dates[i], "%B %Y"))
    tmp <- tempfile(fileext = ".png")
    ggplot2::ggsave(tmp, p, width = 10, height = 8, dpi = 90)
    frames[[i]] <- magick::image_read(tmp)
    file.remove(tmp)
  }
  
  gif <- magick::image_animate(magick::image_join(frames), fps = GIF_FPS)
  out <- file.path(GIF_DIR, sprintf("%s_%d_animation.gif", label, scale))
  magick::image_write(gif, out)
  cat(sprintf("  ✓ %s\n", basename(out)))
}

basin_shp <- terra::vect(BASIN_SHP)

cat("\n── SPI ──\n")
for (sc in SPI_SCALES) render_gif("spi", sc, SPI_SEAS_DIR, find_seasonal_nc_files, basin_shp)

cat("\n── SPEI ──\n")
for (sc in SPEI_SCALES) render_gif("spei", sc, SPEI_SEAS_DIR, find_seasonal_nc_files, basin_shp)

cat("\n── SWEI ──\n")
render_gif("swei", SWEI_SCALE, SWEI_SEAS_DIR, find_swei_seasonal_files, basin_shp)

cat("\n✓ 04_gif_animations.R complete\n")