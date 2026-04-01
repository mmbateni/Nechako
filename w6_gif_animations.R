####################################################################################
# w6_gif_animations.R · ANIMATED GIFS FOR ALL INDICES (Optimised & Enhanced)
####################################################################################

setwd("D:/Nechako_Drought/Nechako/")
source("DROUGHT_ANALYSIS_utils.R")
utils_load_packages(c(
  # Core
  "terra", "ggplot2", "lubridate", "scales",
  # PNG renderer — faster than base png(), no gganimate needed
  "ragg",
  # GIF encoder — called once after all PNGs are written
  "gifski",
  # Colour
  "scico",
  # Cartographic elements
  "ggspatial",
  # Progress + logging
  "cli"
))

if (!dir.exists(WD_PATH)) stop("Working directory not found: ", WD_PATH)
setwd(WD_PATH)
dir.create(GIF_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Constants ──────────────────────────────────────────────────────────────────
GIF_START_YEAR <- 1981
GIF_FPS <- 4
DROUGHT_LIMITS <- c(-4, 4)
FRAME_W <- 900L # px
FRAME_H <- 720L # px
FRAME_RES <- 110L # dpi

# ── 1. Build raster stack ─────────────────────────────────────────────────────
build_stack <- function(data_dir, index_type, scale, find_fn) {
  files <- find_fn(data_dir, index_type, scale)
  if (!length(files)) return(NULL)
  stacks <- lapply(files, terra::rast)
  dates <- do.call(c, lapply(seq_along(stacks), function(i)
    extract_dates_from_nc(files[i], terra::nlyr(stacks[[i]]))))
  list(stack = do.call(c, stacks), dates = dates)
}

# (build_frame removed — we will build a single reusable plot and mutate data)

# ── 3. Render all frames then encode GIF ──────────────────────────────────────
render_gif <- function(index_type, scale, data_dir, find_fn, basin_shp,
                       gif_idx = NA, n_gifs = NA, session_start = NA) {
  
  label <- toupper(index_type)
  tag <- if (is.na(scale)) label else sprintf("%s-%d", label, scale)
  title_str <- tag
  
  # ── Load data ──────────────────────────────────────────────────────────────
  seasonal <- build_stack(data_dir, index_type, scale, find_fn)
  if (is.null(seasonal)) {
    cli::cli_alert_warning("{tag}: no files found \u2014 skipping.")
    return(invisible(NULL))
  }
  r <- seasonal$stack
  dates <- seasonal$dates
  n <- terra::nlyr(r)
  
  # Guard: trim if date count mismatches layer count
  if (length(dates) != n) {
    cli::cli_alert_warning(
      "{tag}: date count ({length(dates)}) != nlyr ({n}). Trimming to min."
    )
    n <- min(length(dates), n)
    dates <- dates[seq_len(n)]
  }
  
  cli::cli_h2("{tag} [{n} frames] \u2500 GIF {gif_idx}/{n_gifs}")
  if (!is.na(session_start)) {
    overall_elapsed <- round(proc.time()["elapsed"] - session_start, 1)
    gifs_done <- gif_idx - 1L
    if (gifs_done > 0L) {
      overall_eta <- round((overall_elapsed / gifs_done) * (n_gifs - gifs_done), 0)
      cli::cli_alert_info(
        "Overall: {gifs_done}/{n_gifs} GIFs done \u00b7 {overall_elapsed}s elapsed \u00b7 ~{overall_eta}s remaining"
      )
    } else {
      cli::cli_alert_info("Overall: first GIF \u2014 no overall ETA yet")
    }
  }
  
  # ── OPTIONAL: crop/mask to basin extent to reduce pixels rendered ──────────
  # (Keeps the same information if you only need the basin. Commented to preserve
  # current framing exactly.)
  # r <- terra::mask(terra::crop(r, terra::ext(basin_shp)), basin_shp)
  
  # ── Extract coordinates and values ONCE for the entire stack ───────────────
  cli::cli_alert_info("Extracting raster values\u2026")
  vals <- terra::values(r, mat = TRUE)  # ncell x nlyr
  if (n < ncol(vals)) vals <- vals[, seq_len(n), drop = FALSE]
  xy <- terra::crds(r, df = TRUE, na.rm = FALSE)  # x, y (matches values())
  
  # Basin outline extracted once, reused for every frame
  basin_geom <- data.frame(
    x = terra::geom(basin_shp)[, "x"],
    y = terra::geom(basin_shp)[, "y"],
    group = 1L
  )
  
  # ── Precompute immutable UI bits ────────────────────────────────────────────
  date_lbls <- format(dates, "%B %Y")
  
  # Reusable data frame backing both raster + contour layers
  # We mutate only 'value' per frame
  df0 <- data.frame(x = xy$x, y = xy$y, value = vals[, 1])
  
  # Build the base plot ONCE; identical appearance and scales
  p_base <-
    ggplot2::ggplot() +
    ggplot2::geom_raster(
      data = df0,
      mapping = ggplot2::aes(x = x, y = y, fill = value)
    ) +
    ggplot2::geom_contour(
      data = df0,
      mapping = ggplot2::aes(x = x, y = y, z = value),
      bins = 6,                    # unchanged from original (no method change)
      colour = "grey25",
      linewidth = 0.35,
      na.rm = TRUE
    ) +
    ggplot2::geom_path(
      data = basin_geom,
      mapping = ggplot2::aes(x = x, y = y, group = group),
      colour = "#1B5E20",
      linewidth = 1.2,
      lineend = "round"
    ) +
    scico::scale_fill_scico(
      palette = "vik",
      limits = DROUGHT_LIMITS,
      midpoint = 0,
      name = title_str,
      guide = ggplot2::guide_colorbar(
        barwidth = 0.8,
        barheight = 10,
        title.position = "top"
      )
    ) +
    ggspatial::annotation_scale(
      location = "bl",
      width_hint = 0.3,
      text_cex = 0.75
    ) +
    ggspatial::annotation_north_arrow(
      location = "tr",
      which_north = "true",
      style = ggspatial::north_arrow_fancy_orienteering(
        fill = c("grey80", "grey20"),
        text_size = 8
      ),
      height = ggplot2::unit(1.2, "cm"),
      width  = ggplot2::unit(1.2, "cm")
    ) +
    ggplot2::labs(
      title = title_str,
      subtitle = date_lbls[1],
      caption = "Nechako Watershed \u00b7 Drought Analysis"
    ) +
    shared_ts_theme(11) +
    ggplot2::theme(
      plot.subtitle = ggplot2::element_text(size = 12, face = "bold", colour = "grey30"),
      legend.position = "right",
      panel.grid.major = ggplot2::element_line(colour = "grey92"),
      panel.grid.minor = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "#EEF4F8", colour = NA)
    ) +
    ggplot2::coord_equal(expand = FALSE)
  
  # ── Temp directory for PNGs (auto-deleted on function exit) ────────────────
  tmp_dir <- tempfile("gifframes_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  png_paths <- file.path(tmp_dir, sprintf("frame_%05d.png", seq_len(n)))
  
  # ── Per-frame rendering loop (mutating df0$value only) ─────────────────────
  pb <- cli::cli_progress_bar(
    name = "Rendering frames",
    total = n,
    format = paste0(
      " {cli::pb_spin} {cli::pb_name} ",
      "{cli::pb_bar} ",
      "{cli::pb_current}/{cli::pb_total} ",
      "[{cli::pb_elapsed} elapsed · ETA {cli::pb_eta}]"
    ),
    clear = FALSE
  )
  
  frame_start <- proc.time()["elapsed"]
  for (i in seq_len(n)) {
    # Update only the values vector
    df0$value <- vals[, i]
    
    # Update only the subtitle; everything else is reused
    p <- p_base + ggplot2::labs(subtitle = date_lbls[i])
    
    # ragg::agg_png: fast, high-quality rasterizer (unchanged)
    ragg::agg_png(
      filename = png_paths[i],
      width = FRAME_W,
      height = FRAME_H,
      res = FRAME_RES,
      units = "px"
      # background = "white" # (optional explicit, white by default)
    )
    print(p)
    grDevices::dev.off()
    
    cli::cli_progress_update(id = pb)
    # If console overhead is noticeable, update every 3–5 frames instead:
    # if (i %% 5L == 0L) cli::cli_progress_update(id = pb)
  }
  cli::cli_progress_done(id = pb)
  
  frame_secs <- round(proc.time()["elapsed"] - frame_start, 1)
  cli::cli_alert_info(
    "Frames rendered in {frame_secs}s ({round(frame_secs / n, 2)}s / frame)"
  )
  
  # ── Encode all PNGs into GIF with gifski (unchanged) ───────────────────────
  out <- file.path(
    GIF_DIR,
    if (is.na(scale)) sprintf("%s_animation.gif", label)
    else sprintf("%s_%d_animation.gif", label, scale)
  )
  
  cli::cli_alert_info("Encoding GIF with gifski\u2026")
  gifski::gifski(
    png_files = png_paths,
    gif_file  = out,
    width     = FRAME_W,
    height    = FRAME_H,
    delay     = 1 / GIF_FPS  # seconds per frame
  )
  size_mb <- round(file.size(out) / 1e6, 1)
  
  if (!is.na(session_start)) {
    overall_elapsed <- round(proc.time()["elapsed"] - session_start, 1)
    gifs_left <- n_gifs - gif_idx
    if (gif_idx > 1L) {
      overall_eta <- round((overall_elapsed / gif_idx) * gifs_left, 0)
      eta_str <- if (gifs_left > 0L) sprintf(" \u00b7 ~%ds remaining", overall_eta) else " \u00b7 last GIF"
    } else {
      eta_str <- ""
    }
    cli::cli_alert_success(
      "Saved: {.file {basename(out)}} ({size_mb} MB) \u00b7 overall {overall_elapsed}s elapsed{eta_str}"
    )
  } else {
    cli::cli_alert_success("Saved: {.file {basename(out)}} ({size_mb} MB)")
  }
}

# ── 4. Run for all indices ─────────────────────────────────────────────────────
basin_shp <- terra::vect(BASIN_SHP)
jobs <- c(
  lapply(SPI_SCALES,  function(sc) list("spi",  sc, SPI_SEAS_DIR,  find_seasonal_nc_files)),
  lapply(SPEI_SCALES, function(sc) list("spei", sc, SPEI_SEAS_DIR, find_seasonal_nc_files)),
  list( list("swei", NA_integer_, SWEI_SEAS_DIR, find_swei_seasonal_files) )
)
n_jobs <- length(jobs)

cli::cli_rule(left = "Nechako Drought \u00b7 GIF Animation")
cli::cli_alert_info("{n_jobs} index/scale combinations to render")

# Outer bar — one tick per completed GIF
pb_outer <- cli::cli_progress_bar(
  name  = "Overall progress",
  total = n_jobs,
  format = paste0(
    "{cli::pb_name} {cli::pb_bar} ",
    "{cli::pb_current}/{cli::pb_total} GIFs ",
    "[elapsed {cli::pb_elapsed} · ETA {cli::pb_eta}]"
  ),
  clear = FALSE
)

session_start <- proc.time()["elapsed"]
for (i in seq_along(jobs)) {
  job <- jobs[[i]]
  render_gif(job[[1]], job[[2]], job[[3]], job[[4]], basin_shp,
             gif_idx = i, n_gifs = n_jobs, session_start = session_start)
  cli::cli_progress_update(id = pb_outer)
}
cli::cli_progress_done(id = pb_outer)

total_min <- round((proc.time()["elapsed"] - session_start) / 60, 1)
cli::cli_rule()
cli::cli_alert_success(
  "All {n_jobs} GIFs complete \u00b7 total time: {total_min} min"
)
cli::cli_rule()
