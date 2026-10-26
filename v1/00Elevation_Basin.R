####################################################################################
# NECHAKO RIVER BASIN — COMPREHENSIVE OVERVIEW MAP  (v4, MERGED)
#
# Produces BOTH map variants from a single run (data are prepared only once):
#   • "original" : legend row "Paved road" (swatch lw 0.70), thin separator line
#                  under the elevation swatch, semi-transparent (92%) legend box,
#                  legend inset at left = 0.74 / bottom = 0.683
#   • "modified" : legend row "Main highways" (swatch lw = map road lw), no
#                  separator, fully opaque legend box, legend inset at
#                  left = 0.685 / bottom = 0.665
#
# Outputs (in out_dir):
#   nechako_overview_map_v2_original.pdf / .png
#   nechako_overview_map_v2_modified.pdf / .png
#
# Data sources:
#   • DEM       : AWS Terrain Tiles via elevatr (z = 8, ~300 m)
#   • Rivers    : BC Freshwater Atlas via bcdata  [PRIMARY]
#   • Lakes     : BC Freshwater Atlas via bcdata  [PRIMARY]
#   • Fraser WS : BC Freshwater Atlas watershed groups via bcdata
#   • Roads     : BC Digital Road Atlas via bcdata [PRIMARY]
#   • Sub-basins: user-supplied shapefiles (Nechako_shps)
#   • Fallback  : HydroSHEDS (rivers) + NaturalEarth lakes + NaturalEarth roads
####################################################################################

# ── 0. PACKAGES ─────────────────────────────────────────────────────────────────
required_pkgs <- c(
  "terra", "sf", "ggplot2", "elevatr",
  "bcdata",
  "geodata",
  "ggspatial", "ggrepel", "ggnewscale",
  "rnaturalearth", "rnaturalearthdata",
  "patchwork", "dplyr", "scales"
)
new_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(new_pkgs) > 0) {
  message("Installing missing packages: ", paste(new_pkgs, collapse = ", "))
  install.packages(new_pkgs, quiet = TRUE)
}
invisible(lapply(required_pkgs, library, character.only = TRUE))
if (!requireNamespace("rnaturalearthhires", quietly = TRUE)) {
  install.packages("rnaturalearthhires",
                   repos = "https://ropensci.r-universe.dev", quiet = TRUE)
}
has_bcdata <- requireNamespace("bcdata", quietly = TRUE)

# ── 1. USER PATHS ────────────────────────────────────────────────────────────────
basin_shp    <- "D:/Nechako_Drought/Nechako/Spatial/nechakoBound_dissolve.shp"
out_dir      <- "D:/Nechako_Drought/Nechako/elevation"
subbasin_dir <- "D:/Nechako_Drought/Nechako/Spatial/Nechako_shps"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# DEM cache — saved after first successful download
dem_cache <- file.path(out_dir, "nechako_dem_z8.tif")

# ── 2. BASIN BOUNDARY ────────────────────────────────────────────────────────────
cat("── Loading basin boundary ──\n")
basin_sf   <- st_read(basin_shp, quiet = TRUE) |> st_transform(4326)
basin_v    <- vect(basin_sf)
basin_3005 <- st_transform(basin_sf, 3005)

# Simplified geometry for bcdata queries (< 500 KB threshold)
basin_3005_query <- st_simplify(basin_3005, dTolerance = 1000,
                                preserveTopology = TRUE)
basin_bbox <- st_bbox(basin_sf)
cat(sprintf("  Basin extent: %.3f to %.3f W, %.3f to %.3f N\n",
            -basin_bbox["xmax"], -basin_bbox["xmin"],
            basin_bbox["ymin"],  basin_bbox["ymax"]))

# Fixed map extent
map_xlim <- c(-128.5, -122.0)
map_ylim <- c(  52.9,   56.2)

# Expanded query bbox in BC Albers for road + Fraser watershed queries
query_bbox_3005 <- st_as_sfc(
  st_bbox(c(xmin = map_xlim[1] - 0.5, xmax = map_xlim[2] + 0.5,
            ymin = map_ylim[1] - 0.5, ymax = map_ylim[2] + 0.5),
          crs = 4326)
) |> st_transform(3005)

# ── 3. DEM ───────────────────────────────────────────────────────────────────────
options(timeout = 300)

if (file.exists(dem_cache)) {
  cat("── Loading cached DEM ──\n")
  cat(sprintf("  Cache: %s\n", dem_cache))
  dem <- rast(dem_cache)
} else {
  cat("── Downloading DEM via elevatr (z = 8) ──\n")
  cat("  (result will be cached to avoid re-downloading)\n")
  dir.create(tempdir(), showWarnings = FALSE, recursive = TRUE)

  dem_raw <- tryCatch(
    get_elev_raster(locations = basin_sf, z = 8, clip = "bbox",
                    src = "aws", neg_to_na = FALSE, verbose = FALSE),
    error = function(e) {
      stop(paste0(
        "elevatr download failed: ", conditionMessage(e), "\n",
        "  Possible causes:\n",
        "    1. No internet connection — check network access.\n",
        "    2. Firewall / proxy blocking AWS (s3.amazonaws.com).\n",
        "       Try: options(httr_config = httr::use_proxy('<host>', <port>))\n",
        "       before re-running.\n",
        "    3. Transient AWS outage — wait a few minutes and retry.\n",
        "  Once a successful download completes, it is cached at:\n",
        "    ", dem_cache, "\n",
        "  and the network call will never be needed again."
      ))
    }
  )

  if (is.null(dem_raw)) {
    stop(paste0(
      "get_elev_raster() returned NULL (no error thrown).\n",
      "This usually means there is no internet access or the AWS endpoint\n",
      "is unreachable. See proxy note above.\n",
      "Cache path (once download succeeds): ", dem_cache
    ))
  }

  dem <- rast(dem_raw)
  writeRaster(dem, dem_cache, overwrite = TRUE)
  cat(sprintf("  DEM cached to: %s\n", dem_cache))
}

dem <- mask(dem, basin_v, touches = TRUE)
names(dem) <- "elevation"
elev_range <- as.numeric(global(dem, c("min","max"), na.rm = TRUE))
cat(sprintf("  Elevation range: %.0f – %.0f m\n", elev_range[1], elev_range[2]))

# ── 4. HILLSHADE ─────────────────────────────────────────────────────────────────
cat("── Computing hillshade ──\n")
dem_alb  <- project(dem, "EPSG:3005", method = "bilinear")
slope_r  <- terrain(dem_alb, "slope",  unit = "radians")
aspect_r <- terrain(dem_alb, "aspect", unit = "radians")
hill     <- shade(slope_r, aspect_r, angle = 40, direction = 315, normalize = TRUE)
hill     <- project(hill, "EPSG:4326", method = "bilinear")
hill     <- resample(hill, dem, method = "bilinear")
names(hill) <- "hillshade"

# ── 5. RIVERS & LAKES ────────────────────────────────────────────────────────────
cat("── Fetching rivers and lakes ──\n")
rivers_sf      <- NULL
rivers_poly_sf <- NULL
lakes_sf       <- NULL

## ── 5a. PRIMARY: BC Freshwater Atlas ─────────────────────────────────────────
if (has_bcdata) {
  tryCatch({
    cat("  Source: BC Freshwater Atlas (bcdata)\n")
    rivers_sf <- bcdc_query_geodata("freshwater-atlas-stream-network") |>
      filter(STREAM_ORDER >= 5) |>
      filter(INTERSECTS(basin_3005_query)) |>
      select(GNIS_NAME, STREAM_ORDER) |>
      collect() |> st_transform(4326) |>
      st_intersection(st_union(basin_sf))

    rivers_poly_sf <- tryCatch({
      bcdc_query_geodata("freshwater-atlas-rivers") |>
        filter(INTERSECTS(basin_3005_query)) |>
        select(GNIS_NAME) |>
        collect() |> st_transform(4326) |>
        st_intersection(st_union(basin_sf))
    }, error = function(e) { message("  rivers poly: ", e$message); NULL })

    lakes_sf <- bcdc_query_geodata("freshwater-atlas-lakes") |>
      filter(AREA_HA >= 5) |>
      filter(INTERSECTS(basin_3005_query)) |>
      select(GNIS_NAME_1, AREA_HA) |>
      collect() |> st_transform(4326) |>
      st_intersection(st_union(basin_sf))

    cat(sprintf("  Rivers: %d | River polys: %d | Lakes: %d\n",
                nrow(rivers_sf),
                if (!is.null(rivers_poly_sf)) nrow(rivers_poly_sf) else 0L,
                nrow(lakes_sf)))
  }, error = function(e) {
    message("  bcdata failed: ", conditionMessage(e), "\n  Falling back ...")
    rivers_sf <<- NULL; rivers_poly_sf <<- NULL; lakes_sf <<- NULL
  })
}

## ── 5b. FALLBACK ──────────────────────────────────────────────────────────────
if (is.null(rivers_sf)) {
  cat("  Source: geodata + NaturalEarth\n")
  tmp_dir <- file.path(tempdir(), "geodata_cache")
  if (!dir.exists(tmp_dir)) dir.create(tmp_dir, recursive = TRUE)
  tryCatch({
    hydro     <- geodata::rivers(path = tmp_dir)
    rivers_sf <- st_as_sf(hydro) |> st_transform(4326) |> st_make_valid() |>
      st_intersection(st_union(basin_sf))
    names(rivers_sf) <- tolower(names(rivers_sf))
    rivers_sf <- rivers_sf |>
      mutate(STREAM_ORDER = if ("ord_stra" %in% names(rivers_sf)) ord_stra
             else NA_integer_,
             GNIS_NAME    = if ("river"    %in% names(rivers_sf)) river
             else NA_character_)
  }, error = function(e) {
    message("  Rivers failed: ", e$message)
    rivers_sf <<- sf::st_sf(STREAM_ORDER = integer(0), GNIS_NAME = character(0),
                            geometry = sf::st_sfc(crs = 4326))
  })
  tryCatch({
    ne_l <- ne_download(scale = 10, type = "lakes", category = "physical",
                        returnclass = "sf") |>
      st_transform(4326) |> st_make_valid() |>
      st_intersection(st_union(basin_sf))
    lakes_sf <- ne_l |>
      mutate(GNIS_NAME_1 = if ("name" %in% names(ne_l)) name else NA_character_,
             AREA_HA     = if ("area" %in% names(ne_l)) area else NA_real_)
  }, error = function(e) {
    message("  Lakes failed: ", e$message)
    lakes_sf <<- sf::st_sf(GNIS_NAME_1 = character(0), AREA_HA = numeric(0),
                           geometry = sf::st_sfc(crs = 4326))
  })
}

## ── 5c. Stream-order line widths ─────────────────────────────────────────────
if (nrow(rivers_sf) == 0) {
  rivers_sf$line_width <- numeric(0)
} else if ("STREAM_ORDER" %in% names(rivers_sf) &&
           !all(is.na(rivers_sf$STREAM_ORDER))) {
  rivers_sf <- rivers_sf |>
    mutate(line_width = case_when(
      STREAM_ORDER >= 8 ~ 1.10, STREAM_ORDER >= 7 ~ 0.75,
      STREAM_ORDER >= 6 ~ 0.50, STREAM_ORDER >= 5 ~ 0.35, TRUE ~ 0.25))
} else {
  rivers_sf$line_width <- 0.5
}

# ── 5.5 FRASER WATERSHED BOUNDARY ────────────────────────────────────────────────
cat("── Fetching Fraser watershed boundary ──\n")
fraser_sf <- NULL
if (has_bcdata) {
  fraser_sf <- tryCatch({
    cat("  Trying freshwater-atlas-major-watersheds ...\n")
    bcdc_query_geodata("freshwater-atlas-major-watersheds") |>
      filter(MAJOR_WATERSHED_SYSTEM == "FRASER RIVER") |>
      collect() |> st_transform(4326) |> st_union()
  }, error = function(e1) {
    tryCatch({
      cat("  Trying watershed-groups union ...\n")
      bcdc_query_geodata("freshwater-atlas-watershed-groups") |>
        filter(WATERSHED_GROUP_CODE %in% c(
          "BRID", "CANA", "CHIL", "CROL", "DEAD", "ENCL", "FRAN",
          "GREE", "HARR", "HORS", "ILEC", "LCHL", "LILL", "LNIC",
          "MAST", "NAZO", "NECH", "NICL", "SALA", "SALM",
          "SARA", "SHUA", "STUL", "THOM", "WIDG"
          # "QUES" removed: Quesnel watershed rendered as a confusing dark polygon
        )) |>
        collect() |> st_transform(4326) |> st_union()
    }, error = function(e2) {
      message("  Fraser watershed failed: ", e2$message)
      NULL
    })
  })
}
if (!is.null(fraser_sf))
  cat(sprintf("  Fraser watershed: %.0f km²\n",
              as.numeric(st_area(fraser_sf)) / 1e6))

# ── 5.6 ROADS ────────────────────────────────────────────────────────────────────
cat("── Fetching paved highways / arterials ──\n")
roads_sf <- NULL
if (has_bcdata) {
  roads_sf <- tryCatch({
    cat("  Source: BC Digital Road Atlas (bcdata)\n")
    bcdc_query_geodata("digital-road-atlas") |>
      filter(ROAD_SURFACE == "paved") |>
      filter(ROAD_CLASS   %in% c("highway", "arterial")) |>
      filter(INTERSECTS(query_bbox_3005)) |>
      select(ROAD_CLASS) |>
      collect() |> st_transform(4326)
  }, error = function(e) {
    message("  bcdata roads failed: ", e$message, "\n  Trying NaturalEarth ...")
    NULL
  })
}
if (is.null(roads_sf)) {
  roads_sf <- tryCatch({
    ne_download(scale = 10, type = "roads",
                category = "cultural", returnclass = "sf") |>
      st_transform(4326) |> st_make_valid() |>
      filter(!grepl("Ferry", type, ignore.case = TRUE)) |>
      st_crop(st_bbox(c(xmin = map_xlim[1] - 0.3, xmax = map_xlim[2] + 0.3,
                        ymin = map_ylim[1] - 0.3, ymax = map_ylim[2] + 0.3),
                      crs = 4326))
  }, error = function(e) {
    message("  NaturalEarth roads failed: ", e$message); NULL
  })
}
if (!is.null(roads_sf))
  cat(sprintf("  Roads: %d segments\n", nrow(roads_sf)))

# ── 5.7 FRASER RIVER LINE ────────────────────────────────────────────────────────
cat("── Fetching Fraser River line ──\n")
fraser_river_sf <- NULL

fraser_query_bbox <- st_as_sfc(
  st_bbox(c(xmin = map_xlim[1] - 0.3,  xmax = -119.0,
            ymin = map_ylim[1] - 0.3,  ymax = map_ylim[2] + 0.3),
          crs = 4326)
) |> st_transform(3005)

if (has_bcdata) {
  fraser_river_sf <- tryCatch({
    cat("  Source: BC FWA stream network (Fraser River main stem)\n")
    bcdc_query_geodata("freshwater-atlas-stream-network") |>
      filter(GNIS_NAME == "Fraser River") |>
      filter(STREAM_ORDER >= 8) |>
      filter(INTERSECTS(fraser_query_bbox)) |>
      select(GNIS_NAME, STREAM_ORDER) |>
      collect() |>
      st_transform(4326) |>
      st_crop(st_bbox(c(xmin = map_xlim[1], xmax = map_xlim[2],
                        ymin = map_ylim[1], ymax = map_ylim[2]), crs = 4326))
  }, error = function(e) {
    message("  BC FWA Fraser River failed: ", e$message)
    NULL
  })
}

if (is.null(fraser_river_sf)) {
  fraser_river_sf <- tryCatch({
    cat("  Fallback: NaturalEarth rivers_lake_centerlines\n")
    ne_download(scale = 10, type = "rivers_lake_centerlines",
                category = "physical", returnclass = "sf") |>
      st_transform(4326) |>
      st_make_valid() |>
      filter(grepl("Fraser", name, ignore.case = TRUE)) |>
      st_crop(st_bbox(c(xmin = map_xlim[1], xmax = map_xlim[2],
                        ymin = map_ylim[1], ymax = map_ylim[2]), crs = 4326))
  }, error = function(e) {
    message("  NaturalEarth Fraser River failed: ", e$message)
    NULL
  })
}

if (!is.null(fraser_river_sf))
  cat(sprintf("  Fraser River: %d segments\n", nrow(fraser_river_sf)))

# Smooth the Fraser River geometry (union → simplify in metres → back to WGS84)
if (!is.null(fraser_river_sf) && nrow(fraser_river_sf) > 0) {
  fraser_river_sf <- fraser_river_sf |>
    st_transform(3005) |>
    st_union() |>
    st_simplify(dTolerance = 150, preserveTopology = TRUE) |>
    st_transform(4326) |>
    st_as_sf() |>
    rename(geometry = x)
  cat(sprintf("  Fraser River after smoothing: %d feature(s)\n",
              nrow(fraser_river_sf)))
}

# ── 6. SUB-BASINS (user-supplied shapefiles) ─────────────────────────────────────
cat("── Loading Nechako sub-basins ──\n")

subbasin_files <- list.files(subbasin_dir, pattern = "\\.shp$", full.names = TRUE)
if (length(subbasin_files) == 0) {
  warning(sprintf(
    "No .shp files found in subbasin_dir (%s) — sub-basin overlay will be skipped.",
    subbasin_dir))
}

subbasins_sf <- NULL
if (length(subbasin_files) > 0) {
  subbasin_list <- lapply(subbasin_files, function(f) {
    nm <- tools::file_path_sans_ext(basename(f))
    tryCatch({
      s <- st_read(f, quiet = TRUE) |> st_make_valid() |> st_transform(4326)
      s_union <- st_union(st_geometry(s))   # dissolve multi-part sub-basins
      st_sf(name = gsub("_", " ", nm), geometry = s_union)
    }, error = function(e) {
      message("  Failed to read ", basename(f), ": ", conditionMessage(e))
      NULL
    })
  })
  subbasin_list <- subbasin_list[!sapply(subbasin_list, is.null)]

  if (length(subbasin_list) > 0) {
    subbasins_sf <- do.call(rbind, subbasin_list) |>
      st_make_valid() |>
      st_intersection(st_union(basin_sf)) |>
      st_make_valid()
    cat(sprintf("  Sub-basins loaded: %d (%s)\n",
                nrow(subbasins_sf), paste(subbasins_sf$name, collapse = ", ")))
  }
}

# Label anchor points guaranteed inside each polygon
if (!is.null(subbasins_sf) && nrow(subbasins_sf) > 0) {
  subbasin_label_pts <- st_point_on_surface(subbasins_sf)
  subbasin_labels_df <- data.frame(
    name = subbasins_sf$name,
    lon  = st_coordinates(subbasin_label_pts)[, 1],
    lat  = st_coordinates(subbasin_label_pts)[, 2],
    stringsAsFactors = FALSE
  )
} else {
  subbasin_labels_df <- NULL
}

# Endako pinned with geom_text; all others use ggrepel
if (!is.null(subbasin_labels_df)) {
  endako_label_df <- subbasin_labels_df[subbasin_labels_df$name == "Endako", ]
  other_labels_df <- subbasin_labels_df[subbasin_labels_df$name != "Endako", ]
} else {
  endako_label_df <- NULL
  other_labels_df <- NULL
}

# ── 7. MANUAL RIVER LABELS ───────────────────────────────────────────────────────
manual_labels_df <- data.frame(
  lon   = c(-124.35, -122.82),
  lat   = c(  53.92,   53.75),
  label = c("Nechako River", "Fraser River"),
  angle = c(  0, 60),
  stringsAsFactors = FALSE
)

# ── 8. COLOUR PALETTES ──────────────────────────────────────────────────────────
n_col <- 512
elev_pal <- c(
  colorRampPalette(c("#2d6a2d","#4a9b4a","#8ab87a"))(round(n_col * 0.18)),
  colorRampPalette(c("#c8d878","#d4b050","#b88030"))(round(n_col * 0.27)),
  colorRampPalette(c("#b88030","#c8a060","#d8b888"))(round(n_col * 0.25)),
  colorRampPalette(c("#d8b888","#ece4c8","#f8f2e0","#ffffff"))(round(n_col * 0.30))
)

water_fill  <- "#7ec8e3"
water_line  <- "#4a90d9"
river_col   <- "#4a90d9"
bg_land     <- "#eeeeee"
fraser_fill <- "#d0d0d0"
fraser_line <- "#aaaaaa"
basin_col   <- "grey10"
road_col    <- "#888888"
road_lw     <- 0.25      # road line width on the map (both variants)

# ── 9. RASTER DATA FRAMES ──────────────────────────────────────────────────────
cat("── Preparing raster data frames ──\n")
dem_df  <- as.data.frame(dem, xy = TRUE, na.rm = TRUE)
hill_df <- as.data.frame(hill, xy = TRUE, na.rm = FALSE) |>
  mutate(shadow = ifelse(is.na(hillshade), NA_real_,
                         pmax(0, pmin(1, (1 - hillshade) * 0.58))))

# ── 10. VARIANT SETTINGS ─────────────────────────────────────────────────────────
# Everything that differs between the two original scripts lives here.
variants <- list(
  original = list(
    suffix         = "original",
    road_label     = "Paved road",
    road_swatch_lw = 0.70,
    show_separator = TRUE,
    legend_bg      = alpha("white", 0.92),
    legend_pos     = c(left = 0.74,  bottom = 0.683, right = 0.995, top = 0.995)
  ),
  modified = list(
    suffix         = "modified",
    road_label     = "Main highways",
    road_swatch_lw = road_lw,          # swatch matches map line width
    show_separator = FALSE,
    legend_bg      = "white",          # opaque: stops roads bleeding through
    legend_pos     = c(left = 0.685, bottom = 0.665, right = 0.995, top = 0.995)
  )
)

# ── 11. LEGEND BUILDER ───────────────────────────────────────────────────────────
leg_ymin     <- 0.25
leg_ymax     <- 3.25
swatch_xctr  <- 0.32
swatch_width <- 0.52
n_tiles      <- 60
tile_h       <- (leg_ymax - leg_ymin) / n_tiles

elev_grad_df <- data.frame(
  x    = swatch_xctr,
  y    = seq(leg_ymin + tile_h / 2, leg_ymax - tile_h / 2, length.out = n_tiles),
  fill = seq(0, 1, length.out = n_tiles)
)

make_legend <- function(road_label, road_swatch_lw, show_separator, legend_bg) {
  ggplot() +
    # Elevation swatch
    annotate("text", x = 1.02, y = 3.62, hjust = 0.5, vjust = 0, size = 6.00,
             fontface = "bold", color = "grey10", label = "Elevation (m)") +
    geom_tile(data = elev_grad_df, aes(x = x, y = y, fill = fill),
              width = swatch_width, height = tile_h) +
    scale_fill_gradientn(colours = elev_pal, guide = "none") +
    annotate("rect",
             xmin = swatch_xctr - swatch_width / 2,
             xmax = swatch_xctr + swatch_width / 2,
             ymin = leg_ymin, ymax = leg_ymax,
             fill = NA, color = "black", linewidth = 0.4) +
    annotate("text", x = swatch_xctr + swatch_width / 2 + 0.05, y = leg_ymax,
             label = paste0(round(elev_range[2]), " m"),
             hjust = 0, vjust = 1.0, size = 4.84, color = "grey15") +
    annotate("text", x = swatch_xctr + swatch_width / 2 + 0.05, y = leg_ymin,
             label = paste0(round(elev_range[1]), " m"),
             hjust = 0, vjust = 0.0, size = 4.84, color = "grey15") +

    # Optional separator (original variant only)
    { if (show_separator)
        annotate("segment", x = 0.05, xend = 2.05, y = 3.90, yend = 3.90,
                 color = "grey65", linewidth = 0.25) } +

    # Roads row (y = 4.50)
    annotate("segment", x = 0.10, xend = 0.52, y = 4.50, yend = 4.50,
             color = road_col, linewidth = road_swatch_lw) +
    annotate("text", x = 0.65, y = 4.50, hjust = 0, size = 5.88, color = "grey10",
             label = road_label) +

    # Nechako Watershed (y = 5.30)
    annotate("rect", xmin = 0.10, xmax = 0.52, ymin = 4.92, ymax = 5.68,
             fill = NA, color = basin_col, linewidth = 0.85) +
    annotate("text", x = 0.65, y = 5.30, hjust = 0, size = 5.88, color = "grey10",
             label = "Nechako Watershed") +

    # Sub-basin boundary (y = 6.10)
    annotate("segment", x = 0.10, xend = 0.52, y = 6.10, yend = 6.10,
             color = "black", linewidth = 0.90, linetype = "dashed") +
    annotate("text", x = 0.65, y = 6.10, hjust = 0, size = 5.88, color = "grey10",
             label = "Sub-basin boundary") +

    xlim(0.05, 2.90) +
    ylim(0.05, 6.52) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = legend_bg,
                                     color = "grey35", linewidth = 0.40),
      plot.margin     = margin(4, 5, 4, 5)
    )
}

# ── 12. INSET LOCATOR MAP (shared) ───────────────────────────────────────────────
cat("── Building inset locator map ──\n")

bc_sf_inset <- ne_states(country = "Canada", returnclass = "sf") |>
  filter(name == "British Columbia") |>
  st_transform(4326) |>
  st_cast("POLYGON")
bc_sf_inset <- bc_sf_inset[which.max(st_area(bc_sf_inset)), ]

canada_sf <- ne_countries(country = "Canada", scale = "small",
                          returnclass = "sf") |> st_transform(4326)
usa_sf    <- ne_countries(country = "United States of America", scale = "small",
                          returnclass = "sf") |> st_transform(4326)

basin_sf_inset <- basin_sf |>
  st_transform(3005) |>
  st_simplify(dTolerance = 2000, preserveTopology = TRUE) |>
  st_transform(4326)

inset <- ggplot() +
  geom_sf(data = usa_sf,         fill = "gray92", color = "gray70", linewidth = 0.20) +
  geom_sf(data = canada_sf,      fill = "gray88", color = "gray65", linewidth = 0.20) +
  geom_sf(data = bc_sf_inset,    fill = "gray75", color = "gray45", linewidth = 0.35) +
  geom_sf(data = basin_sf_inset, fill = NA, color = "#c0392b",
          linewidth = 0.85, alpha = 0.9) +
  coord_sf(xlim = c(-140, -113), ylim = c(47, 61), expand = FALSE) +
  theme_void(base_size = 7) +
  theme(
    panel.background = element_rect(fill = "#c8dff0", color = NA),
    panel.border     = element_rect(fill = NA, color = "grey25", linewidth = 0.6),
    plot.background  = element_rect(fill = "white", color = NA)
  )

# ── 13. MAIN MAP (shared) ────────────────────────────────────────────────────────
cat("── Rendering main map ──\n")

main_map <- ggplot() +
  theme(panel.background = element_rect(fill = bg_land, color = NA)) +

  {if (!is.null(fraser_sf))
    geom_sf(data = fraser_sf, fill = fraser_fill, color = fraser_line,
            linewidth = 0.3, inherit.aes = FALSE)} +

  geom_raster(data = dem_df, aes(x = x, y = y, fill = elevation),
              interpolate = TRUE) +
  scale_fill_gradientn(
    colours = elev_pal, limits = elev_range,
    breaks  = pretty(elev_range, n = 6),
    name    = "Elevation (m)", na.value = NA,
    guide   = guide_colorbar(
      barwidth = unit(10, "cm"), barheight = unit(0.4, "cm"),
      title.position = "top", title.hjust = 0.5,
      ticks.colour = "black", frame.colour = "black",
      frame.linewidth = 0.4, label.position = "bottom"
    )
  ) +

  ggnewscale::new_scale_fill() +
  geom_raster(data = filter(hill_df, !is.na(shadow)),
              aes(x = x, y = y, alpha = shadow),
              fill = "black", interpolate = TRUE) +
  scale_alpha_identity(guide = "none") +

  {if (!is.null(rivers_poly_sf) && nrow(rivers_poly_sf) > 0)
    geom_sf(data = rivers_poly_sf, fill = water_fill, color = water_line,
            linewidth = 0.10, inherit.aes = FALSE)} +

  geom_sf(data = lakes_sf, fill = water_fill, color = water_line,
          linewidth = 0.10, inherit.aes = FALSE) +

  geom_sf(data = rivers_sf, aes(linewidth = line_width),
          color = river_col, inherit.aes = FALSE) +
  scale_linewidth_identity(guide = "none") +

  {if (!is.null(roads_sf) && nrow(roads_sf) > 0)
    geom_sf(data = roads_sf, color = road_col, linewidth = road_lw,
            inherit.aes = FALSE)} +

  {if (!is.null(fraser_river_sf) && nrow(fraser_river_sf) > 0)
    geom_sf(data = fraser_river_sf, color = water_line, linewidth = 0.65,
            lineend = "round", linejoin = "round", inherit.aes = FALSE)} +

  {if (!is.null(subbasins_sf) && nrow(subbasins_sf) > 0)
    geom_sf(data = subbasins_sf, color = "black", fill = NA,
            linewidth = 0.75, linetype = "dashed", inherit.aes = FALSE)} +

  geom_sf(data = basin_sf, fill = NA, color = basin_col,
          linewidth = 1.2, inherit.aes = FALSE) +

  geom_text(data = manual_labels_df,
            aes(x = lon, y = lat, label = label, angle = angle),
            fontface = "italic", colour = "#1a6699",
            size = 3.56, lineheight = 0.85, inherit.aes = FALSE) +

  {if (!is.null(other_labels_df) && nrow(other_labels_df) > 0)
    ggrepel::geom_text_repel(
      data = other_labels_df, aes(x = lon, y = lat, label = name),
      color = "black", fontface = "bold", size = 4.6,
      nudge_x = 0, nudge_y = 0,
      box.padding = 0.15, point.padding = 0.05,
      min.segment.length = 0.3, max.overlaps = 30,
      force = 0.3, force_pull = 2, seed = 11,
      show.legend = FALSE, inherit.aes = FALSE
    )} +
  {if (!is.null(endako_label_df) && nrow(endako_label_df) > 0)
    geom_text(data = endako_label_df, aes(x = lon, y = lat, label = name),
              color = "black", fontface = "bold", size = 4.6,
              inherit.aes = FALSE)} +

  coord_sf(crs = 4326, expand = FALSE, xlim = map_xlim, ylim = map_ylim) +

  annotation_scale(location = "bl", width_hint = 0.18, style = "bar",
                   text_cex = 0.70,
                   pad_x = unit(0.35, "cm"), pad_y = unit(0.35, "cm")) +

  labs(
    title   = "Nechako Watershed",
    x       = "Longitude (\u00b0W)",
    y       = "Latitude (\u00b0N)",
    caption = paste0(
      "Elevation: AWS Terrain Tiles via elevatr (z=8, ~300 m)  |  ",
      "Rivers & Lakes: BC Freshwater Atlas via bcdata\n",
      "Fraser Watershed: BC FWA watershed groups  |  ",
      "Roads: BC Digital Road Atlas / NaturalEarth  |  ",
      "Basin: Nechako watershed (WGS84)  |  ",
      "Sub-basins: user-supplied shapefiles (Nechako_shps)"
    )
  ) +

  theme_minimal(base_size = 12, base_family = "sans") +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = bg_land, color = NA),
    panel.grid.major = element_line(color = alpha("white", 0.40),
                                    linetype = "dashed", linewidth = 0.20),
    panel.border     = element_rect(fill = NA, color = "grey20", linewidth = 0.7),
    plot.title       = element_text(face = "bold", size = 19, hjust = 0.5,
                                    margin = margin(b = 8)),
    plot.caption     = element_text(size = 7.0, color = "grey50", hjust = 0,
                                    margin = margin(t = 6)),
    axis.title       = element_text(size = 15.6, color = "grey30"),
    axis.text        = element_text(size = 14.4, color = "grey30"),
    legend.position  = "none"
  )

# ── 14. ASSEMBLE + SAVE BOTH VARIANTS ────────────────────────────────────────────
saved_files <- character(0)

for (v in variants) {
  cat(sprintf("── Assembling & saving variant: %s ──\n", v$suffix))

  legend_plot <- make_legend(
    road_label     = v$road_label,
    road_swatch_lw = v$road_swatch_lw,
    show_separator = v$show_separator,
    legend_bg      = v$legend_bg
  )

  lp <- v$legend_pos
  final_map <- main_map +
    inset_element(legend_plot,
                  left = lp[["left"]], bottom = lp[["bottom"]],
                  right = lp[["right"]], top = lp[["top"]],
                  align_to = "panel") +
    inset_element(inset,
                  left = 0.00, bottom = 0.235, right = 0.228, top = 0.736,
                  align_to = "panel")

  pdf_path <- file.path(out_dir, sprintf("nechako_overview_map_v2_%s.pdf", v$suffix))
  ggsave(pdf_path, plot = final_map, width = 14, height = 11, units = "in",
         device = cairo_pdf)
  cat(sprintf("  PDF: %s\n", pdf_path))

  png_path <- file.path(out_dir, sprintf("nechako_overview_map_v2_%s.png", v$suffix))
  ggsave(png_path, plot = final_map, width = 14, height = 11, units = "in",
         dpi = 300, device = "png", bg = "white")
  cat(sprintf("  PNG: %s\n", png_path))

  saved_files <- c(saved_files, pdf_path, png_path)
}

cat("\n══ Done! ══\n")
cat(paste0("  ", basename(saved_files), collapse = "\n"), "\n")