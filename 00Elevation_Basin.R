####################################################################################
# NECHAKO RIVER BASIN — COMPREHENSIVE OVERVIEW MAP
# Elevation · Hillshade · Rivers · Lakes · Prince George · Kenney Dam
#
# Data sources (automatic, no manual downloads needed):
#   • DEM       : AWS Terrain Tiles via elevatr (z = 8, ~300 m)
#   • Rivers    : BC Freshwater Atlas via bcdata  [PRIMARY]
#   • Lakes     : BC Freshwater Atlas via bcdata  [PRIMARY]
#   • Fallback  : HydroSHEDS (rivers) + NHD-derived lakes via geodata package
#
# Packages required:
#   terra, sf, tidyterra, ggplot2, elevatr, bcdata, geodata,
#   ggspatial, ggrepel, ggnewscale, rnaturalearth, rnaturalearthdata,
#   patchwork, dplyr, scales
####################################################################################

# ── 0. PACKAGES ─────────────────────────────────────────────────────────────────
required_pkgs <- c(
  "terra", "sf", "tidyterra", "ggplot2", "elevatr",
  "bcdata",                            # BC Freshwater Atlas — now in auto-install
  "geodata",                           # HydroSHEDS fallback rivers/lakes
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
  install.packages(
    "rnaturalearthhires",
    repos = "https://ropensci.r-universe.dev",  # CRAN doesn't host it
    quiet = TRUE
  )
}
# NOTE: osmdata is no longer used as a fallback — it hangs on large bounding boxes.
# geodata::hydrosheds() and geodata::waterbody() are used instead (fast, tiled).
has_bcdata <- requireNamespace("bcdata", quietly = TRUE)

# ── 1. USER PATHS ────────────────────────────────────────────────────────────────
basin_shp <- "D:/Nechako_Drought/Nechako/Spatial/nechakoBound_dissolve.shp"
out_dir   <- "D:/Nechako_Drought/Nechako/elevation"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ── 2. BASIN BOUNDARY ────────────────────────────────────────────────────────────
cat("── Loading basin boundary ──\n")
basin_sf   <- st_read(basin_shp, quiet = TRUE) |> st_transform(4326)
basin_v    <- vect(basin_sf)
basin_3005 <- st_transform(basin_sf, 3005)

# Simplified geometry for bcdata queries (keeps object < 500 KB threshold)
# dTolerance in metres; 1000 m is sufficient for a river/lake query filter
basin_3005_query <- st_simplify(basin_3005, dTolerance = 1000, preserveTopology = TRUE)
cat(sprintf("  Basin geometry size (query): %.0f bytes\n",
            as.numeric(object.size(basin_3005_query))))

basin_bbox <- st_bbox(basin_sf)
cat(sprintf("  Basin extent: %.3f – %.3f W, %.3f – %.3f N\n",
            -basin_bbox["xmax"], -basin_bbox["xmin"],
            basin_bbox["ymin"],  basin_bbox["ymax"]))

# ── 3. DEM — download via elevatr (AWS Terrain Tiles) ───────────────────────────
cat("── Downloading DEM via elevatr (z = 8, ~300 m resolution) ──\n")
# ── Fix: ensure R's session temp directory exists (Windows elevatr bug) ──
dir.create(tempdir(), showWarnings = FALSE, recursive = TRUE)
options(timeout = 300)  # generous timeout for tile downloads
dem_raw <- get_elev_raster(
  locations = basin_sf,
  z         = 8,        # zoom 8 ≈ 300 m; increase to 9 for ~150 m (slower)
  clip      = "bbox",
  src       = "aws",
  neg_to_na = FALSE,
  verbose   = FALSE
)

dem <- rast(dem_raw)
dem <- mask(dem, basin_v, touches = TRUE)
names(dem) <- "elevation"

elev_range <- as.numeric(global(dem, c("min","max"), na.rm = TRUE))
cat(sprintf("  Elevation range: %.0f – %.0f m\n", elev_range[1], elev_range[2]))

# ── 4. HILLSHADE (computed in BC Albers for correct distances) ──────────────────
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

## ── 5a. PRIMARY: BC Freshwater Atlas via bcdata ───────────────────────────────
if (has_bcdata) {
  tryCatch({
    cat("  Source: BC Freshwater Atlas (bcdata)\n")
    
    rivers_sf <- bcdc_query_geodata("freshwater-atlas-stream-network") |>
      filter(STREAM_ORDER >= 6) |>   # order 6+ for cleaner network
      filter(INTERSECTS(basin_3005_query)) |>
      select(GNIS_NAME, STREAM_ORDER) |>
      collect() |>
      st_transform(4326) |>
      st_intersection(st_union(basin_sf))
    
    # Double-banked river polygons (shows Nechako R. as filled blue channel)
    rivers_poly_sf <- tryCatch({
      bcdc_query_geodata("freshwater-atlas-rivers") |>
        filter(INTERSECTS(basin_3005_query)) |>
        select(GNIS_NAME) |>
        collect() |>
        st_transform(4326) |>
        st_intersection(st_union(basin_sf))
    }, error = function(e) { message("  rivers poly failed: ", e$message); NULL })
    
    lakes_sf <- bcdc_query_geodata("freshwater-atlas-lakes") |>
      filter(AREA_HA >= 5) |>   # lowered from 50 to include reservoir arms
      filter(INTERSECTS(basin_3005_query)) |>
      select(GNIS_NAME_1, AREA_HA) |>
      collect() |>
      st_transform(4326) |>
      st_intersection(st_union(basin_sf))
    
    cat(sprintf("  Rivers: %d segments | River polys: %d | Lakes: %d polygons\n",
                nrow(rivers_sf),
                if (!is.null(rivers_poly_sf)) nrow(rivers_poly_sf) else 0L,
                nrow(lakes_sf)))
    
  }, error = function(e) {
    message("  bcdata query failed: ", conditionMessage(e),
            "\n  Falling back to HydroSHEDS / geodata ...")
    rivers_sf      <<- NULL
    rivers_poly_sf <<- NULL
    lakes_sf       <<- NULL
  })
}

## ── 5b. FALLBACK: rivers via geodata + NaturalEarth lakes ─────────────────────
if (is.null(rivers_sf)) {
  cat("  Source: geodata rivers + NaturalEarth lakes\n")
  
  tmp_dir <- file.path(tempdir(), "geodata_cache")
  if (!dir.exists(tmp_dir)) dir.create(tmp_dir, recursive = TRUE)
  
  # Bug 1 fix: geodata::hydrosheds() doesn't exist — correct function is geodata::rivers()
  tryCatch({
    hydro <- geodata::rivers(path = tmp_dir)   # downloads major rivers globally
    rivers_sf <- st_as_sf(hydro) |>
      st_transform(4326) |>
      st_make_valid() |>
      st_intersection(st_union(basin_sf))
    # Normalise column names to lowercase after intersection
    names(rivers_sf) <- tolower(names(rivers_sf))
    # Bug 2 fix: check column existence before mutating, avoid bare `.`
    rivers_sf <- rivers_sf |>
      mutate(
        STREAM_ORDER = if ("ord_stra" %in% names(rivers_sf)) ord_stra else NA_integer_,
        GNIS_NAME    = if ("river"    %in% names(rivers_sf)) river    else NA_character_
      )
    cat(sprintf("  Rivers: %d segments\n", nrow(rivers_sf)))
  }, error = function(e) {
    message("  Rivers failed (", conditionMessage(e), ") — skipping.")
    rivers_sf <<- sf::st_sf(
      STREAM_ORDER = integer(0), GNIS_NAME = character(0),
      geometry = sf::st_sfc(crs = 4326)
    )
  })
  
  # NaturalEarth lakes — Bug 2 fix: avoid `.` reference in mutate
  tryCatch({
    ne_lakes_raw <- ne_download(
      scale = 10, type = "lakes", category = "physical",
      returnclass = "sf"
    ) |>
      st_transform(4326) |>
      st_make_valid() |>
      st_intersection(st_union(basin_sf))
    # Rename columns after intersection — no bare `.` needed
    has_name_col  <- "name"  %in% names(ne_lakes_raw)
    has_area_col  <- "area"  %in% names(ne_lakes_raw)
    lakes_sf <- ne_lakes_raw |>
      mutate(
        GNIS_NAME_1 = if (has_name_col) name else NA_character_,
        AREA_HA     = if (has_area_col) area else NA_real_
      )
    cat(sprintf("  NaturalEarth lakes: %d polygons\n", nrow(lakes_sf)))
  }, error = function(e) {
    message("  Lakes failed (", conditionMessage(e), ") — skipping.")
    lakes_sf <<- sf::st_sf(
      GNIS_NAME_1 = character(0), AREA_HA = numeric(0),
      geometry = sf::st_sfc(crs = 4326)
    )
  })
}

## ── 5c. Stream-order line widths ─────────────────────────────────────────────
# Bug 3 fix: guard against empty sf (0-row data frame crashes $<- assignment)
if (nrow(rivers_sf) == 0) {
  rivers_sf$line_width <- numeric(0)
} else if ("STREAM_ORDER" %in% names(rivers_sf) && !all(is.na(rivers_sf$STREAM_ORDER))) {
  rivers_sf <- rivers_sf |>
    mutate(line_width = case_when(
      STREAM_ORDER >= 8 ~ 1.10,
      STREAM_ORDER >= 7 ~ 0.75,
      STREAM_ORDER >= 6 ~ 0.50,
      TRUE              ~ 0.35
    ))
} else {
  rivers_sf$line_width <- 0.5
}

# ── 6. KEY LABELS (major lakes & rivers to annotate) ─────────────────────────────
# Centroids of the largest lakes for labelling
if ("AREA_HA" %in% names(lakes_sf)) {
  top_lakes <- lakes_sf |>
    arrange(desc(AREA_HA)) |>
    slice_head(n = 8) |>
    filter(!is.na(GNIS_NAME_1), GNIS_NAME_1 != "") |>
    # Keep only the largest polygon per lake name to avoid duplicate labels
    # when st_intersection() splits a lake across the basin boundary
    group_by(GNIS_NAME_1) |>
    slice_max(AREA_HA, n = 1, with_ties = FALSE) |>
    ungroup() |>
    st_centroid()
} else {
  top_lakes <- NULL
}

# ── 7. POINTS OF INTEREST ───────────────────────────────────────────────────────
cat("── Building points of interest ──\n")
poi_df <- data.frame(
  name = c("Prince George", "Vanderhoof",
           "Kenney Dam", "Skins Lake Spillway"),
  lon  = c(-122.7497, -124.0116,
           -125.020,  -125.970),
  lat  = c(  53.9171,   54.0117,
             53.730,    53.771),
  type = c("City",  "Town",
           "Dam",   "Spillway")
)
poi_sf <- st_as_sf(poi_df, coords = c("lon", "lat"), crs = 4326)

cities_sf    <- filter(poi_sf, type %in% c("City", "Town"))
dams_sf      <- filter(poi_sf, type == "Dam")
spillways_sf <- filter(poi_sf, type == "Spillway")

# ── 8. COLOUR PALETTES ──────────────────────────────────────────────────────────
n_col <- 512
elev_pal <- c(
  colorRampPalette(c("#2d6a2d", "#4a9b4a", "#8ab87a"))(round(n_col * 0.18)),
  colorRampPalette(c("#c8d878", "#d4b050", "#b88030"))(round(n_col * 0.27)),
  colorRampPalette(c("#b88030", "#c8a060", "#d8b888"))(round(n_col * 0.25)),
  colorRampPalette(c("#d8b888", "#ece4c8", "#f8f2e0", "#ffffff"))(round(n_col * 0.30))
)
n_col_actual <- length(elev_pal)

water_fill  <- "#5baee0"
water_line  <- "#3a85c0"
river_col   <- "#3a85c0"
bg_water    <- "#c8dff0"   # ocean / background
basin_col   <- "grey15"
city_fill   <- "#e74c3c"
town_fill   <- "#f39c12"
dam_shape   <- 22  # square for dam
city_shape  <- 21  # circle for city
town_shape  <- 21

# ── 9. RASTER DATA FRAMES ──────────────────────────────────────────────────────
cat("── Preparing raster data frames ──\n")
dem_df  <- as.data.frame(dem,  xy = TRUE, na.rm = TRUE)
hill_df <- as.data.frame(hill, xy = TRUE, na.rm = FALSE) |>
  mutate(shadow = ifelse(is.na(hillshade), NA_real_,
                         pmax(0, pmin(1, (1 - hillshade) * 0.58))))

# ── 10. INSET LOCATOR MAP ───────────────────────────────────────────────────────
cat("── Building inset locator map ──\n")

bc_sf     <- ne_states(country = "Canada", returnclass = "sf") |>
  filter(name == "British Columbia") |>
  st_transform(4326)
canada_sf <- ne_countries(country = "Canada", scale = "medium", returnclass = "sf") |>
  st_transform(4326)
usa_sf    <- ne_countries(country = "United States of America",
                          scale = "medium", returnclass = "sf") |>
  st_transform(4326)

inset <- ggplot() +
  geom_sf(data = usa_sf,    fill = "gray88", color = "gray60", linewidth = 0.25) +
  geom_sf(data = canada_sf, fill = "gray82", color = "gray55", linewidth = 0.25) +
  geom_sf(data = bc_sf,     fill = "gray65", color = "gray35", linewidth = 0.40) +
  geom_sf(data = basin_sf,  fill = "#e74c3c", color = "#a93226",
          linewidth = 0.5, alpha = 0.80) +
  coord_sf(xlim = c(-140, -113), ylim = c(47, 61), expand = FALSE) +
  theme_void(base_size = 7) +
  theme(
    panel.background = element_rect(fill = bg_water,  color = NA),
    panel.border     = element_rect(fill = NA, color = "black", linewidth = 0.6),
    plot.background  = element_rect(fill = "white",   color = NA)
  )

# ── 11. MAIN MAP ─────────────────────────────────────────────────────────────────
cat("── Rendering main map ──\n")

main_map <- ggplot() +
  
  # ── Background (ocean/outside basin) ──────────────────────────────────────────
  theme(panel.background = element_rect(fill = bg_water, color = NA)) +
  
  # ── Elevation ────────────────────────────────────────────────────────────────
  geom_raster(data = dem_df,
              aes(x = x, y = y, fill = elevation),
              interpolate = TRUE) +
  scale_fill_gradientn(
    colours  = elev_pal,
    name     = "Elevation (m)",
    limits   = elev_range,
    breaks   = pretty(elev_range, n = 6),
    na.value = bg_water,
    guide    = guide_colorbar(
      barwidth       = unit(0.55, "cm"),
      barheight      = unit(7.5,  "cm"),
      title.position = "top",
      title.hjust    = 0.5,
      ticks.colour   = "black",
      frame.colour   = "black",
      frame.linewidth = 0.4
    )
  ) +
  
  # ── Hillshade overlay (transparency-blend shadow) ─────────────────────────────
  ggnewscale::new_scale_fill() +
  geom_raster(data = filter(hill_df, !is.na(shadow)),
              aes(x = x, y = y, alpha = shadow),
              fill = "black", interpolate = TRUE) +
  scale_alpha_identity(guide = "none") +
  
  # ── River polygons (double-banked channels — wider rivers filled blue) ─────────
  {if (!is.null(rivers_poly_sf) && nrow(rivers_poly_sf) > 0)
    geom_sf(data = rivers_poly_sf,
            fill = water_fill, color = water_line,
            linewidth = 0.10, inherit.aes = FALSE)
  } +
  
  # ── Lakes ─────────────────────────────────────────────────────────────────────
  geom_sf(data  = lakes_sf,
          fill  = water_fill, color = water_line,
          linewidth = 0.10, inherit.aes = FALSE) +
  
  # ── Rivers (width by stream order) ───────────────────────────────────────────
  geom_sf(data      = rivers_sf,
          aes(linewidth = line_width),
          color     = river_col,
          inherit.aes = FALSE) +
  scale_linewidth_identity(guide = "none") +
  
  # ── Basin boundary ────────────────────────────────────────────────────────────
  geom_sf(data  = basin_sf,
          fill  = NA, color = basin_col,
          linewidth = 1.3, inherit.aes = FALSE) +
  
  # ── Towns ─────────────────────────────────────────────────────────────────────
  geom_sf(data  = filter(cities_sf, type == "Town"),
          shape = town_shape, fill = town_fill, color = "black",
          size  = 2.5, stroke = 0.8, inherit.aes = FALSE) +
  
  # ── Prince George (City) ──────────────────────────────────────────────────────
  geom_sf(data  = filter(cities_sf, type == "City"),
          shape = city_shape, fill = city_fill, color = "white",
          size  = 4.5, stroke = 1.2, inherit.aes = FALSE) +
  
  # ── Kenney Dam ────────────────────────────────────────────────────────────────
  geom_sf(data  = dams_sf,
          shape = 22, fill = "#e67e22", color = "black",
          size  = 3.5, stroke = 1.0, inherit.aes = FALSE) +
  
  # ── Skins Lake Spillway ───────────────────────────────────────────────────────
  geom_sf(data  = spillways_sf,
          shape = 24, fill = "#9b59b6", color = "black",
          size  = 3.2, stroke = 1.0, inherit.aes = FALSE) +
  
  # ── POI labels (white fill — constant, never mapped, no legend) ──────────────
  ggrepel::geom_label_repel(
    data          = poi_df,
    aes(x = lon, y = lat, label = name),
    fontface      = ifelse(poi_df$type == "City", "bold", "plain"),
    fill          = alpha("white", 0.88),
    size          = 3.0,
    label.size    = 0.18,
    label.r       = unit(0.11, "lines"),
    color         = "grey15",
    box.padding   = 0.65,
    point.padding = 0.35,
    min.segment.length = 0,
    max.overlaps  = 25,
    force         = 4,
    seed          = 42,
    show.legend   = FALSE,
    inherit.aes   = FALSE
  ) +
  
  # ── Lake labels (blue fill — constant, never mapped, no legend) ───────────────
  {if (!is.null(top_lakes) && nrow(top_lakes) > 0) {
    lake_pts <- data.frame(
      lon  = sf::st_coordinates(top_lakes)[, 1],
      lat  = sf::st_coordinates(top_lakes)[, 2],
      name = top_lakes$GNIS_NAME_1,
      stringsAsFactors = FALSE
    )
    ggrepel::geom_label_repel(
      data          = lake_pts,
      aes(x = lon, y = lat, label = name),
      fontface      = "italic",
      fill          = alpha(water_fill, 0.75),
      size          = 2.6,
      label.size    = 0.15,
      label.r       = unit(0.10, "lines"),
      color         = "grey10",
      box.padding   = 0.50,
      point.padding = 0.30,
      min.segment.length = 0,
      max.overlaps  = 20,
      force         = 4,
      seed          = 77,
      show.legend   = FALSE,
      inherit.aes   = FALSE
    )
  }} +
  
  # ── Coordinate system ────────────────────────────────────────────────────────
  coord_sf(crs = 4326, expand = FALSE) +
  
  # ── Scale bar ─────────────────────────────────────────────────────────────────
  annotation_scale(
    location   = "br",
    width_hint = 0.22,
    style      = "ticks",
    text_cex   = 0.75,
    pad_x      = unit(0.35, "cm"),
    pad_y      = unit(0.35, "cm")
  ) +
  
  # ── North arrow ───────────────────────────────────────────────────────────────
  annotation_north_arrow(
    location = "tr",
    pad_x    = unit(0.35, "cm"),
    pad_y    = unit(0.35, "cm"),
    height   = unit(1.1, "cm"),
    width    = unit(0.85, "cm"),
    style    = north_arrow_orienteering(
      fill     = c("grey15", "white"),
      text_col = "grey15",
      text_size = 8,
      line_col  = "grey15"
    )
  ) +
  
  # ── Custom legend for points ──────────────────────────────────────────────────
  # (manual legend via annotate - cleaner than scale overrides)
  
  # ── Labels ────────────────────────────────────────────────────────────────────
  labs(
    title    = "Nechako River Basin",
    x        = "Longitude (°W)",
    y        = "Latitude (°N)",
    caption  = paste0(
      "Elevation: AWS Terrain Tiles via elevatr (z = 8, ~300 m)",
      "  |  Rivers & Lakes: BC Freshwater Atlas via bcdata (or OpenStreetMap)",
      "\nHillshade: azimuth 315°, altitude 40°",
      "  |  Basin boundary: Nechako watershed (WGS84 / EPSG:4326)"
    )
  ) +
  
  # ── Theme ─────────────────────────────────────────────────────────────────────
  theme_minimal(base_size = 12, base_family = "sans") +
  theme(
    plot.background    = element_rect(fill = "white",  color = NA),
    panel.background   = element_rect(fill = bg_water, color = NA),
    panel.grid.major   = element_line(color = alpha("white", 0.55),
                                      linetype = "dashed", linewidth = 0.25),
    panel.border       = element_rect(fill = NA, color = "grey20", linewidth = 0.7),
    plot.title         = element_text(face = "bold", size = 20, hjust = 0.5,
                                      margin = margin(b = 10)),
    plot.caption       = element_text(size = 7.5, color = "grey50", hjust = 0,
                                      margin = margin(t = 6)),
    axis.title         = element_text(size = 9.5, color = "grey30"),
    axis.text          = element_text(size = 8.5, color = "grey30"),
    legend.position    = "right",
    legend.justification = "center",
    legend.title       = element_text(size = 9.5, face = "bold", hjust = 0.5),
    legend.text        = element_text(size = 8.5),
    legend.background  = element_rect(fill = alpha("white", 0.85), color = NA),
    legend.key.height  = unit(0.4, "cm"),
    plot.margin        = margin(8, 8, 8, 8)
  )

# ── 12. COMBINE MAIN MAP + INSET ─────────────────────────────────────────────────
cat("── Assembling final layout ──\n")

final_map <- main_map +
  inset_element(
    inset,
    left   = 0.00,   # fraction of main panel
    bottom = 0.65,
    right  = 0.195,
    top    = 0.985,
    align_to = "plot"
  )

# ── 14. SAVE ─────────────────────────────────────────────────────────────────────
cat("── Saving outputs ──\n")

pdf_path <- file.path(out_dir, "nechako_overview_map.pdf")
ggsave(pdf_path, plot = final_map,
       width = 13, height = 10.5, units = "in",
       device = cairo_pdf)          # cairo_pdf → better font/transparency rendering
cat(sprintf("  PDF: %s\n", pdf_path))

png_path <- file.path(out_dir, "nechako_overview_map.png")
ggsave(png_path, plot = final_map,
       width = 13, height = 10.5, units = "in",
       dpi = 300, device = "png", bg = "white")
cat(sprintf("  PNG: %s\n", png_path))

cat("\n══ Done! ══\n")
cat("  nechako_overview_map.pdf  — vector, print-quality\n")
cat("  nechako_overview_map.png  — 300 dpi raster\n")