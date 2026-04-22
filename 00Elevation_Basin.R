####################################################################################
# NECHAKO RIVER BASIN — COMPREHENSIVE OVERVIEW MAP  (v2)
# Matches the reference publication figure (journal_pwat_0000263_g002)
#
# NEW in v2:
#   • Fraser watershed background (gray fill)
#   • Paved roads (BC Digital Road Atlas / NaturalEarth fallback)
#   • First Nation communities (yellow dots) — 19 communities
#   • Full municipality set (red dots) — inside + outside basin
#   • Expanded water body labels (Stuart Lake, Trembleur Lake, Francois Lake,
#     Stuart River, Fraser River, Nechako Reservoir (Ootsa Lake))
#   • Comprehensive legend inset replacing ggplot2 default guide
#   • Updated inset locator showing Fraser watershed context
#   • Light-gray panel background (approximates ESRI World Light Gray Base)
#
# Data sources:
#   • DEM       : AWS Terrain Tiles via elevatr (z = 8, ~300 m)
#   • Rivers    : BC Freshwater Atlas via bcdata  [PRIMARY]
#   • Lakes     : BC Freshwater Atlas via bcdata  [PRIMARY]
#   • Fraser WS : BC Freshwater Atlas watershed groups via bcdata
#   • Roads     : BC Digital Road Atlas via bcdata [PRIMARY]
#   • Fallback  : HydroSHEDS (rivers) + NaturalEarth lakes + NaturalEarth roads
#
# Packages required:
#   terra, sf, ggplot2, elevatr, bcdata, geodata,
#   ggspatial, ggrepel, ggnewscale, rnaturalearth, rnaturalearthdata,
#   patchwork, dplyr, scales
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
basin_shp <- "D:/Nechako_Drought/Nechako/Spatial/nechakoBound_dissolve.shp"
out_dir   <- "D:/Nechako_Drought/Nechako/elevation"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

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

# Fixed map extent — wide enough for Kitimat (W), Mackenzie (NE), Quesnel (S)
map_xlim <- c(-130.0, -121.0)
map_ylim <- c(  52.3,   56.3)

# Expanded query bbox in BC Albers for road + Fraser watershed queries
query_bbox_3005 <- st_as_sfc(
  st_bbox(c(xmin = map_xlim[1] - 0.5, xmax = map_xlim[2] + 0.5,
            ymin = map_ylim[1] - 0.5, ymax = map_ylim[2] + 0.5),
          crs = 4326)
) |> st_transform(3005)

# ── 3. DEM ───────────────────────────────────────────────────────────────────────
cat("── Downloading DEM via elevatr (z = 8) ──\n")
dir.create(tempdir(), showWarnings = FALSE, recursive = TRUE)
options(timeout = 300)
dem_raw <- get_elev_raster(locations = basin_sf, z = 8, clip = "bbox",
                           src = "aws", neg_to_na = FALSE, verbose = FALSE)
dem <- rast(dem_raw)
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
      filter(STREAM_ORDER >= 6) |>
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
      STREAM_ORDER >= 6 ~ 0.50, TRUE              ~ 0.35))
} else {
  rivers_sf$line_width <- 0.5
}

# ── 5.5 FRASER WATERSHED BOUNDARY ────────────────────────────────────────────────
cat("── Fetching Fraser watershed boundary ──\n")
fraser_sf <- NULL
if (has_bcdata) {
  # Try 1: major watersheds record
  fraser_sf <- tryCatch({
    cat("  Trying freshwater-atlas-major-watersheds ...\n")
    bcdc_query_geodata("freshwater-atlas-major-watersheds") |>
      filter(MAJOR_WATERSHED_SYSTEM == "FRASER RIVER") |>
      collect() |> st_transform(4326) |> st_union()
  }, error = function(e1) {
    # Try 2: union of known Fraser drainage watershed group codes
    tryCatch({
      cat("  Trying watershed-groups union ...\n")
      bcdc_query_geodata("freshwater-atlas-watershed-groups") |>
        filter(WATERSHED_GROUP_CODE %in% c(
          "BRID", "CANA", "CHIL", "CROL", "DEAD", "ENCL", "FRAN",
          "GREE", "HARR", "HORS", "ILEC", "LCHL", "LILL", "LNIC",
          "MAST", "NAZO", "NECH", "NICL", "QUES", "SALA", "SALM",
          "SARA", "SHUA", "STUL", "THOM", "WIDG"
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
cat("── Fetching paved roads ──\n")
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

# ── 5.7 FRASER RIVER LINE (full visible extent — south of basin AND upstream N) ──
# The existing rivers_sf is clipped to the Nechako basin, so the Fraser main
# stem (which runs outside the basin) must be fetched separately.
# The river flows: source in Rocky Mts (NE) → northwest → Prince George → south.
# We need both the southern reach (Prince George → Quesnel) AND the northern/
# upstream reach (Prince George → Fort George / Giscome area along the E edge).
# A dedicated INTERSECTS bbox extended east to lon = -119.5 captures those
# upstream segments that would otherwise fall outside query_bbox_3005.
cat("── Fetching Fraser River line ──\n")

fraser_river_sf <- NULL

# Dedicated query bbox: full map width + extra east to catch upstream bend
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
      filter(STREAM_ORDER >= 8) |>          # 8+ captures full main stem
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
  # Fallback: NaturalEarth rivers_lake_centerlines (single named polyline)
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

# Smooth the Fraser River geometry:
# Many short FWA segments rendered with a thick line produce jagged "spikes".
# Fix: union all segments into one geometry, simplify in projected metres (BC
# Albers), then reproject back to WGS84.
if (!is.null(fraser_river_sf) && nrow(fraser_river_sf) > 0) {
  fraser_river_sf <- fraser_river_sf |>
    st_transform(3005) |>
    st_union() |>                              # merge all segments
    st_simplify(dTolerance = 150,              # 150 m tolerance — smooth but accurate
                preserveTopology = TRUE) |>
    st_transform(4326) |>
    st_as_sf() |>
    rename(geometry = x)                       # st_as_sf names the col "x" after union
  cat(sprintf("  Fraser River after smoothing: %d feature(s)\n",
              nrow(fraser_river_sf)))
}

# ── 6. KEY LAKES FOR LABELLING ───────────────────────────────────────────────────
if ("AREA_HA" %in% names(lakes_sf)) {
  top_lakes <- lakes_sf |>
    arrange(desc(AREA_HA)) |>
    slice_head(n = 12) |>
    filter(!is.na(GNIS_NAME_1), GNIS_NAME_1 != "") |>
    group_by(GNIS_NAME_1) |>
    slice_max(AREA_HA, n = 1, with_ties = FALSE) |>
    ungroup() |>
    st_centroid()
} else {
  top_lakes <- NULL
}

# ── 7. POINTS OF INTEREST ───────────────────────────────────────────────────────
cat("── Building points of interest ──\n")

# ── 7a. Municipalities — inside basin only ────────────────────────────────────────
muni_inside_df <- data.frame(
  name = c("Vanderhoof", "Burns Lake", "Fraser Lake",
           "Fort St. James", "Granisle", "Houston"),
  lon  = c(-124.0116, -125.759, -124.847, -124.253, -126.220, -126.668),
  lat  = c(  54.0117,   54.232,   54.054,   54.434,   54.920,   54.398),
  type = "Municipality", stringsAsFactors = FALSE
)

# Prince George — outside basin but included by request (city symbol, bold label)
pg_df <- data.frame(
  name = "Prince George", lon = -122.7497, lat = 53.9171,
  type = "City", stringsAsFactors = FALSE
)
pg_sf <- st_as_sf(pg_df, coords = c("lon","lat"), crs = 4326)

# ── 7b. Infrastructure (Dam, Spillway, Powerhouse) ───────────────────────────────
# Kemano Powerhouse is outside the basin but kept by request.
infra_df <- data.frame(
  name    = c("Kenney Dam",  "Skins Lake Spillway", "Kemano Powerhouse"),
  lon     = c(-125.020,      -125.970,               -127.882),
  lat     = c(  53.730,        53.771,                  53.558),
  type    = c("Dam",          "Spillway",              "Powerhouse"),
  nudge_x = c( 0.00,          -0.45,                    0.00),
  nudge_y = c( 0.00,           0.20,                    0.00),
  stringsAsFactors = FALSE
)

# ── 7c. Combined label data frame (municipalities + Prince George + infrastructure)
poi_labels_df <- bind_rows(
  muni_inside_df |> mutate(nudge_x = 0, nudge_y = 0),
  pg_df          |> mutate(nudge_x = -1.20, nudge_y = 0.25),
  infra_df       |> select(name, lon, lat, type, nudge_x, nudge_y)
)

# sf objects
muni_inside_sf <- st_as_sf(muni_inside_df, coords = c("lon","lat"), crs = 4326)
dams_sf        <- st_as_sf(filter(infra_df, type == "Dam"),
                           coords = c("lon","lat"), crs = 4326)
spillways_sf   <- st_as_sf(filter(infra_df, type == "Spillway"),
                           coords = c("lon","lat"), crs = 4326)
powerhouses_sf <- st_as_sf(filter(infra_df, type == "Powerhouse"),
                           coords = c("lon","lat"), crs = 4326)

# ── 7.5 MANUAL WATER BODY LABELS ─────────────────────────────────────────────────
# lon/lat = label anchor; angle = text rotation (degrees)
manual_labels_df <- data.frame(
  lon   = c(-126.20, -124.35, -124.25, -124.90, -125.20, -122.82),
  lat   = c(  53.22,   53.92,   54.55,   54.05,   54.40,   53.75),
  label = c("Nechako Reservoir\n(Ootsa Lake)",
            "Nechako River",
            "Stuart Lake",
            "Francois Lake",
            "Trembleur Lake",
            "Fraser River"),
  angle = c(  0,  0, 75,  0,  0, 60),
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

water_fill    <- "#7ec8e3"  # waterbody fill (light blue)
water_line    <- "#4a90d9"  # waterbody border
river_col     <- "#4a90d9"  # river lines
bg_land       <- "#eeeeee"  # panel background ~ ESRI World Light Gray Base
fraser_fill   <- "#d0d0d0"  # Fraser watershed fill (medium gray)
fraser_line   <- "#aaaaaa"  # Fraser watershed border
basin_col     <- "grey10"   # Nechako boundary outline
muni_col      <- "#e74c3c"  # Municipality (red)
road_col      <- "#888888"  # Paved road
dam_fill      <- "#e67e22"  # Kenney Dam (orange)
spillway_fill <- "#9b59b6"  # Skins Lake Spillway (purple)
power_fill    <- "#f1c40f"  # Kemano Powerhouse (yellow)

# ── 9. RASTER DATA FRAMES ──────────────────────────────────────────────────────
cat("── Preparing raster data frames ──\n")
dem_df  <- as.data.frame(dem, xy = TRUE, na.rm = TRUE)
hill_df <- as.data.frame(hill, xy = TRUE, na.rm = FALSE) |>
  mutate(shadow = ifelse(is.na(hillshade), NA_real_,
                         pmax(0, pmin(1, (1 - hillshade) * 0.58))))

# ── 10. LEGEND INSET ─────────────────────────────────────────────────────────────
cat("── Building legend inset ──\n")

# Abstract y-axis: 10 categorical rows (spacing 0.95) above a separator,
# then "Elevation (m)" title and colour swatch below.
# Row centres (top → bottom):
#   City(13.1) Municipality(12.15) Nechako WS(11.2) Fraser WS(10.25)
#   Waterbodies(9.3) Lake label(8.35) Paved road(7.4)
#   Kenney Dam(6.45) Skins Spillway(5.5) Kemano Powerhouse(4.55)
#   separator(4.0) Elevation title(3.75) swatch(0.30–3.40)

leg_ymin <- 0.30
leg_ymax <- 3.40
swatch_xctr  <- 0.375
swatch_width <- 0.55
n_tiles      <- 60
tile_h       <- (leg_ymax - leg_ymin) / n_tiles

elev_grad_df <- data.frame(
  x    = swatch_xctr,
  y    = seq(leg_ymin + tile_h / 2, leg_ymax - tile_h / 2, length.out = n_tiles),
  fill = seq(0, 1, length.out = n_tiles)
)

legend_plot <- ggplot() +
  
  # ── City / Municipality — two dot sizes on one row ───────────────────────
  # Large dot (city) at x=0.22, small dot (municipality) at x=0.52,
  # single label to the right.
  annotate("point", x=0.22, y=12.60,
           shape=21, fill=muni_col, color="grey20", size=3.8, stroke=0.9) +
  annotate("point", x=0.52, y=12.60,
           shape=21, fill=muni_col, color="grey20", size=2.5, stroke=0.7) +
  annotate("text",  x=0.78, y=12.60, hjust=0, size=2.45, color="grey10",
           label="City / Municipality") +
  
  # ── Nechako Watershed ────────────────────────────────────────────────────
  annotate("rect",  xmin=0.10, xmax=0.65, ymin=11.32, ymax=12.08,
           fill=NA, color=basin_col, linewidth=0.85) +
  annotate("text",  x=0.78, y=11.70, hjust=0, size=2.45, color="grey10",
           label="Nechako Watershed") +
  
  # ── Fraser Watershed ─────────────────────────────────────────────────────
  annotate("rect",  xmin=0.10, xmax=0.65, ymin=10.37, ymax=11.13,
           fill=fraser_fill, color=fraser_line, linewidth=0.4) +
  annotate("text",  x=0.78, y=10.75, hjust=0, size=2.45, color="grey10",
           label="Fraser Watershed") +
  
  # ── Waterbodies ──────────────────────────────────────────────────────────
  annotate("rect",  xmin=0.10, xmax=0.65, ymin=9.42, ymax=10.18,
           fill=water_fill, color=water_line, linewidth=0.4) +
  annotate("text",  x=0.78, y=9.80, hjust=0, size=2.45, color="grey10",
           label="Waterbodies") +
  
  # ── Lake name label (italic blue box) ────────────────────────────────────
  annotate("label", x=0.375, y=8.85,
           label="Lake", fontface="italic",
           fill=alpha(water_fill, 0.70), color="grey10",
           size=2.0, label.size=0.12, label.r=unit(0.10,"lines"),
           label.padding=unit(0.15,"lines")) +
  annotate("text",  x=0.78, y=8.85, hjust=0, size=2.45, color="grey10",
           label="Lake name") +
  
  # ── Paved road ───────────────────────────────────────────────────────────
  annotate("segment", x=0.10, xend=0.65, y=7.90, yend=7.90,
           color=road_col, linewidth=0.75) +
  annotate("text",  x=0.78, y=7.90, hjust=0, size=2.45, color="grey10",
           label="Paved road") +
  
  # ── Kenney Dam (orange filled square) ────────────────────────────────────
  annotate("point", x=0.375, y=6.95,
           shape=22, fill=dam_fill, color="black", size=3.0, stroke=0.8) +
  annotate("text",  x=0.78, y=6.95, hjust=0, size=2.45, color="grey10",
           label="Kenney Dam") +
  
  # ── Skins Lake Spillway (purple upward triangle) ──────────────────────────
  annotate("point", x=0.375, y=6.00,
           shape=24, fill=spillway_fill, color="black", size=2.8, stroke=0.8) +
  annotate("text",  x=0.78, y=6.00, hjust=0, size=2.45, color="grey10",
           label="Skins Lake Spillway") +
  
  # ── Kemano Powerhouse (yellow diamond) ───────────────────────────────────
  annotate("point", x=0.375, y=5.05,
           shape=23, fill=power_fill, color="black", size=3.0, stroke=0.8) +
  annotate("text",  x=0.78, y=5.05, hjust=0, size=2.45, color="grey10",
           label="Kemano Powerhouse") +
  
  # ── Separator ────────────────────────────────────────────────────────────
  annotate("segment", x=0.05, xend=2.60, y=4.50, yend=4.50,
           color="grey65", linewidth=0.28) +
  
  # ── Elevation swatch ─────────────────────────────────────────────────────
  annotate("text",  x=0.375, y=4.28, hjust=0.5, vjust=1, size=2.45,
           fontface="bold", color="grey10", label="Elevation (m)") +
  geom_tile(data=elev_grad_df, aes(x=x, y=y, fill=fill),
            width=swatch_width, height=tile_h) +
  scale_fill_gradientn(colours=elev_pal, guide="none") +
  annotate("rect",
           xmin = swatch_xctr - swatch_width / 2,
           xmax = swatch_xctr + swatch_width / 2,
           ymin = leg_ymin, ymax = leg_ymax,
           fill=NA, color="black", linewidth=0.4) +
  annotate("text",  x=0.70, y=leg_ymax,
           label=paste0(round(elev_range[2]), " m"),
           hjust=0, vjust=1.0, size=2.1, color="grey15") +
  annotate("text",  x=0.70, y=leg_ymin,
           label=paste0(round(elev_range[1]), " m"),
           hjust=0, vjust=0.0, size=2.1, color="grey15") +
  
  # ── Axes ─────────────────────────────────────────────────────────────────
  xlim(0.05, 2.65) +
  ylim(0.05, 13.20) +
  theme_void() +
  theme(
    plot.background = element_rect(fill=alpha("white", 0.92),
                                   color="grey35", linewidth=0.45),
    plot.margin     = margin(4, 5, 4, 5)
  )

# ── 11. INSET LOCATOR MAP ────────────────────────────────────────────────────────
cat("── Building inset locator map ──\n")

# Keep ONLY the single largest polygon from BC = the mainland.
# Avoids the dot-pronoun issue inside mutate() by assigning first,
# then indexing with which.max(st_area(...)).
bc_sf_inset <- ne_states(country="Canada", returnclass="sf") |>
  filter(name=="British Columbia") |>
  st_transform(4326) |>
  st_cast("POLYGON")
bc_sf_inset <- bc_sf_inset[which.max(st_area(bc_sf_inset)), ]

# scale="small" (1:110M) — no small island fragments by design
canada_sf <- ne_countries(country="Canada",
                          scale="small", returnclass="sf") |>
  st_transform(4326)
usa_sf    <- ne_countries(country="United States of America",
                          scale="small", returnclass="sf") |>
  st_transform(4326)

# Simplify basin_sf for the inset (no need for full resolution at small scale)
basin_sf_inset <- basin_sf |>
  st_transform(3005) |>
  st_simplify(dTolerance = 2000, preserveTopology = TRUE) |>
  st_transform(4326)

inset <- ggplot() +
  geom_sf(data=usa_sf,       fill="gray92", color="gray70", linewidth=0.20) +
  geom_sf(data=canada_sf,    fill="gray88", color="gray65", linewidth=0.20) +
  geom_sf(data=bc_sf_inset,  fill="gray75", color="gray45", linewidth=0.35) +
  geom_sf(data=basin_sf_inset, fill=NA, color="#c0392b",
          linewidth=0.85, alpha=0.9) +
  annotate("text", x=-135.5, y=54.0, label="British\nColumbia",
           size=1.9, color="grey20", fontface="italic", lineheight=0.85) +
  coord_sf(xlim=c(-140,-113), ylim=c(47,61), expand=FALSE) +
  theme_void(base_size=7) +
  theme(
    panel.background = element_rect(fill="#c8dff0", color=NA),
    panel.border     = element_rect(fill=NA, color="grey25", linewidth=0.6),
    plot.background  = element_rect(fill="white", color=NA)
  )

# ── 12. MAIN MAP ─────────────────────────────────────────────────────────────────
cat("── Rendering main map ──\n")

main_map <- ggplot() +
  
  # Background: light gray land (ESRI World Light Gray Base approximation)
  theme(panel.background = element_rect(fill=bg_land, color=NA)) +
  
  # Fraser watershed fill (behind elevation)
  {if (!is.null(fraser_sf))
    geom_sf(data=fraser_sf,
            fill=fraser_fill, color=fraser_line,
            linewidth=0.3, inherit.aes=FALSE)} +
  
  # Elevation raster (masked to Nechako basin)
  geom_raster(data=dem_df,
              aes(x=x, y=y, fill=elevation),
              interpolate=TRUE) +
  scale_fill_gradientn(
    colours  = elev_pal, limits = elev_range,
    breaks   = pretty(elev_range, n = 6),
    name     = "Elevation (m)",
    na.value = NA,
    guide    = guide_colorbar(
      barwidth       = unit(0.5, "cm"),
      barheight      = unit(6.5, "cm"),
      title.position = "top",
      title.hjust    = 0.5,
      ticks.colour   = "black",
      frame.colour   = "black",
      frame.linewidth = 0.4
    )
  ) +
  
  # Hillshade overlay
  ggnewscale::new_scale_fill() +
  geom_raster(data=filter(hill_df, !is.na(shadow)),
              aes(x=x, y=y, alpha=shadow),
              fill="black", interpolate=TRUE) +
  scale_alpha_identity(guide="none") +
  
  # River channel polygons (wide rivers as blue fill)
  {if (!is.null(rivers_poly_sf) && nrow(rivers_poly_sf) > 0)
    geom_sf(data=rivers_poly_sf,
            fill=water_fill, color=water_line,
            linewidth=0.10, inherit.aes=FALSE)} +
  
  # Lakes
  geom_sf(data=lakes_sf,
          fill=water_fill, color=water_line,
          linewidth=0.10, inherit.aes=FALSE) +
  
  # River lines (width by stream order)
  geom_sf(data=rivers_sf,
          aes(linewidth=line_width),
          color=river_col, inherit.aes=FALSE) +
  scale_linewidth_identity(guide="none") +
  
  # Paved roads
  {if (!is.null(roads_sf) && nrow(roads_sf) > 0)
    geom_sf(data=roads_sf,
            color=road_col, linewidth=0.25,
            inherit.aes=FALSE)} +
  
  # Fraser River main stem (outside Nechako basin)
  {if (!is.null(fraser_river_sf) && nrow(fraser_river_sf) > 0)
    geom_sf(data=fraser_river_sf,
            color=water_line, linewidth=0.65,
            lineend="round", linejoin="round",
            inherit.aes=FALSE)} +
  
  # Nechako basin boundary (thick black outline)
  geom_sf(data=basin_sf,
          fill=NA, color=basin_col,
          linewidth=1.2, inherit.aes=FALSE) +
  
  # Infrastructure markers (Dam, Spillway, Powerhouse)
  geom_sf(data=dams_sf,
          shape=22, fill=dam_fill, color="black",
          size=3.2, stroke=0.9, inherit.aes=FALSE) +
  geom_sf(data=spillways_sf,
          shape=24, fill=spillway_fill, color="black",
          size=3.0, stroke=0.9, inherit.aes=FALSE) +
  geom_sf(data=powerhouses_sf,
          shape=23, fill=power_fill, color="black",
          size=3.2, stroke=0.9, inherit.aes=FALSE) +
  
  # Municipalities — inside basin only
  geom_sf(data=muni_inside_sf,
          shape=21, fill=muni_col, color="grey20",
          size=3.5, stroke=0.8, inherit.aes=FALSE) +
  
  # Prince George — outside basin, larger city symbol
  geom_sf(data=pg_sf,
          shape=21, fill=muni_col, color="grey20",
          size=4.5, stroke=1.0, inherit.aes=FALSE) +
  
  # Water body italic labels (Nechako R., Stuart Lake, Francois Lake, etc.)
  geom_text(
    data=manual_labels_df,
    aes(x=lon, y=lat, label=label, angle=angle),
    fontface="italic", colour="#1a6699",
    size=2.7, lineheight=0.85, inherit.aes=FALSE
  ) +
  
  # Municipality + infrastructure labels
  ggrepel::geom_label_repel(
    data          = poi_labels_df,
    aes(x = lon, y = lat, label = name),
    fontface      = ifelse(poi_labels_df$name == "Prince George", "bold", "plain"),
    fill          = alpha("white", 0.85),
    size          = 2.8,
    label.size    = 0.15,
    label.r       = unit(0.10, "lines"),
    color         = "grey10",
    nudge_x       = poi_labels_df$nudge_x,
    nudge_y       = poi_labels_df$nudge_y,
    box.padding   = 0.55,
    point.padding = 0.30,
    min.segment.length = 0,
    max.overlaps  = 30,
    force         = 5,
    seed          = 42,
    show.legend   = FALSE,
    inherit.aes   = FALSE
  ) +
  
  # Lake labels (blue background, italic)
  {if (!is.null(top_lakes) && nrow(top_lakes) > 0) {
    lake_pts <- data.frame(
      lon  = sf::st_coordinates(top_lakes)[, 1],
      lat  = sf::st_coordinates(top_lakes)[, 2],
      name = top_lakes$GNIS_NAME_1,
      stringsAsFactors = FALSE
    )
    # Per-lake nudges:
    # Oosta Lake is already labelled by manual_labels_df → exclude here
    # to avoid a duplicate label on the same water body.
    lake_pts <- lake_pts[!lake_pts$name %in% "Oosta Lake", ]
    
    # Per-lake nudges:
    # • Takla Lake  — nudge_x=0.85 °: shorter connector (~60 km), label stays
    #                  inside basin, far from Kemano at lon -127.9
    # • Tahtsa Lake — nudge_x=0.65 °: clears the Kemano Powerhouse marker
    lake_nudge <- data.frame(
      name    = c("Takla Lake", "Tahtsa Lake"),
      nudge_x = c( 0.85,         0.65),
      nudge_y = c( 0.00,         0.10),
      stringsAsFactors = FALSE
    )
    lake_pts <- merge(lake_pts, lake_nudge, by="name", all.x=TRUE)
    lake_pts$nudge_x[is.na(lake_pts$nudge_x)] <- 0
    lake_pts$nudge_y[is.na(lake_pts$nudge_y)] <- 0
    
    ggrepel::geom_label_repel(
      data          = lake_pts,
      aes(x=lon, y=lat, label=name),
      fontface      = "italic",
      fill          = alpha(water_fill, 0.70),
      size          = 2.5,
      label.size    = 0.12,
      label.r       = unit(0.10, "lines"),
      color         = "grey10",
      nudge_x       = lake_pts$nudge_x,
      nudge_y       = lake_pts$nudge_y,
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
  
  coord_sf(crs=4326, expand=FALSE,
           xlim=map_xlim, ylim=map_ylim) +
  
  # Scale bar
  annotation_scale(
    location="br", width_hint=0.18,
    style="bar", text_cex=0.70,
    pad_x=unit(0.35,"cm"), pad_y=unit(0.35,"cm")
  ) +
  
  # Labels
  labs(
    title   = "Nechako Watershed",
    x       = "Longitude (\u00b0W)",
    y       = "Latitude (\u00b0N)",
    caption = paste0(
      "Elevation: AWS Terrain Tiles via elevatr (z=8, ~300 m)  |  ",
      "Rivers & Lakes: BC Freshwater Atlas via bcdata\n",
      "Fraser Watershed: BC FWA watershed groups  |  ",
      "Roads: BC Digital Road Atlas / NaturalEarth  |  ",
      "Basin: Nechako watershed (WGS84)\n",
      "\u25a0 Orange = Kenney Dam   \u25b2 Purple = Skins Lake Spillway   ",
      "\u25c6 Yellow = Kemano Powerhouse"
    )
  ) +
  
  # Theme
  theme_minimal(base_size=11, base_family="sans") +
  theme(
    plot.background  = element_rect(fill="white",   color=NA),
    panel.background = element_rect(fill=bg_land,   color=NA),
    panel.grid.major = element_line(color=alpha("white", 0.40),
                                    linetype="dashed", linewidth=0.20),
    panel.border     = element_rect(fill=NA, color="grey20", linewidth=0.7),
    plot.title       = element_text(face="bold", size=18, hjust=0.5,
                                    margin=margin(b=8)),
    plot.caption     = element_text(size=7.0, color="grey50", hjust=0,
                                    margin=margin(t=6)),
    axis.title       = element_text(size=9.0, color="grey30"),
    axis.text        = element_text(size=8.0, color="grey30"),
    legend.position  = "none"   # all legend items in legend_plot inset
  )

# ── 13. COMBINE MAIN MAP + INSETS ────────────────────────────────────────────────
cat("── Assembling final layout ──\n")

final_map <- main_map +
  # Locator inset — top-left (inside panel)
  inset_element(
    inset,
    left=0.00, bottom=0.68, right=0.19, top=0.985,
    align_to="plot"
  ) +
  # Legend inset — top-right (lon ≈ -122.5 to -121, outside Nechako basin)
  inset_element(
    legend_plot,
    left=0.835, bottom=0.28, right=0.995, top=0.985,
    align_to="panel"
  )

# ── 14. SAVE ─────────────────────────────────────────────────────────────────────
cat("── Saving outputs ──\n")

pdf_path <- file.path(out_dir, "nechako_overview_map_v2.pdf")
ggsave(pdf_path, plot=final_map, width=14, height=11, units="in",
       device=cairo_pdf)
cat(sprintf("  PDF: %s\n", pdf_path))

png_path <- file.path(out_dir, "nechako_overview_map_v2.png")
ggsave(png_path, plot=final_map, width=14, height=11, units="in",
       dpi=300, device="png", bg="white")
cat(sprintf("  PNG: %s\n", png_path))

cat("\n══ Done! ══\n")
cat("  nechako_overview_map_v2.pdf  — vector, print-quality\n")
cat("  nechako_overview_map_v2.png  — 300 dpi raster\n")