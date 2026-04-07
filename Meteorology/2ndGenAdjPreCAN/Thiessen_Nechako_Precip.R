##############################################################################
#  Thiessen_Nechako_Precip.R
#
#  PURPOSE
#    Compute Thiessen-polygon-weighted, basin-average monthly and annual
#    precipitation for the Nechako River Basin using AHCCD adjusted
#    station data (second-generation, updated to Dec 2017).
#
#  METHOD
#    1. All candidate AHCCD stations near/inside the basin are projected to
#       BC Albers Equal-Area (EPSG:3005) — an equal-area conic projection
#       with standard parallels at 50 N and 58.5 N, well-suited for BC.
#       Geographic (EPSG:4326) and UTM coordinates are AVOIDED because
#       geodesic area distortion over a basin this size (~27 000 km2,
#       spanning ~5 degrees of latitude) would produce systematically
#       incorrect Thiessen area weights.
#    2. Voronoi (Thiessen) polygons are built from the full station network
#       surrounding the basin, clipped to the basin boundary, and each
#       polygon's fraction of the basin area becomes that station's weight.
#       Stations with zero overlap are automatically excluded.
#    3. Monthly basin-average precipitation = sum(P_i * w_i) / sum(w_i),
#       re-normalised over stations with available data for that month.
#       A month is set to NA if less than 50 % of basin weight is covered.
#    4. Annual totals are summed from monthly values; years with fewer than
#       10 valid months are set to NA.
#
#  INPUTS
#    PREC_DIR   : folder containing the AHCCD mt*.txt files (from the zip)
#    BASIN_FILE : the doc.kml extracted from Nechako_Basin.kmz
#
#  OUTPUTS  (written to OUT_DIR)
#    Nechako_Thiessen_weights.csv          – station weights
#    Nechako_Thiessen_precip_monthly.csv   – monthly basin-average series
#    Nechako_Thiessen_precip_annual.csv    – annual basin-average series
#    Nechako_Thiessen_map.png              – Thiessen polygon map
#
#  AUTHOR   [your name]
#  DATE     [date]
##############################################################################


# ── 0.  Libraries ─────────────────────────────────────────────────────────────
required_pkgs <- c("sf", "dplyr", "tidyr", "ggplot2", "readr", "lubridate")
invisible(lapply(required_pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}))


# ── 1.  USER PATHS  (edit these) ──────────────────────────────────────────────
PREC_DIR   <- "prec"                      # folder with mt*.txt AHCCD files
BASIN_FILE <- "Nechako_Basin__1_.kmz"     # .kmz file — sf/GDAL reads it directly
OUT_DIR    <- "output"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)


# ── 2.  Projection ────────────────────────────────────────────────────────────
#  BC Albers Equal-Area Conic (EPSG:3005)
#    - Standard parallels: 50 N, 58.5 N  (bracket the Nechako basin)
#    - Equal-area: every polygon's computed area is physically correct
#    - This is the projection used by the BC government for all provincial
#      area calculations and is the right choice here.
CRS_BCALB <- 3005


# ── 3.  Station catalogue ─────────────────────────────────────────────────────
#  All AHCCD stations within ~200 km of the Nechako basin boundary.
#  Including surrounding stations is essential: they define the Voronoi
#  cell boundaries near the basin edge, preventing any single inside-basin
#  station from "claiming" area that should belong to a closer outside station.
#
#  Coordinates from Environment Canada station inventory (WGS84).

stations <- tibble::tribble(
  ~id,          ~name,                         ~lat,      ~lon,
  # ── Inside basin ──────────────────────────────────────────────────────
  "1085835",   "OOTSA L SKINS L SPILLWAY",     53.733,  -125.550,
  "1092970",   "FORT ST JAMES",                54.450,  -124.250,
  # ── Near outside (<150 km) ────────────────────────────────────────────
  "1096450",   "PRINCE GEORGE",                53.883,  -122.683,   # ~4 km
  "1064020",   "KEMANO",                       53.567,  -127.983,   # ~20 km
  "1096630",   "QUESNEL",                      53.033,  -122.517,   # ~59 km
  "1060841",   "BELLA COOLA",                  52.367,  -126.767,   # ~62 km
  "1077500",   "SMITHERS",                     54.817,  -127.183,   # ~69 km
  "1183090",   "GERMANSEN LANDING",            55.783,  -124.717,   # ~73 km
  "1073347",   "HAZELTON TEMLEHAN",            55.267,  -127.667,   # ~99 km
  "1068130",   "TERRACE",                      54.467,  -128.567,   # ~113 km
  "1080870",   "BIG CREEK",                    51.967,  -123.033,   # ~134 km
  "1090660",   "BARKERVILLE",                  53.067,  -121.517,   # ~140 km
  "1088010",   "TATLAYOKO LAKE",               51.667,  -124.400,   # ~147 km
  "1098940",   "WILLIAMS LAKE",                52.167,  -122.050,   # ~155 km
  "1183000",   "FORT ST JOHN",                 56.233,  -120.850    # ~289 km (north anchor)
)


# ── 4.  Read basin boundary ─────────────────────────────────────────────────
message("Reading basin boundary ...")

# ── Resolve absolute path first (avoids working-directory pitfalls) ────────
abs_basin <- normalizePath(BASIN_FILE, mustWork = FALSE)
if (!file.exists(abs_basin))
  stop("Basin file not found: ", abs_basin,
       "\n  Check BASIN_FILE at the top of this script.")

# ── Build GDAL-readable path via /vsizip/ ──────────────────────────────────
# R's unzip() requires the working directory to contain the file, which is
# fragile.  GDAL's /vsizip/ virtual filesystem reads KMZ/ZIP archives
# directly from their absolute path without any extraction:
#
#   /vsizip//absolute/path/to/file.kmz/doc.kml
#
# All standard KMZ files (Google Earth, ArcGIS) put the KML at "doc.kml"
# inside the archive root.  We try that first, then fall back to probing.

kml_path <- if (grepl("\\.kmz$", abs_basin, ignore.case = TRUE)) {
  
  candidate <- paste0("/vsizip/", abs_basin, "/doc.kml")
  
  gdal_ok <- tryCatch({ sf::st_layers(candidate); TRUE }, error = function(e) FALSE)
  
  if (gdal_ok) {
    message("  Using GDAL vsizip path: .../doc.kml")
    candidate
    
  } else {
    # Non-standard internal name: probe the archive root as a GDAL directory
    root_path   <- paste0("/vsizip/", abs_basin)
    root_layers <- tryCatch(sf::st_layers(root_path)$name,
                            error = function(e) character(0))
    kml_hits <- grep("\\.kml$", root_layers, value = TRUE, ignore.case = TRUE)
    
    if (!length(kml_hits))
      stop("No .kml found inside: ", abs_basin,
           "\n  Archive contents seen by GDAL: ",
           if (length(root_layers)) paste(root_layers, collapse = ", ") else "(none — check file integrity)")
    
    chosen <- paste0(root_path, "/", kml_hits[1])
    message(sprintf("  Non-standard KML entry: %s", kml_hits[1]))
    chosen
  }
  
} else {
  message("  Reading KML directly (no vsizip needed).")
  abs_basin
}

# ── List layers and read ───────────────────────────────────────────────────
kmz_layers <- sf::st_layers(kml_path)$name
message(sprintf("  KML layers: %s", paste(kmz_layers, collapse = ", ")))

basin_raw <- st_read(kml_path, layer = kmz_layers[1], quiet = TRUE)

poly_types <- c("POLYGON", "MULTIPOLYGON", "GEOMETRYCOLLECTION")
basin_raw  <- basin_raw[st_geometry_type(basin_raw) %in% poly_types, ]

if (nrow(basin_raw) == 0)
  stop("No polygon features in layer '", kmz_layers[1], "'.",
       "\n  Available layers: ", paste(kmz_layers, collapse = ", "))

basin_sf <- basin_raw  |>
  st_union()            |>
  st_make_valid()       |>
  st_transform(CRS_BCALB)

basin_area_km2 <- as.numeric(st_area(basin_sf)) / 1e6
message(sprintf("  Basin area (BC Albers): %.0f km2", basin_area_km2))


# ── 5.  Project stations to BC Albers ─────────────────────────────────────────
stations_sf <- stations |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) |>
  st_transform(CRS_BCALB)


# ── 6.  Build Thiessen (Voronoi) polygons ─────────────────────────────────────
message("Building Voronoi polygons ...")

# Envelope: rectangular bbox 500 km around the basin in projected coords.
# A large envelope prevents any border-station Voronoi cell from being
# prematurely clipped before it intersects the basin.
envelope_sfc <- st_as_sfc(
  st_bbox(st_buffer(basin_sf, 500000))
)

# Compute Voronoi tessellation on the merged (multi)point geometry
vor_geom <- st_voronoi(
  st_union(stations_sf),
  envelope = envelope_sfc
)

# Extract individual polygons from the geometry collection
vor_polys <- st_collection_extract(vor_geom, "POLYGON") |> st_sf()

# --- Match each Voronoi polygon to the station it contains -----------------
#  st_voronoi does NOT preserve point order, so we must do a spatial join.
vor_polys$station_id <- NA_character_

for (i in seq_len(nrow(vor_polys))) {
  # which station point falls inside this polygon?
  inside <- as.logical(
    st_within(stations_sf, vor_polys[i, ], sparse = FALSE)[, 1]
  )
  if (any(inside)) vor_polys$station_id[i] <- stations$id[which(inside)]
}

# Drop any unmatched polygons (shouldn't occur, but defensive)
vor_polys <- filter(vor_polys, !is.na(station_id))


# ── 7.  Clip Voronoi polygons to basin and compute area weights ───────────────
message("Clipping to basin and computing area weights ...")

# Suppress sf warnings about attribute-geometry relationships during intersection
suppressWarnings(
  vor_clipped <- st_intersection(vor_polys, basin_sf)
)

# Area in m2 (BC Albers is equal-area, so st_area is exact)
# st_area() on the sf object uses the active geometry column regardless of its name
vor_clipped$area_m2 <- as.numeric(st_area(vor_clipped))

# Some intersection slivers may be near-zero; clean them up
vor_clipped <- filter(vor_clipped, area_m2 > 1e4)  # > 0.01 km2

total_area <- sum(vor_clipped$area_m2)
vor_clipped$weight <- vor_clipped$area_m2 / total_area

# Join station metadata
weights_tbl <- vor_clipped |>
  st_drop_geometry() |>
  select(station_id, area_m2, weight) |>
  left_join(stations |> select(id, name, lat, lon),
            by = c("station_id" = "id")) |>
  arrange(desc(weight)) |>
  mutate(
    area_km2     = round(area_m2 / 1e6, 1),
    weight_pct   = round(weight * 100, 2)
  )

# Report
message("\n── Thiessen Weights (influential stations only) ──────────────────────")
print(
  weights_tbl |>
    select(station_id, name, lat, lon, area_km2, weight_pct) |>
    rename(`Weight (%)` = weight_pct, `Area (km2)` = area_km2)
)

influential_ids <- weights_tbl$station_id
weights_named   <- setNames(weights_tbl$weight, weights_tbl$station_id)


# ── 8.  Parse AHCCD monthly total-precipitation files ─────────────────────────
#
#  File format (comma-delimited, encoding latin1):
#    Line 1  : station metadata (English)
#    Line 2  : station metadata (French)
#    Line 3  : column headers (English)
#    Line 4  : column headers (French)
#    Lines 5+: data rows
#
#  Data row structure — 36 columns:
#    Col  1       : Year (integer)
#    Cols 2,3     : January  value (mm), flag
#    Cols 4,5     : February value (mm), flag
#    ...
#    Cols 24,25   : December value (mm), flag
#    Cols 26,27   : Annual   value (mm), flag
#    Cols 28,29   : Winter   value (mm), flag
#    Cols 30,31   : Spring   value (mm), flag
#    Cols 32,33   : Summer   value (mm), flag
#    Cols 34,35   : Autumn   value (mm), flag
#    Col  36      : trailing empty (artefact of trailing comma)
#
#  Missing values : -9999.9 (flag = "M")

MONTHS     <- c("Jan","Feb","Mar","Apr","May","Jun",
                "Jul","Aug","Sep","Oct","Nov","Dec")
# Column indices of the 12 monthly VALUES in each data row (1-based)
MON_COLS   <- seq(2, 24, by = 2)   # 2,4,6,...,24

parse_ahccd <- function(filepath, station_id) {
  
  lines <- tryCatch(
    readLines(filepath, encoding = "latin1", warn = FALSE),
    error = function(e) { warning("Cannot read: ", filepath); return(NULL) }
  )
  if (is.null(lines)) return(NULL)
  
  # Data starts at line 5; drop blank lines
  data_lines <- lines[5:length(lines)]
  data_lines <- data_lines[nzchar(trimws(data_lines))]
  if (!length(data_lines)) return(NULL)
  
  # Read all 36 columns as character to avoid type-conversion issues,
  # then coerce numeric columns explicitly.
  raw <- tryCatch(
    read.csv(
      text            = paste(data_lines, collapse = "\n"),
      header          = FALSE,
      colClasses      = "character",
      strip.white     = TRUE,
      na.strings      = c("", " ")
    ),
    error = function(e) {
      warning("Parse error in ", basename(filepath), ": ", e$message)
      return(NULL)
    }
  )
  if (is.null(raw)) return(NULL)
  
  # Year column
  year_vec <- suppressWarnings(as.integer(raw[[1]]))
  
  # Extract and coerce the 12 monthly value columns
  monthly_vals <- raw[, MON_COLS, drop = FALSE]
  monthly_vals <- lapply(monthly_vals, function(x) {
    v <- suppressWarnings(as.numeric(x))
    v[v <= -9999] <- NA   # replace -9999.9 sentinel
    v
  })
  names(monthly_vals) <- MONTHS
  
  result <- bind_cols(
    tibble(Year = year_vec),
    as_tibble(monthly_vals)
  ) |>
    filter(!is.na(Year)) |>   # drop any stray header rows
    pivot_longer(
      cols      = all_of(MONTHS),
      names_to  = "Month_name",
      values_to = "Precip_mm"
    ) |>
    mutate(
      Month      = match(Month_name, MONTHS),
      Date       = as.Date(sprintf("%04d-%02d-01", Year, Month)),
      station_id = station_id
    ) |>
    select(station_id, Date, Year, Month, Precip_mm)
  
  return(result)
}

# Load data for all influential stations
message("\nReading AHCCD station files ...")
precip_list <- lapply(influential_ids, function(sid) {
  fpath <- file.path(PREC_DIR, paste0("mt", sid, ".txt"))
  if (!file.exists(fpath)) {
    warning("  File not found: ", fpath)
    return(NULL)
  }
  df <- parse_ahccd(fpath, sid)
  if (!is.null(df)) {
    yr_range <- range(df$Year, na.rm = TRUE)
    n_obs    <- sum(!is.na(df$Precip_mm))
    message(sprintf("  %-10s  %-32s  %d–%d  (%d valid months)",
                    sid,
                    weights_tbl$name[weights_tbl$station_id == sid],
                    yr_range[1], yr_range[2], n_obs))
  }
  df
})
names(precip_list) <- influential_ids
precip_all <- bind_rows(precip_list)


# ── 9.  Compute weighted monthly basin average ────────────────────────────────
message("\nComputing weighted monthly series ...")

monthly_avg <- precip_all |>
  left_join(
    tibble(station_id = names(weights_named), weight = weights_named),
    by = "station_id"
  ) |>
  group_by(Date, Year, Month) |>
  summarise(
    # Sum of weights for stations that have data this month
    w_avail    = sum(weight[!is.na(Precip_mm)], na.rm = TRUE),
    # Weighted mean, re-normalised to the available weight fraction
    Precip_mm  = if (all(is.na(Precip_mm))) {
      NA_real_
    } else {
      sum(Precip_mm * weight, na.rm = TRUE) /
        sum(weight[!is.na(Precip_mm)], na.rm = TRUE)
    },
    N_stations = as.integer(sum(!is.na(Precip_mm))),
    .groups    = "drop"
  ) |>
  mutate(
    # Fraction of basin area covered by non-missing stations (0–1)
    Coverage_frac = round(w_avail, 4),
    # Flag months where data coverage < 50 % of basin area
    Precip_mm     = ifelse(Coverage_frac < 0.50, NA_real_, Precip_mm),
    Precip_mm     = round(Precip_mm, 2)
  ) |>
  arrange(Date) |>
  select(Date, Year, Month, Precip_mm, N_stations, Coverage_frac)


# ── 10.  Compute annual basin-average totals ───────────────────────────────────
#   Annual total = sum of 12 monthly values.
#   Rule: at least 10 valid months required; otherwise NA.
#   When 10–11 months are valid, the total is prorated to 12 months
#   (multiply by 12 / N_valid) with a flag noting the adjustment.

message("Computing annual series ...")

annual_avg <- monthly_avg |>
  group_by(Year) |>
  summarise(
    N_months_valid = sum(!is.na(Precip_mm)),
    Precip_sum_obs = sum(Precip_mm, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    # Prorated estimate when 10 or 11 months available
    Precip_mm_annual = case_when(
      N_months_valid == 12 ~ round(Precip_sum_obs, 1),
      N_months_valid >= 10 ~ round(Precip_sum_obs * 12 / N_months_valid, 1),
      TRUE                 ~ NA_real_
    ),
    Record_complete  = N_months_valid == 12,
    Note = case_when(
      N_months_valid == 12 ~ "complete",
      N_months_valid >= 10 ~ sprintf("prorated from %d months", N_months_valid),
      TRUE                 ~ "insufficient data (< 10 months)"
    )
  ) |>
  select(Year, Precip_mm_annual, N_months_valid, Record_complete, Note)


# ── 11.  Write CSV outputs ─────────────────────────────────────────────────────
message("\nWriting output files ...")

# (a) Thiessen weights
write_csv(
  weights_tbl |>
    select(station_id, name, lat, lon, area_km2, weight, weight_pct),
  file.path(OUT_DIR, "Nechako_Thiessen_weights.csv"),
  na = "NA"
)

# (b) Monthly series
write_csv(
  monthly_avg,
  file.path(OUT_DIR, "Nechako_Thiessen_precip_monthly.csv"),
  na = "NA"
)

# (c) Annual series
write_csv(
  annual_avg,
  file.path(OUT_DIR, "Nechako_Thiessen_precip_annual.csv"),
  na = "NA"
)

message(sprintf("  → %s", file.path(OUT_DIR, "Nechako_Thiessen_weights.csv")))
message(sprintf("  → %s", file.path(OUT_DIR, "Nechako_Thiessen_precip_monthly.csv")))
message(sprintf("  → %s", file.path(OUT_DIR, "Nechako_Thiessen_precip_annual.csv")))


# ── 12.  Diagnostic map ───────────────────────────────────────────────────────
message("Generating Thiessen map ...")

# Attach station names and weights to clipped polygons for plotting
vor_plot <- vor_clipped |>
  left_join(
    weights_tbl |> select(station_id, name, weight_pct),
    by = "station_id"
  )

# Clip stations to a plotting extent
basin_bbox   <- st_bbox(st_buffer(basin_sf, 80000))

map_plot <- ggplot() +
  # Thiessen polygons (fill by weight)
  geom_sf(data    = vor_plot,
          aes(fill = weight_pct),
          colour  = "white", linewidth = 0.6, alpha = 0.85) +
  scale_fill_distiller(
    name   = "Basin\nweight (%)",
    palette = "YlOrRd",
    direction = 1
  ) +
  # Basin boundary
  geom_sf(data      = basin_sf,
          fill      = NA,
          colour    = "black",
          linewidth = 1.0) +
  # All candidate stations (open circles)
  geom_sf(data   = stations_sf,
          shape  = 21, fill = "grey90", colour = "black",
          size   = 2.2, stroke = 0.8) +
  # Influential stations (filled circles, larger)
  geom_sf(data   = filter(stations_sf, id %in% influential_ids),
          shape  = 21, fill = "steelblue", colour = "black",
          size   = 3.5, stroke = 1.0) +
  # Labels for influential stations
  geom_sf_text(
    data  = filter(stations_sf, id %in% influential_ids),
    aes(label = name),
    size  = 2.8, fontface = "bold",
    nudge_y = 8000, check_overlap = TRUE
  ) +
  coord_sf(
    xlim = c(basin_bbox["xmin"], basin_bbox["xmax"]),
    ylim = c(basin_bbox["ymin"], basin_bbox["ymax"])
  ) +
  labs(
    title    = "Nechako Basin — Thiessen Polygon Weights",
    subtitle = "BC Albers Equal-Area (EPSG:3005) | Open circles = candidate stations",
    x        = NULL, y = NULL,
    caption  = paste0("Basin area: ", round(basin_area_km2, 0), " km²")
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position   = "right",
    plot.title        = element_text(face = "bold", size = 13),
    plot.subtitle     = element_text(size = 9, colour = "grey40"),
    panel.grid.major  = element_line(colour = "grey85", linewidth = 0.3)
  )

ggsave(
  filename = file.path(OUT_DIR, "Nechako_Thiessen_map.png"),
  plot     = map_plot,
  width    = 10, height = 7, dpi = 150
)
message(sprintf("  → %s", file.path(OUT_DIR, "Nechako_Thiessen_map.png")))


# ── 13.  Quick sanity summary ─────────────────────────────────────────────────
message("\n── Summary ───────────────────────────────────────────────────────────")
message(sprintf("  Monthly series: %d months  (%d–%d)",
                nrow(monthly_avg),
                min(monthly_avg$Year), max(monthly_avg$Year)))
message(sprintf("  Valid monthly values: %d / %d  (%.1f%%)",
                sum(!is.na(monthly_avg$Precip_mm)),
                nrow(monthly_avg),
                100 * mean(!is.na(monthly_avg$Precip_mm))))
message(sprintf("  Annual series:  %d years  (%d complete)",
                nrow(annual_avg),
                sum(annual_avg$Record_complete, na.rm = TRUE)))
message(sprintf("  Mean annual precip (complete years): %.0f mm",
                mean(annual_avg$Precip_mm_annual[annual_avg$Record_complete],
                     na.rm = TRUE)))
message("\nDone.\n")