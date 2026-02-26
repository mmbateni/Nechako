####################################################################################
# NECHAKO BASIN — ELEVATION MAP  (v4 — terra::plot native rendering)
#
# Root cause of all previous white maps:
#   v1/v2: rasterImage() silently ignores a character hex matrix -> white
#   v3:    plotRGB() with integer 0-255 bands and scale=255 was double-scaling
#
# THIS VERSION:
#   - terra::plot(dem, col=elev_cols)   <- terra handles orientation internally
#   - hillshade layered with semi-transparent colours via plot(add=TRUE)
#   - Manual legend bar drawn after
#   - Zero custom pixel manipulation
####################################################################################

library(terra)
library(sf)

# ===== PATHS =====
dem_dir      <- "D:/Nechako_Drought/Nechako/Spatial/nechako-dem"
basin_shp    <- "D:/Nechako_Drought/Nechako/Spatial/nechakoBound_dissolve.shp"
stations_csv <- "D:/Nechako_Drought/Nechako/Hydrology/data_retrievalWaterSurveyofCanada/stations.csv"  # <-- UPDATE THIS PATH

out_dir      <- "D:/Nechako_Drought/Nechako/elevation"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

bg_col <- "#d6eaf8"

# ===== LOAD DEM =====
cat("Loading DEM...\n")
dem <- rast(dem_dir)
cat(sprintf("  Extent    : %s\n", as.character(ext(dem))))
cat(sprintf("  Resolution: %.6f deg\n", res(dem)[1]))

# ===== LOAD BASIN BOUNDARY =====
cat("Loading basin boundary...\n")
basin_sf <- st_read(basin_shp, quiet = TRUE)
basin_v  <- vect(basin_shp)
if (!same.crs(basin_v, crs(dem))) {
  basin_v  <- project(basin_v, crs(dem))
  basin_sf <- st_transform(basin_sf, st_crs(crs(dem)))
}

# Re-mask in R so the NA edge matches the vector boundary exactly
cat("Masking DEM to basin boundary...\n")
dem <- mask(dem, basin_v, touches = TRUE)

elev_range <- as.numeric(global(dem, c("min","max"), na.rm = TRUE))
cat(sprintf("  Elevation : %.0f - %.0f m\n", elev_range[1], elev_range[2]))

# Sanity check: print a sample of actual pixel values
samp <- spatSample(dem, 10, na.rm = TRUE)
cat(sprintf("  Sample values (m): %s\n",
            paste(round(samp[[1]]), collapse = ", ")))

# ===== HILLSHADE =====
cat("Computing hillshade...\n")
dem_m    <- project(dem, "EPSG:3005", method = "bilinear")
slope_r  <- terrain(dem_m, "slope",  unit = "radians")
aspect_r <- terrain(dem_m, "aspect", unit = "radians")
hill     <- shade(slope_r, aspect_r, angle = 40, direction = 315, normalize = TRUE)
hill     <- project(hill, crs(dem), method = "bilinear")
hill     <- resample(hill, dem, method = "bilinear")

# ===== COLOUR PALETTES =====
n_col <- 512

# Hypsometric: green lowlands -> brown mid -> white peaks
elev_cols <- c(
  colorRampPalette(c("#336600","#669933","#99cc66"))(round(n_col * 0.20)),
  colorRampPalette(c("#cccc66","#cc9933","#996633"))(round(n_col * 0.35)),
  colorRampPalette(c("#996633","#cc9966","#d4b896"))(round(n_col * 0.25)),
  colorRampPalette(c("#d4b896","#e8dcc8","#f5f0e8","#ffffff"))(round(n_col * 0.20))
)

# Hillshade overlay colours: shadows are semi-transparent dark grey,
# highlights are semi-transparent white, neutral (0.5) is fully transparent
n_hs <- 256
hill_overlay_cols <- sapply(seq(0, 1, length.out = n_hs), function(h) {
  if (h < 0.5) rgb(0, 0, 0, alpha = (0.5 - h) * 0.9)   # shadow
  else         rgb(1, 1, 1, alpha = (h - 0.5) * 0.6)   # highlight
})

# ===== MAIN ELEVATION MAP PDF =====
pdf_path <- file.path(out_dir, "nechako_elevation_map.pdf")
cat(sprintf("Rendering -> %s\n", basename(pdf_path)))

pdf(pdf_path, width = 11, height = 9, family = "Helvetica")

par(mar = c(5.5, 4.5, 3.5, 5.5),
    oma = c(0, 0, 0, 0),
    bg  = bg_col)

# Layer 1 — elevation (terra handles row/col orientation correctly)
plot(dem,
     col    = elev_cols,
     legend = FALSE,
     axes   = TRUE,
     box    = FALSE,
     xlab   = "Longitude (W)",
     ylab   = "Latitude (N)",
     colNA  = bg_col,
     main   = "")

# Layer 2 — semi-transparent hillshade overlay
plot(hill,
     col    = hill_overlay_cols,
     legend = FALSE,
     axes   = FALSE,
     box    = FALSE,
     colNA  = NA,
     add    = TRUE)

# Layer 3 — basin boundary (border removed)
# plot(st_geometry(basin_sf), add = TRUE,
#      col = NA, border = "black", lwd = 2.2)

# Grid lines
grid(nx = NULL, ny = NULL, col = "white", lty = 2, lwd = 0.35)

# ===== MANUAL ELEVATION LEGEND =====
legend_x  <- grconvertX(0.930, "npc", "user")
legend_x2 <- grconvertX(0.965, "npc", "user")
legend_y0 <- grconvertY(0.10,  "npc", "user")
legend_y1 <- grconvertY(0.90,  "npc", "user")
n_leg     <- 150
leg_y     <- seq(legend_y0, legend_y1, length.out = n_leg + 1)
leg_cols  <- colorRampPalette(elev_cols)(n_leg)
for (i in seq_len(n_leg))
  rect(legend_x, leg_y[i], legend_x2, leg_y[i+1],
       col = leg_cols[i], border = NA)
rect(legend_x, legend_y0, legend_x2, legend_y1, border = "black", lwd = 0.8)

leg_elevs <- pretty(elev_range, n = 6)
leg_elevs <- leg_elevs[leg_elevs >= elev_range[1] & leg_elevs <= elev_range[2]]
leg_ypos  <- legend_y0 +
  (leg_elevs - elev_range[1]) / diff(elev_range) *
  (legend_y1 - legend_y0)
tick_len  <- (legend_x2 - legend_x) * 0.6
for (i in seq_along(leg_elevs)) {
  segments(legend_x2, leg_ypos[i], legend_x2 + tick_len, leg_ypos[i], lwd = 0.7)
  text(legend_x2 + tick_len * 1.2, leg_ypos[i],
       paste0(leg_elevs[i], " m"), adj = 0, cex = 0.78)
}
text((legend_x + legend_x2) / 2,
     legend_y1 + (legend_y1 - legend_y0) * 0.05,
     "Elevation", adj = c(0.5, 0), cex = 0.9, font = 2)

# ===== TITLE, FOOTNOTE, NORTH ARROW =====
title(main = "Nechako Basin - Digital Elevation Model",
      line = 1.5, cex.main = 1.4, font.main = 2)
mtext(sprintf(
  "WGS84 | Resolution: %.4f deg | Elevation: %.0f - %.0f m | Source: n50w130_dem",
  res(dem)[1], elev_range[1], elev_range[2]),
  side = 1, line = 3, cex = 0.8, col = "gray30")

dx <- diff(c(xmin(dem), xmax(dem)))
dy <- diff(c(ymin(dem), ymax(dem)))
arrows(xmin(dem) + dx * 0.06, ymin(dem) + dy * 0.05,
       xmin(dem) + dx * 0.06, ymin(dem) + dy * 0.10,
       length = 0.12, lwd = 2, col = "black")
text(xmin(dem) + dx * 0.06, ymin(dem) + dy * 0.115,
     "N", cex = 1.0, font = 2)

dev.off()
cat(sprintf("Elevation map saved: %s\n", pdf_path))

# ===== SLOPE AND ASPECT (unchanged — already working) =====
cat("Rendering slope + aspect map...\n")
slope_deg  <- terrain(dem_m, "slope",  unit = "degrees")
aspect_deg <- terrain(dem_m, "aspect", unit = "degrees")
slope_deg  <- project(slope_deg,  crs(dem))
aspect_deg <- project(aspect_deg, crs(dem))

pdf(file.path(out_dir, "nechako_slope_aspect.pdf"),
    width = 12, height = 6, family = "Helvetica")
par(mfrow = c(1, 2), mar = c(3, 3, 2.5, 1), oma = c(0, 0, 2, 0))

plot(slope_deg, main = "Slope (degrees)",
     col = hcl.colors(101, "YlOrRd", rev = FALSE),
     axes = TRUE, box = FALSE,
     plg = list(title = "deg", cex = 0.85), cex.main = 1.1)
# plot(st_geometry(basin_sf), add = TRUE, col = NA, border = "black", lwd = 2)

plot(aspect_deg, main = "Aspect (degrees from North)",
     col = hcl.colors(101, "RdYlBu", rev = FALSE),
     axes = TRUE, box = FALSE,
     plg = list(title = "deg", cex = 0.85), cex.main = 1.1)
# plot(st_geometry(basin_sf), add = TRUE, col = NA, border = "black", lwd = 2)

mtext("Nechako Basin - Terrain Derivatives",
      outer = TRUE, cex = 1.2, font = 2, line = 0.3)
dev.off()
cat(sprintf("Slope + aspect saved: %s\n",
            file.path(out_dir, "nechako_slope_aspect.pdf")))

####################################################################################
# MAP 3 — ELEVATION + HYDROMETRY STATIONS
####################################################################################
cat("Rendering elevation + hydrometry stations map...\n")

library(readr)

stations_df <- read_csv(stations_csv)

# Keep only stations with valid coordinates
stations_df <- stations_df[!is.na(stations_df$LATITUDE) & !is.na(stations_df$LONGITUDE), ]

pdf(file.path(out_dir, "nechako_elevation_hydrostations.pdf"),
    width = 11, height = 9, family = "Helvetica")

# Explicitly reset multi-panel layout that may carry over from slope/aspect map
par(mfrow = c(1, 1),
    mar   = c(5.5, 4.5, 3.5, 5.5),
    oma   = c(0, 0, 0, 0),
    bg    = bg_col)

# Anchor plot extent explicitly from DEM extent
dem_ext <- as.vector(ext(dem))   # xmin, xmax, ymin, ymax
xlims   <- dem_ext[1:2]
ylims   <- dem_ext[3:4]

# Layer 1 — elevation
plot(dem,
     col    = elev_cols,
     legend = FALSE,
     axes   = TRUE,
     box    = FALSE,
     xlim   = xlims,
     ylim   = ylims,
     xlab   = "Longitude (W)",
     ylab   = "Latitude (N)",
     colNA  = bg_col,
     main   = "")

# Layer 2 — hillshade
plot(hill,
     col    = hill_overlay_cols,
     legend = FALSE,
     axes   = FALSE,
     box    = FALSE,
     colNA  = NA,
     add    = TRUE)

# Grid
grid(nx = NULL, ny = NULL, col = "white", lty = 2, lwd = 0.35)

# Layer 3 — hydrometry station points
points(stations_df$LONGITUDE, stations_df$LATITUDE,
       pch = 21, bg = "yellow", col = "black", cex = 1.0, lwd = 0.8)

# Labels — single pass, check.overlap suppresses colliding labels automatically
text(stations_df$LONGITUDE, stations_df$LATITUDE,
     labels = stations_df$STATION_NUMBER,
     pos = 3, cex = 0.50, col = "black", font = 1,
     offset = 0.35, check.overlap = TRUE)

# ===== LEGEND BAR =====
legend_x  <- grconvertX(0.930, "npc", "user")
legend_x2 <- grconvertX(0.965, "npc", "user")
legend_y0 <- grconvertY(0.10,  "npc", "user")
legend_y1 <- grconvertY(0.90,  "npc", "user")
n_leg     <- 150
leg_y     <- seq(legend_y0, legend_y1, length.out = n_leg + 1)
leg_cols  <- colorRampPalette(elev_cols)(n_leg)
for (i in seq_len(n_leg))
  rect(legend_x, leg_y[i], legend_x2, leg_y[i+1],
       col = leg_cols[i], border = NA)
rect(legend_x, legend_y0, legend_x2, legend_y1, border = "black", lwd = 0.8)

leg_elevs <- pretty(elev_range, n = 6)
leg_elevs <- leg_elevs[leg_elevs >= elev_range[1] & leg_elevs <= elev_range[2]]
leg_ypos  <- legend_y0 +
  (leg_elevs - elev_range[1]) / diff(elev_range) * (legend_y1 - legend_y0)
tick_len  <- (legend_x2 - legend_x) * 0.6
for (i in seq_along(leg_elevs)) {
  segments(legend_x2, leg_ypos[i], legend_x2 + tick_len, leg_ypos[i], lwd = 0.7)
  text(legend_x2 + tick_len * 1.2, leg_ypos[i],
       paste0(leg_elevs[i], " m"), adj = 0, cex = 0.78)
}
text((legend_x + legend_x2) / 2,
     legend_y1 + (legend_y1 - legend_y0) * 0.05,
     "Elevation", adj = c(0.5, 0), cex = 0.9, font = 2)

# Station legend
legend("bottomleft",
       legend = "Hydrometry Station",
       pch = 21, pt.bg = "yellow", col = "black",
       pt.cex = 1.1, cex = 0.8, bg = "white", box.lwd = 0.5,
       inset = c(0.01, 0.02))

title(main = "Nechako Basin - Hydrometry Stations",
      line = 1.5, cex.main = 1.4, font.main = 2)
mtext(sprintf(
  "WGS84 | Resolution: %.4f deg | Elevation: %.0f - %.0f m | Source: n50w130_dem",
  res(dem)[1], elev_range[1], elev_range[2]),
  side = 1, line = 3, cex = 0.8, col = "gray30")

# North arrow
dx <- diff(c(xmin(dem), xmax(dem)))
dy <- diff(c(ymin(dem), ymax(dem)))
arrows(xmin(dem) + dx * 0.06, ymin(dem) + dy * 0.05,
       xmin(dem) + dx * 0.06, ymin(dem) + dy * 0.10,
       length = 0.12, lwd = 2, col = "black")
text(xmin(dem) + dx * 0.06, ymin(dem) + dy * 0.115,
     "N", cex = 1.0, font = 2)

dev.off()
cat(sprintf("Hydrometry stations map saved: %s\n",
            file.path(out_dir, "nechako_elevation_hydrostations.pdf")))

####################################################################################
# MAP 4 — ELEVATION + SNOW PILLOWS (real observed locations)
####################################################################################
cat("Rendering elevation + snow pillows map...\n")

# Real (observed) coordinates from station metadata
snow_pillows <- data.frame(
  id   = c("1B08P",      "1B01P",       "1B02P"),
  name = c("Mt. Pondosy", "Mount Wells", "Tahtsa Lake"),
  lat  = c(53.15,         53.72,         53.59),
  lon  = c(-126.88,       -126.42,       -127.64)
)

pdf(file.path(out_dir, "nechako_elevation_snowpillows.pdf"),
    width = 11, height = 9, family = "Helvetica")

par(mfrow = c(1, 1),
    mar   = c(5.5, 4.5, 3.5, 5.5),
    oma   = c(0, 0, 0, 0),
    bg    = bg_col)

# Layer 1 — elevation
plot(dem,
     col    = elev_cols,
     legend = FALSE,
     axes   = TRUE,
     box    = FALSE,
     xlim   = xlims,
     ylim   = ylims,
     xlab   = "Longitude (W)",
     ylab   = "Latitude (N)",
     colNA  = bg_col,
     main   = "")

# Layer 2 — hillshade
plot(hill,
     col    = hill_overlay_cols,
     legend = FALSE,
     axes   = FALSE,
     box    = FALSE,
     colNA  = NA,
     add    = TRUE)

# Grid
grid(nx = NULL, ny = NULL, col = "white", lty = 2, lwd = 0.35)

# Layer 3 — snow pillow points
points(snow_pillows$lon, snow_pillows$lat,
       pch = 23, bg = "deepskyblue", col = "black", cex = 1.6, lwd = 1.2)

# Labels: per-point pos to avoid clipping at map edges
# pos: 1=below, 2=left, 3=above, 4=right
# Tahtsa Lake (index 3) is the westernmost point — label to the right (pos=4)
id_pos   <- c(3, 3, 4)   # Mt.Pondosy above, Mt.Wells above, Tahtsa Lake right
name_pos <- c(1, 1, 4)   # names below for first two, right for Tahtsa Lake

# ID labels
for (i in seq_len(nrow(snow_pillows))) {
  text(snow_pillows$lon[i], snow_pillows$lat[i],
       labels = snow_pillows$id[i],
       pos = id_pos[i], cex = 0.85, col = "black", font = 2, offset = 0.5)
}
# Station name labels (slightly smaller offset below ID for non-Tahtsa stations)
for (i in seq_len(nrow(snow_pillows))) {
  text(snow_pillows$lon[i], snow_pillows$lat[i],
       labels = snow_pillows$name[i],
       pos = name_pos[i], cex = 0.75, col = "gray20", font = 3,
       offset = if (name_pos[i] == 4) 3.5 else 0.5)
}

# ===== LEGEND BAR =====
legend_x  <- grconvertX(0.930, "npc", "user")
legend_x2 <- grconvertX(0.965, "npc", "user")
legend_y0 <- grconvertY(0.10,  "npc", "user")
legend_y1 <- grconvertY(0.90,  "npc", "user")
n_leg     <- 150
leg_y     <- seq(legend_y0, legend_y1, length.out = n_leg + 1)
leg_cols  <- colorRampPalette(elev_cols)(n_leg)
for (i in seq_len(n_leg))
  rect(legend_x, leg_y[i], legend_x2, leg_y[i+1],
       col = leg_cols[i], border = NA)
rect(legend_x, legend_y0, legend_x2, legend_y1, border = "black", lwd = 0.8)

leg_elevs <- pretty(elev_range, n = 6)
leg_elevs <- leg_elevs[leg_elevs >= elev_range[1] & leg_elevs <= elev_range[2]]
leg_ypos  <- legend_y0 +
  (leg_elevs - elev_range[1]) / diff(elev_range) * (legend_y1 - legend_y0)
tick_len  <- (legend_x2 - legend_x) * 0.6
for (i in seq_along(leg_elevs)) {
  segments(legend_x2, leg_ypos[i], legend_x2 + tick_len, leg_ypos[i], lwd = 0.7)
  text(legend_x2 + tick_len * 1.2, leg_ypos[i],
       paste0(leg_elevs[i], " m"), adj = 0, cex = 0.78)
}
text((legend_x + legend_x2) / 2,
     legend_y1 + (legend_y1 - legend_y0) * 0.05,
     "Elevation", adj = c(0.5, 0), cex = 0.9, font = 2)

# Snow pillow legend
legend("bottomleft",
       legend = "Snow Pillow",
       pch = 23, pt.bg = "deepskyblue", col = "black",
       pt.cex = 1.4, cex = 0.85, bg = "white", box.lwd = 0.5,
       inset = c(0.01, 0.02))

title(main = "Nechako Basin - Snow Pillow Locations",
      line = 1.5, cex.main = 1.4, font.main = 2)
mtext(sprintf(
  "WGS84 | Resolution: %.4f deg | Elevation: %.0f - %.0f m | Source: n50w130_dem",
  res(dem)[1], elev_range[1], elev_range[2]),
  side = 1, line = 3, cex = 0.8, col = "gray30")

# North arrow
arrows(xmin(dem) + dx * 0.06, ymin(dem) + dy * 0.05,
       xmin(dem) + dx * 0.06, ymin(dem) + dy * 0.10,
       length = 0.12, lwd = 2, col = "black")
text(xmin(dem) + dx * 0.06, ymin(dem) + dy * 0.115,
     "N", cex = 1.0, font = 2)

dev.off()
cat(sprintf("Snow pillows map saved: %s\n",
            file.path(out_dir, "nechako_elevation_snowpillows.pdf")))

cat("\n============================================================\n")
cat("COMPLETE\n")
cat("  nechako_elevation_map.pdf            -- elevation + hillshade\n")
cat("  nechako_slope_aspect.pdf             -- slope and aspect\n")
cat("  nechako_elevation_hydrostations.pdf  -- elevation + hydrometry stations\n")
cat("  nechako_elevation_snowpillows.pdf    -- elevation + snow pillows\n")
cat("============================================================\n")