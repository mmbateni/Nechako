####################################################################################
# NECHAKO BASIN — ELEVATION & CLIMATE MAPS (v12 — UNITS IN TITLE, IMPROVED PRECIP)
####################################################################################
# Fixes in this version:
# 1. Units moved to legend title (e.g., "Elevation (m)" not "500 m")
# 2. Precipitation: finer classes at lower values, coarser at higher
# 3. Legend bars: narrow (5% width), positioned at 72-77% of plot width
# 4. All legends use par("usr") coordinates with visible numbers
# 5. Precipitation: ERA5-Land m/DAY -> mm/year
# 6. Temperature: ERA5-Land Kelvin -> °C
####################################################################################

library(terra)
library(sf)
library(readr)

####################################################################################
# PATHS
####################################################################################
dem_dir      <- "D:/Nechako_Drought/Nechako/Spatial/nechako-dem"
basin_shp    <- "D:/Nechako_Drought/Nechako/Spatial/nechakoBound_dissolve.shp"
stations_csv <- "D:/Nechako_Drought/Nechako/Hydrology/data_retrievalWaterSurveyofCanada/stations.csv"
era5_dir     <- "D:/Nechako_Drought/Nechako/monthly_data_direct"
out_dir      <- "D:/Nechako_Drought/Nechako/elevation"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
bg_col <- "#d6eaf8"

####################################################################################
# LOAD DEM
####################################################################################
cat("Loading DEM...\n")
dem <- rast(dem_dir)
cat(sprintf("  Extent    : %s\n", as.character(ext(dem))))
cat(sprintf("  Resolution: %.6f deg\n", res(dem)[1]))

####################################################################################
# LOAD BASIN BOUNDARY
####################################################################################
cat("Loading basin boundary...\n")
basin_sf <- st_read(basin_shp, quiet = TRUE)
basin_v  <- vect(basin_shp)

if (!same.crs(basin_v, crs(dem))) {
  basin_v  <- project(basin_v, crs(dem))
  basin_sf <- st_transform(basin_sf, st_crs(crs(dem)))
}

cat("Masking DEM to basin boundary...\n")
dem <- mask(dem, basin_v, touches = TRUE)
elev_range <- as.numeric(global(dem, c("min","max"), na.rm = TRUE))
cat(sprintf("  Elevation : %.0f - %.0f m\n", elev_range[1], elev_range[2]))

####################################################################################
# HILLSHADE
####################################################################################
cat("Computing hillshade...\n")
dem_m     <- project(dem, "EPSG:3005", method = "bilinear")
slope_r   <- terrain(dem_m, "slope", unit = "radians")
aspect_r  <- terrain(dem_m, "aspect", unit = "radians")
hill      <- shade(slope_r, aspect_r, angle = 40, direction = 315, normalize = TRUE)
hill      <- project(hill, crs(dem), method = "bilinear")
hill      <- resample(hill, dem, method = "bilinear")

####################################################################################
# COLOUR PALETTES
####################################################################################
n_col <- 512
elev_cols  <- c(
  colorRampPalette(c("#336600", "#669933", "#99cc66"))(round(n_col * 0.20)),
  colorRampPalette(c("#cccc66", "#cc9933", "#996633"))(round(n_col * 0.35)),
  colorRampPalette(c("#996633", "#cc9966", "#d4b896"))(round(n_col * 0.25)),
  colorRampPalette(c("#d4b896", "#e8dcc8", "#f5f0e8", "#ffffff"))(round(n_col * 0.20))
)

n_hs <- 256
hill_overlay_cols <- sapply(seq(0, 1, length.out = n_hs), function(h) {
  if (h < 0.5) rgb(0, 0, 0, alpha = (0.5 - h) * 0.9)
  else         rgb(1, 1, 1, alpha = (h - 0.5) * 0.6)
})

####################################################################################
# HELPER FUNCTION: DRAW LEGEND BAR WITH NUMBERS
# Units in title, NOT with numbers
####################################################################################
draw_legend_bar <- function(data_range, col_palette, label_text, unit_text, 
                            custom_breaks = NULL, custom_labels = NULL) {
  usr <- par("usr")
  
  # Legend position: 72-77% of plot width (narrow bar - 5% width)
  #                  55-90% of plot height (upper portion)
  legend_x   <- usr[1] + (usr[2] - usr[1]) * 0.72
  legend_x2  <- usr[1] + (usr[2] - usr[1]) * 0.77
  legend_y0  <- usr[3] + (usr[4] - usr[3]) * 0.55
  legend_y1  <- usr[3] + (usr[4] - usr[3]) * 0.90
  
  # Draw color bar
  n_leg <- 150
  leg_y <- seq(legend_y0, legend_y1, length.out = n_leg + 1)
  leg_cols <- colorRampPalette(col_palette)(n_leg)
  
  for (i in seq_len(n_leg))
    rect(legend_x, leg_y[i], legend_x2, leg_y[i+1],
         col = leg_cols[i], border = NA)
  
  rect(legend_x, legend_y0, legend_x2, legend_y1, border = "black", lwd = 0.8)
  
  # Draw tick marks and numbers (NO units with numbers)
  if (!is.null(custom_breaks)) {
    leg_vals <- custom_breaks
  } else {
    leg_vals <- pretty(data_range, n = 6)
    leg_vals <- leg_vals[leg_vals >= data_range[1] & leg_vals <= data_range[2]]
  }
  
  leg_ypos <- legend_y0 +
    (leg_vals - data_range[1]) / diff(data_range) *
    (legend_y1 - legend_y0)
  
  # Tick marks (40% of bar width)
  tick_len <- (legend_x2 - legend_x) * 0.40
  
  for (i in seq_along(leg_vals)) {
    segments(legend_x2, leg_ypos[i], legend_x2 + tick_len, leg_ypos[i], lwd = 0.7)
    # Numbers WITHOUT units
    text(legend_x2 + tick_len * 1.05, leg_ypos[i], 
         leg_vals[i], adj = 0, cex = 0.85)
  }
  
  # Label above the color bar WITH UNITS
  text((legend_x + legend_x2) / 2,
       legend_y1 + (legend_y1 - legend_y0) * 0.08,
       paste0(label_text, " (", unit_text, ")"), 
       adj = c(0.5, 0), cex = 0.95, font = 2)
}

####################################################################################
# MAP 1 — MAIN ELEVATION MAP
####################################################################################
cat("Rendering elevation map...\n")
pdf_path <- file.path(out_dir, "nechako_elevation_map.pdf")
cat(sprintf("Rendering -> %s\n", basename(pdf_path)))

pdf(pdf_path, width = 11, height = 9, family = "Helvetica")
par(mar = c(5.5, 4.5, 3.5, 11.0),
    oma = c(0, 0, 0, 0),
    bg  = bg_col)

plot(dem,
     col    = elev_cols,
     legend = FALSE,
     axes   = TRUE,
     box    = FALSE,
     xlab   = "Longitude (W)",
     ylab   = "Latitude (N)",
     colNA  = bg_col,
     main   = "")

plot(hill,
     col    = hill_overlay_cols,
     legend = FALSE,
     axes   = FALSE,
     box    = FALSE,
     colNA  = NA,
     add    = TRUE)

plot(st_geometry(basin_sf), add = TRUE,
     col = NA, border = "black", lwd = 2.2)

grid(nx = NULL, ny = NULL, col = "white", lty = 2, lwd = 0.35)

draw_legend_bar(elev_range, elev_cols, "Elevation", "m")

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

####################################################################################
# SLOPE AND ASPECT
####################################################################################
cat("Rendering slope + aspect map...\n")
slope_deg  <- terrain(dem_m, "slope", unit = "degrees")
aspect_deg <- terrain(dem_m, "aspect", unit = "degrees")
slope_deg  <- project(slope_deg, crs(dem))
aspect_deg <- project(aspect_deg, crs(dem))

pdf(file.path(out_dir, "nechako_slope_aspect.pdf"),
    width = 12, height = 6, family = "Helvetica")
par(mfrow = c(1, 2), mar = c(3, 3, 2.5, 4), oma = c(0, 0, 2, 0))

plot(slope_deg, main = "Slope (degrees)",
     col = hcl.colors(101, "YlOrRd", rev = FALSE),
     axes = TRUE, box = FALSE,
     plg = list(title = "deg", cex = 0.85, width = 0.4), cex.main = 1.1)
plot(st_geometry(basin_sf), add = TRUE, col = NA, border = "black", lwd = 2)

plot(aspect_deg, main = "Aspect (degrees from North)",
     col = hcl.colors(101, "RdYlBu", rev = FALSE),
     axes = TRUE, box = FALSE,
     plg = list(title = "deg", cex = 0.85, width = 0.4), cex.main = 1.1)
plot(st_geometry(basin_sf), add = TRUE, col = NA, border = "black", lwd = 2)

mtext("Nechako Basin - Terrain Derivatives",
      outer = TRUE, cex = 1.2, font = 2, line = 0.3)
dev.off()
cat(sprintf("Slope + aspect saved: %s\n",
            file.path(out_dir, "nechako_slope_aspect.pdf")))

####################################################################################
# MAP 3 — ELEVATION + HYDROMETRY STATIONS
####################################################################################
cat("Rendering elevation + hydrometry stations map...\n")
stations_df <- read_csv(stations_csv)
stations_df <- stations_df[!is.na(stations_df$LATITUDE) & !is.na(stations_df$LONGITUDE), ]

pdf(file.path(out_dir, "nechako_elevation_hydrostations.pdf"),
    width = 11, height = 9, family = "Helvetica")
par(mfrow = c(1, 1),
    mar   = c(5.5, 4.5, 3.5, 11.0),
    oma   = c(0, 0, 0, 0),
    bg    = bg_col)

dem_ext <- as.vector(ext(dem))
xlims   <- dem_ext[1:2]
ylims   <- dem_ext[3:4]

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

plot(hill,
     col    = hill_overlay_cols,
     legend = FALSE,
     axes   = FALSE,
     box    = FALSE,
     colNA  = NA,
     add    = TRUE)

grid(nx = NULL, ny = NULL, col = "white", lty = 2, lwd = 0.35)

points(stations_df$LONGITUDE, stations_df$LATITUDE,
       pch = 21, bg = "yellow", col = "black", cex = 1.0, lwd = 0.8)

text(stations_df$LONGITUDE, stations_df$LATITUDE,
     labels = stations_df$STATION_NUMBER,
     pos = 3, cex = 0.50, col = "black", font = 1,
     offset = 0.35)

draw_legend_bar(elev_range, elev_cols, "Elevation", "m")

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

arrows(xmin(dem) + dx * 0.06, ymin(dem) + dy * 0.05,
       xmin(dem) + dx * 0.06, ymin(dem) + dy * 0.10,
       length = 0.12, lwd = 2, col = "black")
text(xmin(dem) + dx * 0.06, ymin(dem) + dy * 0.115,
     "N", cex = 1.0, font = 2)

dev.off()
cat(sprintf("Hydrometry stations map saved: %s\n",
            file.path(out_dir, "nechako_elevation_hydrostations.pdf")))

####################################################################################
# MAP 4 — ELEVATION + SNOW PILLOWS
####################################################################################
cat("Rendering elevation + snow pillows map...\n")
snow_pillows <- data.frame(
  id   = c("1B08P", "1B01P", "1B02P"),
  name = c("Mt. Pondosy", "Mount Wells", "Tahtsa Lake"),
  lat  = c(53.15, 53.72, 53.59),
  lon  = c(-126.88, -126.42, -127.64)
)

pdf(file.path(out_dir, "nechako_elevation_snowpillows.pdf"),
    width = 11, height = 9, family = "Helvetica")
par(mfrow = c(1, 1),
    mar   = c(5.5, 4.5, 3.5, 11.0),
    oma   = c(0, 0, 0, 0),
    bg    = bg_col)

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

plot(hill,
     col    = hill_overlay_cols,
     legend = FALSE,
     axes   = FALSE,
     box    = FALSE,
     colNA  = NA,
     add    = TRUE)

grid(nx = NULL, ny = NULL, col = "white", lty = 2, lwd = 0.35)

points(snow_pillows$lon, snow_pillows$lat,
       pch = 23, bg = "deepskyblue", col = "black", cex = 1.6, lwd = 1.2)

id_pos   <- c(3, 3, 4)
name_pos <- c(1, 1, 4)

for (i in seq_len(nrow(snow_pillows))) {
  text(snow_pillows$lon[i], snow_pillows$lat[i],
       labels = snow_pillows$id[i],
       pos = id_pos[i], cex = 0.85, col = "black", font = 2, offset = 0.5)
}

for (i in seq_len(nrow(snow_pillows))) {
  text(snow_pillows$lon[i], snow_pillows$lat[i],
       labels = snow_pillows$name[i],
       pos = name_pos[i], cex = 0.75, col = "gray20", font = 3,
       offset = if (name_pos[i] == 4) 3.5 else 0.5)
}

draw_legend_bar(elev_range, elev_cols, "Elevation", "m")

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

arrows(xmin(dem) + dx * 0.06, ymin(dem) + dy * 0.05,
       xmin(dem) + dx * 0.06, ymin(dem) + dy * 0.10,
       length = 0.12, lwd = 2, col = "black")
text(xmin(dem) + dx * 0.06, ymin(dem) + dy * 0.115,
     "N", cex = 1.0, font = 2)

dev.off()
cat(sprintf("Snow pillows map saved: %s\n",
            file.path(out_dir, "nechako_elevation_snowpillows.pdf")))

####################################################################################
# MAP 5 — MEAN ANNUAL PRECIPITATION (ERA5-Land 1950-2025)
# IMPROVED CLASSIFICATION: finer at lower values, coarser at higher
####################################################################################
cat("\nRendering mean annual precipitation map...\n")

cat("Loading precipitation data...\n")
prec_file <- file.path(era5_dir, "total_precipitation_monthly.nc")
prec_rast <- rast(prec_file)

if (!same.crs(crs(prec_rast), crs(dem))) {
  prec_rast <- project(prec_rast, crs(dem), method = "bilinear")
}

prec_rast <- mask(prec_rast, basin_v, touches = TRUE)

# ERA5-Land: meters PER DAY -> convert to mm/year
days_in_month <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

n_months <- nlyr(prec_rast)
n_complete_years <- floor(n_months / 12)
cat(sprintf("  Processing %d months (%d complete years)...\n", n_months, n_complete_years))

annual_prec_list <- list()
for (yr in 1:n_complete_years) {
  start_month <- (yr - 1) * 12 + 1
  end_month <- yr * 12
  
  monthly_mm_list <- list()
  for (m in 1:12) {
    month_layer <- prec_rast[[start_month + m - 1]]
    monthly_mm_list[[m]] <- month_layer * days_in_month[m] * 1000
  }
  
  year_total <- sum(rast(monthly_mm_list))
  annual_prec_list[[yr]] <- year_total
}

annual_prec_stack <- rast(annual_prec_list)
mean_annual_prec <- mean(annual_prec_stack, na.rm = TRUE)

prec_range <- as.numeric(global(mean_annual_prec, c("min", "max"), na.rm = TRUE))
cat(sprintf("  Mean annual precipitation: %.0f - %.0f mm/year\n", 
            prec_range[1], prec_range[2]))

# ===== IMPROVED NON-EQUAL BREAKS FOR PRECIPITATION =====
# Finer classes at LOWER values (where most of basin falls)
# Coarser classes at HIGHER values (small areas with extreme precip)
# Using percentile-based breaks for better spatial distribution
prec_breaks <- c(
  prec_range[1],                          # minimum (~500 mm)
  prec_range[1] + diff(prec_range) * 0.05,   # 5%
  prec_range[1] + diff(prec_range) * 0.10,   # 10%
  prec_range[1] + diff(prec_range) * 0.15,   # 15%
  prec_range[1] + diff(prec_range) * 0.20,   # 20%
  prec_range[1] + diff(prec_range) * 0.30,   # 30%
  prec_range[1] + diff(prec_range) * 0.40,   # 40%
  prec_range[1] + diff(prec_range) * 0.50,   # 50%
  prec_range[1] + diff(prec_range) * 0.65,   # 65%
  prec_range[1] + diff(prec_range) * 0.80,   # 80%
  prec_range[2]                           # maximum (~2800 mm)
)

prec_breaks <- round(prec_breaks, 0)
prec_breaks <- unique(prec_breaks)

cat(sprintf("  Precipitation breaks: %s mm/yr\n", 
            paste(prec_breaks, collapse = ", ")))

n_prec_cols <- length(prec_breaks) - 1
prec_cols <- colorRampPalette(c("#ffffcc", "#ffeda0", "#fed976", "#feb24c", 
                                "#fd8d3c", "#fc4e2a", "#e31a1c", "#bd0026", 
                                "#800026"))(n_prec_cols)

# Classify raster using reclassification matrix
rcl <- cbind(
  prec_breaks[-length(prec_breaks)],
  prec_breaks[-1],
  1:(length(prec_breaks)-1)
)
prec_classified <- classify(mean_annual_prec, rcl = rcl, include.lowest = TRUE)

pdf_path_prec <- file.path(out_dir, "nechako_mean_annual_precipitation.pdf")
cat(sprintf("Rendering -> %s\n", basename(pdf_path_prec)))

pdf(pdf_path_prec, width = 11, height = 9, family = "Helvetica")
par(mar = c(5.5, 4.5, 3.5, 11.0),
    oma = c(0, 0, 0, 0),
    bg  = bg_col)

plot(prec_classified,
     col    = prec_cols,
     legend = FALSE,
     axes   = TRUE,
     box    = FALSE,
     xlab   = "Longitude (W)",
     ylab   = "Latitude (N)",
     colNA  = bg_col,
     main   = "")

plot(hill,
     col    = hill_overlay_cols,
     legend = FALSE,
     axes   = FALSE,
     box    = FALSE,
     colNA  = NA,
     add    = TRUE)

plot(st_geometry(basin_sf), add = TRUE,
     col = NA, border = "black", lwd = 2.2)

grid(nx = NULL, ny = NULL, col = "white", lty = 2, lwd = 0.35)

draw_legend_bar(prec_range, prec_cols, "Mean Annual Precipitation", "mm/yr",
                custom_breaks = prec_breaks)

title(main = "Nechako Basin - Mean Annual Precipitation (1950-2025)",
      line = 1.5, cex.main = 1.4, font.main = 2)

mtext(sprintf(
  "WGS84 | ERA5-Land | Period: 1950-2025 | Range: %.0f - %.0f mm/year",
  prec_range[1], prec_range[2]),
  side = 1, line = 3, cex = 0.8, col = "gray30")

arrows(xmin(dem) + dx * 0.06, ymin(dem) + dy * 0.05,
       xmin(dem) + dx * 0.06, ymin(dem) + dy * 0.10,
       length = 0.12, lwd = 2, col = "black")
text(xmin(dem) + dx * 0.06, ymin(dem) + dy * 0.115,
     "N", cex = 1.0, font = 2)

dev.off()
cat(sprintf("Mean annual precipitation map saved: %s\n", pdf_path_prec))

####################################################################################
# MAP 6 — MEAN ANNUAL TEMPERATURE (ERA5-Land 1950-2025)
####################################################################################
cat("\nRendering mean annual temperature map...\n")

cat("Loading temperature data...\n")
temp_file <- file.path(era5_dir, "2m_temperature_monthly.nc")
temp_rast <- rast(temp_file)

if (!same.crs(crs(temp_rast), crs(dem))) {
  temp_rast <- project(temp_rast, crs(dem), method = "bilinear")
}

temp_rast <- mask(temp_rast, basin_v, touches = TRUE)

temp_rast <- temp_rast - 273.15

n_months <- nlyr(temp_rast)
n_complete_years <- floor(n_months / 12)
cat(sprintf("  Processing %d months (%d complete years)...\n", n_months, n_complete_years))

annual_temp_list <- list()
for (yr in 1:n_complete_years) {
  start_month <- (yr - 1) * 12 + 1
  end_month <- yr * 12
  year_mean <- mean(temp_rast[[start_month:end_month]])
  annual_temp_list[[yr]] <- year_mean
}

annual_temp_stack <- rast(annual_temp_list)
mean_annual_temp <- mean(annual_temp_stack, na.rm = TRUE)

temp_range <- as.numeric(global(mean_annual_temp, c("min", "max"), na.rm = TRUE))
cat(sprintf("  Mean annual temperature: %.1f - %.1f °C\n", 
            temp_range[1], temp_range[2]))

temp_cols <- colorRampPalette(c("#2c7bb6", "#abd9e9", "#ffffbf", 
                                "#fed976", "#d73027"))(256)

pdf_path_temp <- file.path(out_dir, "nechako_mean_annual_temperature.pdf")
cat(sprintf("Rendering -> %s\n", basename(pdf_path_temp)))

pdf(pdf_path_temp, width = 11, height = 9, family = "Helvetica")
par(mar = c(5.5, 4.5, 3.5, 11.0),
    oma = c(0, 0, 0, 0),
    bg  = bg_col)

plot(mean_annual_temp,
     col    = temp_cols,
     legend = FALSE,
     axes   = TRUE,
     box    = FALSE,
     xlab   = "Longitude (W)",
     ylab   = "Latitude (N)",
     colNA  = bg_col,
     main   = "")

plot(hill,
     col    = hill_overlay_cols,
     legend = FALSE,
     axes   = FALSE,
     box    = FALSE,
     colNA  = NA,
     add    = TRUE)

plot(st_geometry(basin_sf), add = TRUE,
     col = NA, border = "black", lwd = 2.2)

grid(nx = NULL, ny = NULL, col = "white", lty = 2, lwd = 0.35)

draw_legend_bar(temp_range, temp_cols, "Mean Annual Temperature", "°C")

title(main = "Nechako Basin - Mean Annual Temperature (1950-2025)",
      line = 1.5, cex.main = 1.4, font.main = 2)

mtext(sprintf(
  "WGS84 | ERA5-Land | Period: 1950-2025 | Range: %.1f - %.1f °C",
  temp_range[1], temp_range[2]),
  side = 1, line = 3, cex = 0.8, col = "gray30")

arrows(xmin(dem) + dx * 0.06, ymin(dem) + dy * 0.05,
       xmin(dem) + dx * 0.06, ymin(dem) + dy * 0.10,
       length = 0.12, lwd = 2, col = "black")
text(xmin(dem) + dx * 0.06, ymin(dem) + dy * 0.115,
     "N", cex = 1.0, font = 2)

dev.off()
cat(sprintf("Mean annual temperature map saved: %s\n", pdf_path_temp))

####################################################################################
# MAP 7 — COMBINED PRECIPITATION AND TEMPERATURE (SIDE-BY-SIDE)
####################################################################################
cat("\nRendering combined precipitation + temperature map...\n")

pdf_path_combined <- file.path(out_dir, "nechako_climate_combined.pdf")
cat(sprintf("Rendering -> %s\n", basename(pdf_path_combined)))

pdf(pdf_path_combined, width = 14, height = 7, family = "Helvetica")
par(mfrow = c(1, 2), 
    mar = c(4, 4, 3, 4),
    oma = c(0, 0, 2, 0),
    bg  = bg_col)

plot(prec_classified,
     col    = prec_cols,
     legend = TRUE,
     axes   = TRUE,
     box    = FALSE,
     xlab   = "Longitude (W)",
     ylab   = "Latitude (N)",
     colNA  = bg_col,
     main   = "Mean Annual Precipitation",
     plg    = list(title = "mm/yr", cex = 0.8, width = 0.4))

plot(hill, col = hill_overlay_cols, legend = FALSE, 
     axes = FALSE, box = FALSE, colNA = NA, add = TRUE)
plot(st_geometry(basin_sf), add = TRUE, col = NA, border = "black", lwd = 2)

plot(mean_annual_temp,
     col    = temp_cols,
     legend = TRUE,
     axes   = TRUE,
     box    = FALSE,
     xlab   = "Longitude (W)",
     ylab   = "",
     colNA  = bg_col,
     main   = "Mean Annual Temperature",
     plg    = list(title = "°C", cex = 0.8, width = 0.4))

plot(hill, col = hill_overlay_cols, legend = FALSE, 
     axes = FALSE, box = FALSE, colNA = NA, add = TRUE)
plot(st_geometry(basin_sf), add = TRUE, col = NA, border = "black", lwd = 2)

mtext("Nechako Basin - Climate Variables (ERA5-Land, 1950-2025)", 
      outer = TRUE, cex = 1.3, font = 2, line = 0.3)

dev.off()
cat(sprintf("Combined climate map saved: %s\n", pdf_path_combined))

####################################################################################
# SUMMARY
####################################################################################
cat("\n============================================================\n")
cat("ALL MAPS COMPLETE\n")
cat("============================================================\n")
cat("ELEVATION MAPS:\n")
cat("  nechako_elevation_map.pdf            -- elevation + hillshade\n")
cat("  nechako_slope_aspect.pdf             -- slope and aspect\n")
cat("  nechako_elevation_hydrostations.pdf  -- elevation + hydrometry stations\n")
cat("  nechako_elevation_snowpillows.pdf    -- elevation + snow pillows\n")
cat("\nCLIMATE MAPS:\n")
cat("  nechako_mean_annual_precipitation.pdf  -- precipitation (mm/yr, improved classes)\n")
cat("  nechako_mean_annual_temperature.pdf    -- temperature (°C)\n")
cat("  nechako_climate_combined.pdf           -- both variables side-by-side\n")
cat("============================================================\n")