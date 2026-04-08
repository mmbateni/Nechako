library(terra)
library(sf)
library(httr2)
library(ncdf4)
library(dplyr)
library(lubridate)

# ══════════════════════════════════════════════════════════════════════════════
#  CONFIG
# ══════════════════════════════════════════════════════════════════════════════
HOST       <- "https://pavics.ouranos.ca"
TOKEN      <- "***"#Use your Pavics API Token here
TARGET_CRS <- "EPSG:3005"
SHP_FILE   <- "D:/Nechako_Drought/Nechako/Spatial/nechakoBound_dissolve.shp"
OUT_DIR    <- "D:/Nechako_Drought/Nechako/RDRS"

out_full    <- file.path(OUT_DIR, "rdrs_tas_pr_nechako_1980_2018.nc")
out_tif_tas <- file.path(OUT_DIR, "rdrs_tas_nechako_3005.tif")
out_tif_pr  <- file.path(OUT_DIR, "rdrs_pr_nechako_3005.tif")
out_csv_mon <- file.path(OUT_DIR, "rdrs_basin_monthly_1980_2018.csv")
out_csv_ann <- file.path(OUT_DIR, "rdrs_basin_annual_1980_2018.csv")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

NCSS_BASE <- paste0(HOST,
                    "/twitcher/ows/proxy/thredds/ncss/grid/",
                    "datasets/reanalyses/day_RDRSv2.1_NAM.ncml"
)

# ══════════════════════════════════════════════════════════════════════════════
#  STEP 1: Convert Nechako bbox → rotated pole projection coords
# ══════════════════════════════════════════════════════════════════════════════
NP_LAT <- 31.758316
NP_LON <- 87.597031

geo_to_rotpole <- function(lon_geo, lat_geo, np_lon = NP_LON, np_lat = NP_LAT) {
  d2r     <- pi / 180
  lon_r   <- lon_geo * d2r;  lat_r   <- lat_geo * d2r
  nplon_r <- np_lon  * d2r;  nplat_r <- np_lat  * d2r
  sin_rlat <- sin(nplat_r)*sin(lat_r) +
    cos(nplat_r)*cos(lat_r)*cos(lon_r - nplon_r)
  rlat <- asin(pmax(-1, pmin(1, sin_rlat))) / d2r
  num  <- cos(lat_r) * sin(lon_r - nplon_r)
  den  <- cos(nplat_r)*sin(lat_r) -
    sin(nplat_r)*cos(lat_r)*cos(lon_r - nplon_r)
  rlon <- (atan2(num, den) / d2r) %% 360
  if (rlon < 300) rlon <- rlon + 360
  c(rlon = rlon, rlat = rlat)
}

corners <- expand.grid(lon = c(-128.0, -122.3), lat = c(52.5, 56.5))
rot     <- t(mapply(geo_to_rotpole, corners$lon, corners$lat))
minx <- min(rot[, "rlon"]) - 0.2;  maxx <- max(rot[, "rlon"]) + 0.2
miny <- min(rot[, "rlat"]) - 0.2;  maxy <- max(rot[, "rlat"]) + 0.2
message(sprintf("Rotated pole bbox — minx:%.3f maxx:%.3f miny:%.3f maxy:%.3f",
                minx, maxx, miny, maxy))

# ══════════════════════════════════════════════════════════════════════════════
#  STEP 2: Download tas + pr in one NCSS request
# ══════════════════════════════════════════════════════════════════════════════
if (file.exists(out_full)) {
  message("NetCDF already on disk — skipping download.")
} else {
  message("Downloading tas + pr  1980–2018 ...")
  resp <- request(NCSS_BASE) |>
    req_headers(Authorization = paste("Bearer", TOKEN)) |>
    req_url_query(
      var        = c("tas", "pr"),
      minx = minx, maxx = maxx, miny = miny, maxy = maxy,
      time_start = "1980-01-01T00:00:00Z",
      time_end   = "2018-12-31T00:00:00Z",
      accept     = "netcdf4-classic",
      .multi     = "explode"
    ) |>
    req_timeout(3600) |>
    req_error(is_error = \(r) FALSE) |>
    req_perform()
  
  message("HTTP status: ", resp_status(resp))
  if (resp_status(resp) != 200)
    stop("Download failed: ", substr(resp_body_string(resp), 1, 500))
  
  writeBin(resp_body_raw(resp), out_full)
  message(sprintf("Saved: %.1f MB", file.size(out_full) / 1e6))
}

# ══════════════════════════════════════════════════════════════════════════════
#  STEP 3: Get geographic extent from 2D lat/lon arrays in file
# ══════════════════════════════════════════════════════════════════════════════
nc        <- nc_open(out_full)
all_names <- c(names(nc$var), names(nc$dim))
message("Variables : ", paste(names(nc$var), collapse = ", "))
message("Dimensions: ", paste(names(nc$dim), collapse = ", "))

lon2d <- NULL;  lat2d <- NULL
for (nm in c("lon", "lon_1", "longitude", "nav_lon")) {
  if (nm %in% all_names) {
    tryCatch({ lon2d <- ncvar_get(nc, nm); message("lon: '", nm, "'"); break },
             error = function(e) NULL)
  }
}
for (nm in c("lat", "lat_1", "latitude", "nav_lat")) {
  if (nm %in% all_names) {
    tryCatch({ lat2d <- ncvar_get(nc, nm); message("lat: '", nm, "'"); break },
             error = function(e) NULL)
  }
}
nc_close(nc)

if (is.null(lon2d) || is.null(lat2d)) {
  message("2D lat/lon not found — using hard bbox")
  lon_min <- -128.2;  lon_max <- -122.1
  lat_min <-   52.3;  lat_max <-   56.7
} else {
  lon_min <- min(lon2d, na.rm=TRUE);  lon_max <- max(lon2d, na.rm=TRUE)
  lat_min <- min(lat2d, na.rm=TRUE);  lat_max <- max(lat2d, na.rm=TRUE)
}
message(sprintf("Geographic extent: lon [%.3f, %.3f]  lat [%.3f, %.3f]",
                lon_min, lon_max, lat_min, lat_max))

# ══════════════════════════════════════════════════════════════════════════════
#  STEP 4: Load watershed polygon
# ══════════════════════════════════════════════════════════════════════════════
if (!file.exists(SHP_FILE)) stop("Shapefile not found: ", SHP_FILE)
nechako_sf   <- st_read(SHP_FILE, quiet = TRUE)
nechako_3005 <- vect(st_transform(nechako_sf, TARGET_CRS))

# ── Helper: assign WGS84 extent, reproject, crop, mask ───────────────────────
process_var <- function(subds_name) {
  message(sprintf("\n── Processing: %s", subds_name))
  r <- rast(out_full, subds = subds_name)
  message(sprintf("  Raw: %d layers | %d x %d cells", nlyr(r), nrow(r), ncol(r)))
  crs(r) <- "EPSG:4326"
  ext(r) <- c(lon_min, lon_max, lat_min, lat_max)
  message("  Reprojecting to BC Albers...")
  r_3005   <- project(r, TARGET_CRS, method = "bilinear")
  r_crop   <- crop(r_3005, ext(nechako_3005), snap = "out")
  r_masked <- mask(r_crop, nechako_3005)
  message(sprintf("  Done: %d x %d cells | Extent: %.0f %.0f %.0f %.0f",
                  nrow(r_masked), ncol(r_masked),
                  xmin(r_masked), xmax(r_masked), ymin(r_masked), ymax(r_masked)))
  r_masked
}

# ══════════════════════════════════════════════════════════════════════════════
#  STEP 5: Process temperature (K → °C) and precipitation (kg m-2 s-1 → mm/day)
# ══════════════════════════════════════════════════════════════════════════════
tas_final <- process_var("tas") - 273.15
pr_final  <- process_var("pr")  * 86400

# Assign meaningful layer names = dates (1980-01-01 … 2018-12-31)
dates_all <- seq(as.Date("1980-01-01"), as.Date("2018-12-31"), by = "day")
stopifnot(nlyr(tas_final) == length(dates_all))
names(tas_final) <- as.character(dates_all)
names(pr_final)  <- as.character(dates_all)

# Quick validation
t1 <- as.numeric(global(tas_final[[1]], "range", na.rm = TRUE))
p1 <- as.numeric(global(pr_final[[1]],  "range", na.rm = TRUE))
message(sprintf("tas Day-1 range: %.1f – %.1f °C", t1[1], t1[2]))
message(sprintf("pr  Day-1 range: %.2f – %.2f mm/day", p1[1], p1[2]))

writeRaster(tas_final, out_tif_tas, overwrite = TRUE)
writeRaster(pr_final,  out_tif_pr,  overwrite = TRUE)
message("Saved: ", out_tif_tas)
message("Saved: ", out_tif_pr)

# ══════════════════════════════════════════════════════════════════════════════
#  STEP 6: Area-weighted basin-average daily time series
#  cellSize() in projected CRS gives true cell areas in m² — correct for
#  irregular grids after reprojection
# ══════════════════════════════════════════════════════════════════════════════
message("\nComputing area-weighted basin averages...")

cell_areas  <- cellSize(tas_final[[1]], unit = "m")   # m² per cell
cell_areas  <- mask(cell_areas, nechako_3005)
total_area  <- global(cell_areas, "sum", na.rm = TRUE)[[1]]

area_weighted_mean <- function(r_stack) {
  # Multiply each layer by cell area, sum, divide by total area
  weighted <- r_stack * cell_areas
  vals     <- global(weighted, "sum", na.rm = TRUE)[[1]] / total_area
  vals
}

tas_daily_mean <- area_weighted_mean(tas_final)   # vector length 14245
pr_daily_mean  <- area_weighted_mean(pr_final)

daily_df <- data.frame(
  date  = dates_all,
  year  = year(dates_all),
  month = month(dates_all),
  tas_C      = tas_daily_mean,
  pr_mm_day  = pr_daily_mean
)

# ══════════════════════════════════════════════════════════════════════════════
#  STEP 7: Monthly aggregation
#  tas  → monthly mean (°C)
#  pr   → monthly total (mm/month) = mean mm/day × days in month
# ══════════════════════════════════════════════════════════════════════════════
monthly_df <- daily_df |>
  group_by(year, month) |>
  summarise(
    date          = as.Date(sprintf("%d-%02d-01", first(year), first(month))),
    tas_mean_C    = mean(tas_C,     na.rm = TRUE),
    pr_total_mm   = sum(pr_mm_day,  na.rm = TRUE),   # sum of mm/day over days = mm/month
    n_days        = n(),
    .groups       = "drop"
  ) |>
  arrange(date) |>
  select(date, year, month, n_days, tas_mean_C, pr_total_mm)

message(sprintf("Monthly rows: %d  (expect %d)", nrow(monthly_df),
                length(seq(as.Date("1980-01-01"), as.Date("2018-12-01"), by = "month"))))

# ══════════════════════════════════════════════════════════════════════════════
#  STEP 8: Annual aggregation
#  tas  → annual mean (°C)
#  pr   → annual total (mm/year) = sum of monthly totals
# ══════════════════════════════════════════════════════════════════════════════
annual_df <- monthly_df |>
  group_by(year) |>
  summarise(
    tas_mean_C   = mean(tas_mean_C,  na.rm = TRUE),
    pr_total_mm  = sum(pr_total_mm,  na.rm = TRUE),
    n_months     = n(),
    .groups      = "drop"
  ) |>
  arrange(year)

message(sprintf("Annual rows: %d  (expect 39)", nrow(annual_df)))

# Sanity check — typical Nechako values from literature:
# mean annual temp ~ 2°C,  mean annual precip ~ 826 mm
message(sprintf("Mean annual temp : %.2f °C  (literature ~2 °C)",
                mean(annual_df$tas_mean_C)))
message(sprintf("Mean annual precip: %.0f mm  (literature ~826 mm)",
                mean(annual_df$pr_total_mm)))

# ══════════════════════════════════════════════════════════════════════════════
#  STEP 9: Write CSVs
# ══════════════════════════════════════════════════════════════════════════════
write.csv(monthly_df, out_csv_mon, row.names = FALSE)
write.csv(annual_df,  out_csv_ann, row.names = FALSE)
message("Saved: ", out_csv_mon)
message("Saved: ", out_csv_ann)

# ══════════════════════════════════════════════════════════════════════════════
#  STEP 10: Plots
# ══════════════════════════════════════════════════════════════════════════════
par(mfrow = c(2, 2))

# Map — temperature day 1
plot(tas_final[[1]],
     main = "RDRS v2.1 — Tas 1980-01-01 (°C)\nNechako, BC Albers")
lines(nechako_3005, col = "black", lwd = 2)

# Map — precipitation day 1
plot(pr_final[[1]],
     main = "RDRS v2.1 — Pr 1980-01-01 (mm/day)\nNechako, BC Albers")
lines(nechako_3005, col = "black", lwd = 2)

# Time series — monthly temperature
plot(monthly_df$date, monthly_df$tas_mean_C,
     type = "l", col = "tomato", lwd = 1.2,
     xlab = "", ylab = "Mean Temp (°C)",
     main = "Basin-Average Monthly Temperature")
abline(h = 0, lty = 2, col = "grey60")

# Time series — monthly precipitation
plot(monthly_df$date, monthly_df$pr_total_mm,
     type = "l", col = "steelblue", lwd = 1.2,
     xlab = "", ylab = "Precipitation (mm/month)",
     main = "Basin-Average Monthly Precipitation")

par(mfrow = c(1, 1))

message("\n══ All done ══")
message("  GeoTIFF temperature  : ", out_tif_tas)
message("  GeoTIFF precipitation: ", out_tif_pr)
message("  Monthly CSV          : ", out_csv_mon)
message("  Annual  CSV          : ", out_csv_ann)
