# ======= Nechako Basin ERA data Retrieval workflow using a LOCAL shapefile =======
# Packages
library(terra)
library(sf)
library(dplyr)
library(ecmwfr)

# 0) READ LOCAL BASIN SHAPEFILE (MULTIPOLYGON -> SINGLE OUTLINE)
# --------------------------------------------------------------
basin_shp <- "Nechako.shp"   
stopifnot(file.exists(basin_shp))
nechako_raw <- st_read(basin_shp, quiet = TRUE)
stopifnot(!is.na(st_geometry(nechako_raw)))

# If the shapefile contains multiple parts, dissolve to one polygon
# st_make_valid helps avoid topology errors before union
nechako_basin <- nechako_raw %>%
  st_make_valid() %>%
  dplyr::summarise(.groups = "drop") %>%   # dissolve
  st_make_valid()

# If CRS is unknown, set it to BC Albers (EPSG:3005) or WGS84 as appropriate.
# st_crs(nechako_basin) <- 3005     # if your local SHP is in BC Albers
# st_crs(nechako_basin) <- 4326     # if your local SHP is already lon/lat

stopifnot(!is.na(st_crs(nechako_basin)))   # we need a known CRS

# Prepare WGS84 version for bbox/plot/export
nechako_basin_ll <- st_transform(nechako_basin, 4326)

plot(st_geometry(nechako_basin_ll), col = "lightblue", border = "black", lwd = 1.5)
title("Nechako Basin (Local SHP → dissolved, WGS84)")

# 1) WRITE OUTPUTS (clean WGS84 shapefile + preview)
# --------------------------------------------------
out_dir   <- "out"
out_layer <- "Nechako_Basin_EPSG4326"
dir.create(out_dir, showWarnings = FALSE)

# Remove old files (optional)
unlink(file.path(out_dir, paste0(out_layer, c(".shp",".shx",".dbf",".prj",".cpg"))))

st_write(nechako_basin_ll, dsn = out_dir, layer = out_layer,
         driver = "ESRI Shapefile", delete_layer = TRUE)

# Quick PNG preview
png_filename <- file.path(out_dir, "Nechako_Basin_preview.png")
png(png_filename, width = 1200, height = 900, res = 120)
plot(nechako_basin_ll, border = "blue", lwd = 2,
     main = "Nechako River Basin – EPSG:4326")
dev.off()
message("Preview map saved: ", png_filename)

# 2) QA SUMMARY: bbox, area, validity
# -----------------------------------
message("\n=== Nechako Basin: Summary ===")
print(nechako_basin_ll)
message("CRS (EPSG): ", st_crs(nechako_basin_ll)$epsg)

bbox_ll <- st_bbox(nechako_basin_ll)
print(bbox_ll)

pad_deg <- 0.1
era5_bbox_padded <- bbox_ll
era5_bbox_padded["xmin"] <- era5_bbox_padded["xmin"] - pad_deg
era5_bbox_padded["ymin"] <- era5_bbox_padded["ymin"] - pad_deg
era5_bbox_padded["xmax"] <- era5_bbox_padded["xmax"] + pad_deg
era5_bbox_padded["ymax"] <- era5_bbox_padded["ymax"] + pad_deg
message("\nERA5-Land padded bbox (deg): ",
        paste0("lon:[", round(era5_bbox_padded["xmin"], 4), ", ",
               round(era5_bbox_padded["xmax"], 4), "]  "),
        paste0("lat:[", round(era5_bbox_padded["ymin"], 4), ", ",
               round(era5_bbox_padded["ymax"], 4), "]"))

# Area in BC Albers (m^2 -> km^2)
nechako_basin_bc <- st_transform(nechako_basin_ll, 3005)
area_km2 <- as.numeric(st_area(nechako_basin_bc)) / 1e6
message("Area (BC Albers): ", format(round(area_km2, 1), big.mark = ","), " km^2")

# Valid?
is_valid <- sf::st_is_valid(nechako_basin_ll)
message("Geometry valid? ", is_valid)

# 3) PREP EXTENTS FOR RASTER WORK AND ERA5 REQUEST
# ------------------------------------------------
# terra::ext() expects a Spat*; convert sf -> SpatVector
nechako_vect <- terra::vect(nechako_basin_ll)
era5_ext <- terra::ext(nechako_vect)
print(era5_ext)

# CDS expects area vector: [North, West, South, East]
area_for_request <- c(
  as.numeric(era5_bbox_padded["ymax"]),
  as.numeric(era5_bbox_padded["xmin"]),
  as.numeric(era5_bbox_padded["ymin"]),
  as.numeric(era5_bbox_padded["xmax"])
)

# 4) ECMWF CDS AUTH + ERA5-LAND HOURLY REQUEST (robust)
# -----------------------------------------------------
options(keyring_backend = "file")

conf_path <- file.path(getwd(), ".cdsapirc")
stopifnot(file.exists(conf_path))
lines <- readLines(conf_path)

# Parse `url:` and `key: uid:token`
url_line <- grep("^\\s*url\\s*:", lines, value = TRUE)
key_line <- grep("^\\s*key\\s*:", lines, value = TRUE)
stopifnot(length(url_line) == 1, length(key_line) == 1)

cds_url <- trimws(sub("^\\s*url\\s*:\\s*", "", url_line))
uid_token <- trimws(sub("^\\s*key\\s*:\\s*", "", key_line))
parts <- strsplit(uid_token, ":", fixed = TRUE)[[1]]
stopifnot(length(parts) == 2)
cds_uid   <- parts[1]
cds_token <- parts[2]

# Register credentials
wf_set_key(user = cds_uid, key = cds_token, overwrite = TRUE)
wf_set_url(cds_url)

# ERA5-Land hourly (recommended; aggregate to daily locally)
request <- list(
  dataset_short_name = "reanalysis-era5-land",
  product_type = "reanalysis",
  format = "netcdf",
  variable = c("2m_temperature", "total_precipitation"),
  year = "2022",                             # <-- adjust years
  month = sprintf("%02d", 1:12),
  day = sprintf("%02d", 1:31),
  time = sprintf("%02d:00", 0:23),           # all 24 hours
  area = area_for_request
)

# Uncomment to download (requires valid credentials & CDS availability):
# nc_file <- wf_request(
#   user = cds_uid,
#   request = request,
#   transfer = TRUE,
#   path = out_dir
# )
# message("Downloaded NetCDF: ", nc_file)

# If you've already downloaded: set path explicitly
nc_file <- "path_to_your_downloaded.nc"   # <-- replace after download
stopifnot(file.exists(nc_file))

# 5) LOAD ERA5, CROP & MASK TO BASIN
# ----------------------------------
era5_rast <- terra::rast(nc_file)
print(era5_rast)

# Set CRS if missing (ERA5-Land is in lon/lat WGS84)
if (is.na(terra::crs(era5_rast))) {
  terra::crs(era5_rast) <- "EPSG:4326"
}

# Fast crop to the vector extent, then mask
era5_crop <- terra::crop(era5_rast, nechako_vect, snap = "out")
era5_nechako <- terra::mask(era5_crop, nechako_vect)

print(era5_nechako)
plot(era5_nechako[[1]], main = "ERA5-Land (first layer)")
plot(nechako_vect, add = TRUE, border = "blue", lwd = 2)

# 6) OPTIONAL: DAILY AGGREGATION (mean for T2m, sum for precip)
# -------------------------------------------------------------
# Extract POSIXct time stamps from the raster
ts <- terra::time(era5_nechako)  # one per layer
stopifnot(!is.null(ts) && length(ts) == terra::nlyr(era5_nechako))
day_index <- as.Date(ts)

# Split variables by name (depends on CDS variable naming)
vnames <- names(era5_nechako)
is_t2m  <- grepl("2m_temperature", vnames, ignore.case = TRUE)
is_tp   <- grepl("total_precipitation", vnames, ignore.case = TRUE)

era5_t2m <- era5_nechako[[which(is_t2m)]]
era5_tp  <- era5_nechako[[which(is_tp)]]

# Daily mean temperature
t2m_daily <- terra::tapp(era5_t2m, index = day_index, fun = mean, na.rm = TRUE)
# Daily precipitation sum (hourly accumulation to daily total)
tp_daily  <- terra::tapp(era5_tp,  index = day_index, fun = sum,  na.rm = TRUE)

# Save daily products (GeoTIFFs)
t2m_tif <- file.path(out_dir, "ERA5L_T2m_daily_mean_Nechako.tif")
tp_tif  <- file.path(out_dir, "ERA5L_TP_daily_sum_Nechako.tif")
terra::writeRaster(t2m_daily, t2m_tif, overwrite = TRUE)
terra::writeRaster(tp_daily,  tp_tif,  overwrite = TRUE)
message("Daily T2m saved: ", t2m_tif)
message

# 5. VERIFY THE RESULT
# ---------------------
print(era5_nechako)
# Plot the first layer (e.g., Day 1 temperature) with the basin outline
# plot(era5_nechako[[1]], main = "ERA5-Land Data Cropped to Nechako Basin")
# lines(nechako_basin, col = "red", lwd = 2)
