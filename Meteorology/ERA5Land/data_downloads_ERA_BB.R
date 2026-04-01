# ======= Nechako Basin ERA data Retrieval workflow using a Bounding Box =======
# Packages
library(terra)
library(sf)
library(dplyr)
library(ecmwfr)

# PREP EXTENTS FOR RASTER WORK AND ERA5 REQUEST

# CDS expects area vector: [North, West, South, East]
area_for_request <- c(
  as.numeric(56.5),
  as.numeric(-128.5),
  as.numeric(51.9),
  as.numeric(-122)
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
