library(terra)
library(ecmwfr)

library(sf)
library(dplyr)

# 1) Read WATERSHED GROUPS (to discover correct group codes)
wfs_groups <- paste0(
  "https://openmaps.gov.bc.ca/geo/pub/WHSE_BASEMAPPING.FWA_WATERSHED_GROUPS_POLY/ows?",
  "service=WFS&version=2.0.0&request=GetFeature&",
  "typeNames=pub:WHSE_BASEMAPPING.FWA_WATERSHED_GROUPS_POLY&",
  "srsName=EPSG:3005&outputFormat=application/json"
)
fwa_groups <- st_read(wfs_groups, quiet = TRUE)
stopifnot(!is.na(st_crs(fwa_groups)))  # 3005

# Find all watershed groups whose names mention your study area keywords
nechako_groups <- fwa_groups %>%
  filter(
    grepl("Nechako",  WATERSHED_GROUP_NAME, ignore.case = TRUE) |
      grepl("Stuart",   WATERSHED_GROUP_NAME, ignore.case = TRUE) |
      grepl("Stellako", WATERSHED_GROUP_NAME, ignore.case = TRUE) |
      grepl("Endako",   WATERSHED_GROUP_NAME, ignore.case = TRUE) |
      grepl("Chilako",  WATERSHED_GROUP_NAME, ignore.case = TRUE) |
      grepl("Nautley",  WATERSHED_GROUP_NAME, ignore.case = TRUE) |
      grepl("Cheslatta",WATERSHED_GROUP_NAME, ignore.case = TRUE)
  )

# Capture the group codes (these are consistent across FWA layers)
codes <- sort(unique(nechako_groups$WATERSHED_GROUP_CODE))
print(codes)   # sanity check: expect codes for Nechako + tributaries

# Build a mask geometry for spatial QC (optional but useful)
nechako_mask <- nechako_groups %>% summarise() %>% st_make_valid()

# 2) Read **ASSESSMENT WATERSHEDS** (20k mesoscale polygons, target 2k–10k ha)
wfs_aw <- paste0(
  "https://openmaps.gov.bc.ca/geo/pub/WHSE_BASEMAPPING.FWA_ASSESSMENT_WATERSHEDS_POLY/ows?",
  "service=WFS&version=2.0.0&request=GetFeature&",
  "typeNames=pub:WHSE_BASEMAPPING.FWA_ASSESSMENT_WATERSHEDS_POLY&",
  "srsName=EPSG:3005&outputFormat=application/json"
)
aw20k <- st_read(wfs_aw, quiet = TRUE)
stopifnot(!is.na(st_crs(aw20k)))  # 3005

# Sanity: AREA_HA should cluster around thousands (not single digits)
summary(aw20k$AREA_HA)

# 3) Filter by watershed group code (attribute filter) and intersect with the mask (spatial)
nechako_aw <- aw20k %>%
  filter(WATERSHED_GROUP_CODE %in% codes) %>%
  st_filter(nechako_mask) %>%     # keeps only polygons inside the group mask
  st_make_valid()

cat("Assessment watersheds selected:", nrow(nechako_aw), "\n")
# Expect **hundreds**, not 28; and AREA_HA median should be thousands.

# 4) Dissolve to one basin outline, then transform to WGS84
nechako_basin <- nechako_aw %>% summarise() %>% st_make_valid()
nechako_basin_ll <- st_transform(nechako_basin, 4326)

plot(st_geometry(nechako_basin_ll), col = "lightblue", border = "black", lwd = 1.5)
title("Nechako Basin (Assessment Watersheds → dissolved, WGS84)")




# 5) Write outputs safely
#    (A) Shapefile: write to a directory (dsn) and use a layer name
out_dir   <- "out"
out_layer <- "Nechako_FWA_Group_EPSG4326"
dir.create(out_dir, showWarnings = FALSE)

# Remove old files if present (optional)
unlink(file.path(out_dir, paste0(out_layer, c(".shp",".shx",".dbf",".prj",".cpg"))))

st_write(nechako_basin_ll, dsn = "out", layer = "Nechako_NRBlike_EPSG4326",
         driver = "ESRI Shapefile", delete_layer = TRUE)


# 2. INSPECT THE SHAPEFILE ON CONSOLE & PLOT
# -------------------------------------------

# -------------------------------------------
# 6) QUICK QA: PRINT SUMMARY, BBOX, AREA, AND PLOT
# -------------------------------------------

message("\n=== Nechako FWA Watershed Group: Summary ===")

# A) Basic sf summary (geometry type, CRS)
print(nechako_basin_ll)
message("CRS (EPSG): ", st_crs(nechako_basin_ll)$epsg)

# B) Geographic bounding box (WGS84) — useful for ERA5-Land requests
bbox_ll <- st_bbox(nechako_basin_ll)  # returns xmin, ymin, xmax, ymax in lon/lat
print(bbox_ll)

# Optional: pad the ERA5 bbox a bit (e.g., 0.05 degrees) to ensure full coverage
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

# C) Extent via terra (if you prefer terra::ext)
#    Note: ext() expects a Spat* object; convert sf -> SpatVector first.
nechako_vect <- terra::vect(nechako_basin_ll)
era5_ext <- terra::ext(nechako_vect)
print(era5_ext)

# D) Area (km^2) computed in BC Albers (EPSG:3005) for accuracy
nechako_basin_bc <- st_transform(nechako_basin_ll, 3005)
area_m2 <- as.numeric(st_area(nechako_basin_bc))
area_km2 <- area_m2 / 1e6
message("Area (BC Albers): ", format(round(area_km2, 1), big.mark = ","), " km^2")

# E) Validity check
is_valid <- sf::st_is_valid(nechako_basin_ll)
message("Geometry valid? ", is_valid)

# F) Quick plot
plot(nechako_basin_ll, border = "blue", lwd = 2,
     main = "Nechako Watershed Group (FWA) – EPSG:4326")

# Optional: also save a small PNG preview
png_filename <- file.path("out", "Nechako_FWA_Group_preview.png")
png(png_filename, width = 1200, height = 900, res = 120)
plot(nechako_basin_ll, border = "blue", lwd = 2,
     main = "Nechako Watershed Group (FWA) – EPSG:4326")
dev.off()
message("Preview map saved: ", png_filename)

# 3. CONFIGURE AND DOWNLOAD ERA5-LAND DAILY DATA
# ----------------------------------------------
# Use the basin's extent for the download area to minimize file size
basin_ext <- ext(nechako_basin)
# Format as North, West, South, East for the CDS request
area_for_request <- c(era5_bbox_padded["ymax"], era5_bbox_padded["xmin"],
                      era5_bbox_padded["ymin"], era5_bbox_padded["xmax"])

# --1) Use the cross-platform encrypted file backend to avoid Windows wincred issues
options(keyring_backend = "file")

# --2) Read the .cdsapirc from the current working directory
conf_path <- file.path(getwd(), ".cdsapirc")
stopifnot(file.exists(conf_path))
lines    <- readLines(conf_path)
key_line <- grep("^\\s*key\\s*:", lines, value = TRUE)
token    <- trimws(sub("^\\s*key\\s*:\\s*", "", key_line))

# --3) Register the token with ecmwfr and verify
wf_set_key(key = token)

request <- list(
  format = "netcdf",
  variable = c("2m_temperature", "total_precipitation"),
  product_type = "reanalysis",
  year = "2022", # Example year; modify your list of years as needed
  month = sprintf("%02d", 1:12),
  day = sprintf("%02d", 1:31),
  time = "00:00",
  statistic = "daily_mean", # Key parameter for daily data
  area = area_for_request,   # Use the dynamic basin extent
  dataset_short_name = "reanalysis-era5-land"
)

wf_request(request, transfer = TRUE, path = "out") # Uncomment and run with your credentials
nc_file <- "path_to_your_downloaded.nc" # Define this after download

# 4. LOAD ERA5 DATA, CROP, AND MASK TO BASIN
# -------------------------------------------
# Load the downloaded NetCDF as a SpatRaster
era5_rast <- rast(nc_file)
print(era5_rast)
crs(era5_rast)        # check projection
time(era5_rast)       # if available, shows time stamps (terra >= 1.6)
# If crs(era5_rast) is empty but the data are lon/lat, set it:
crs(era5_rast) <- "EPSG:4326"
# Crop to the basin extent first (faster)
era5_crop <- crop(era5_rast, nechako_vect, snap = "out")
# Crop and mask the raster to the basin polygon in one step
era5_nechako <- crop(era5_rast, nechako_basin, mask = TRUE)
# Verify
print(era5_nechako)
plot(era5_nechako[[1]], main = "ERA5 layer (first time slice)")
plot(nechako_vect, add = TRUE, border = "blue", lwd = 2)


# 5. VERIFY THE RESULT
# ---------------------
# Print information about the cropped raster
print(era5_nechako)

# Plot the first layer (e.g., Day 1 temperature) with the basin outline
# plot(era5_nechako[[1]], main = "ERA5-Land Data Cropped to Nechako Basin")
# lines(nechako_basin, col = "red", lwd = 2)