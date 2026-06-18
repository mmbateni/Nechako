library(terra)
library(sf)
library(dplyr)

setwd("D:/Nechako_Drought/Nechako/Lakes")

# ============================================================================
# 1. DOWNLOAD CDED DEM DATA (25m resolution)
# ============================================================================
base_url <- "https://pub.data.gov.bc.ca/datasets/175624"
# URL directories: 93e, 93f (NO leading zero)
url_dirs <- c("93e", "93f")
# File prefixes: 093e, 093f (WITH leading zero)
file_prefixes <- c("093e", "093f")
download_dir <- "dem_data_cded"
dir.create(download_dir, showWarnings = FALSE)

# Generate URLs - mix of URL path (no zero) and filename (with zero)
urls <- c()
for (i in seq_along(url_dirs)) {
  url_dir <- url_dirs[i]      # e.g., "93e"
  file_pre <- file_prefixes[i] # e.g., "093e"
  
  for (block in sprintf("%02d", 1:16)) {
    for (side in c("e", "w")) {
      # URL: .../93e/093e01_e.dem.zip
      urls <- c(urls, file.path(base_url, url_dir, paste0(file_pre, block, "_", side, ".dem.zip")))
    }
  }
}

# Download and extract
for (url in urls) {
  zip_file <- file.path(download_dir, basename(url))
  dem_file <- file.path(download_dir, sub("\\.zip$", "", basename(url)))
  
  if (!file.exists(dem_file)) {
    if (!file.exists(zip_file)) {
      message(sprintf("Downloading %s...", basename(zip_file)))
      tryCatch({
        download.file(url, zip_file, mode = "wb", quiet = TRUE)
      }, error = function(e) message("Download failed: ", e$message))
    }
    
    if (file.exists(zip_file)) {
      message(sprintf("Extracting %s...", basename(zip_file)))
      unzip(zip_file, exdir = download_dir, overwrite = FALSE)
    }
  }
}

# ============================================================================
# 2. LOAD AND MERGE DEM TILES
# ============================================================================
dem_files <- list.files(download_dir, pattern = "\\.dem$", 
                        full.names = TRUE, ignore.case = TRUE)

if (length(dem_files) == 0) {
  stop("No .dem files found. Check download directory.")
}

message(sprintf("Loading %d DEM tiles...", length(dem_files)))
dem_list <- lapply(dem_files, rast)
dem <- do.call(merge, dem_list)

message(sprintf("DEM loaded. Resolution: %.1f m", res(dem)[1]))

# ============================================================================
# 3-8. Continue with existing code (polygon loading, clipping, hypsometry)
# ============================================================================
res_poly <- st_read("./FWA_LAKES_POLY/FWLKSPL_polygon.shp") |>
  filter(if_any(starts_with("GNIS_ID"), ~ .x == 16603)) |>
  st_union() |>
  vect()

res_poly <- project(res_poly, dem)
res_poly_buf <- buffer(res_poly, width = 50)

message("Clipping DEM to reservoir boundary...")
dem_clipped <- mask(crop(dem, res_poly_buf), res_poly)

message("Extracting hypsometry...")
hyps <- freq(round(dem_clipped)) |>
  as.data.frame() |>
  mutate(
    elevation = value,
    area_m2 = count * prod(res(dem_clipped)),
    area_km2 = area_m2 / 1e6,
    .keep = "unused"
  ) |>
  arrange(elevation) |>
  filter(!is.na(elevation))

dz <- diff(hyps$elevation)
mean_area <- (hyps$area_m2[-nrow(hyps)] + hyps$area_m2[-1]) / 2
segment_volume <- dz * mean_area

hyps <- hyps |>
  mutate(
    volume_m3 = c(0, cumsum(segment_volume)),
    volume_mcm = volume_m3 / 1e6
  )

cat("\n=== DEM Summary ===\n")
cat(sprintf("Elevation range: %.1f - %.1f m\n", 
            min(hyps$elevation), max(hyps$elevation)))
cat(sprintf("Total surface area: %.2f km²\n", sum(hyps$area_km2)))
cat(sprintf("Total volume: %.2f MCM\n", max(hyps$volume_mcm)))
cat(sprintf("DEM resolution: %.1f m\n", res(dem_clipped)[1]))

par(mfrow = c(1, 2))
plot(hyps$elevation, hyps$area_km2, type = "l", col = "blue", lwd = 2,
     xlab = "Elevation (m)", ylab = "Area (km²)", main = "Area-Elevation Curve")
grid()
plot(hyps$elevation, hyps$volume_mcm, type = "l", col = "red", lwd = 2,
     xlab = "Elevation (m)", ylab = "Volume (MCM)", main = "Volume-Elevation Curve")
grid()

head(hyps, 10)