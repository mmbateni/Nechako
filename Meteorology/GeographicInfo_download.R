library(bcdata)
library(sf)
library(ggplot2)
library(dplyr)
library(purrr)
library(terra)

# 1. DEFINE YOUR REGION OF INTEREST IN WGS84
nechako_bbox_wgs84 <- st_bbox(c(xmin = -128, ymin = 52.5, xmax = -122.0, ymax = 56.5), crs = 4326)
nechako_poly_wgs84 <- st_as_sfc(nechako_bbox_wgs84)

# 2. TRANSFORM TO BC ALBERS (EPSG:3005) FOR DATA DOWNLOAD
nechako_poly_3005 <- st_transform(nechako_poly_wgs84, 3005)
nechako_bbox_3005 <- st_bbox(nechako_poly_3005)

# 3. DEFINE DATA TO DOWNLOAD (COMMENT OUT LAYERS YOU DON'T NEED)
datasets_to_get <- list(
  watersheds = "freshwater-atlas-watersheds",
  dem_tiles = "dem-tile-index",  # We'll use a tile index for DEM
  rivers = "freshwater-atlas-rivers",
  lakes = "freshwater-atlas-lakes",
  obstructions = "freshwater-atlas-obstructions"
)


# 4. DOWNLOAD AND FILTER DATA

map_layers <- list()
cache_dir <- "bcdata_cache"
dir.create(cache_dir, showWarnings = FALSE)

# Function to download with caching
download_with_cache <- function(dataset_id, bbox = NULL, filter = NULL, is_raster = FALSE, roi_poly = NULL) {
  # Create a unique cache file name
  bbox_str <- if(!is.null(bbox)) paste(round(as.numeric(bbox), 3), collapse="_") else "full"
  cache_file <- file.path(cache_dir, paste0(dataset_id, "_", bbox_str, ".rds"))
  
  if(file.exists(cache_file)) {
    message("Loading cached: ", dataset_id)
    return(readRDS(cache_file))
  }
  
  message("Downloading: ", dataset_id)
  
  if(is_raster) {
    # 1. Get the tile index and find tiles within your Region of Interest (roi_poly)
    tile_index <- bcdc_get_data("dem-tile-index")
    tiles_in_roi <- tile_index[st_intersects(tile_index, roi_poly, sparse = FALSE), ]
    
    if(nrow(tiles_in_roi) == 0) {
      stop("No DEM tiles found in the specified bounding box.")
    }
    
    message("  Found ", nrow(tiles_in_roi), " tile(s) to download.")
    
    dem_list <- list()
    for(i in 1:nrow(tiles_in_roi)) {
      tile <- tiles_in_roi[i, ]
      message("  Downloading tile ", i, "/", nrow(tiles_in_roi), ": ", tile$TILE_NAME)
      
      tryCatch({
        tile_data <- bcdc_get_data(tile$TILE_NAME)
        dem_list[[i]] <- rast(tile_data)
      }, error = function(e) {
        message("  Failed to download tile: ", tile$TILE_NAME, " - ", e$message)
      })
    }
    
    # Merge all successfully downloaded tiles
    dem_list <- dem_list[!sapply(dem_list, is.null)]
    if(length(dem_list) > 0) {
      dem_merged <- do.call(merge, dem_list)
      # Crop the merged raster to your precise ROI polygon
      result <- crop(dem_merged, vect(roi_poly))
    } else {
      stop("No DEM tiles were successfully downloaded.")
    }
    
  } else {
    # For vector data: use bcdc_query_geodata with bbox filter
    query <- bcdc_query_geodata(record = dataset_id)
    
    # Apply bounding box filter if provided
    if(!is.null(bbox)) {
      # Convert bbox vector to sfc polygon for the query
      bbox_poly <- st_as_sfc(bbox)
      st_crs(bbox_poly) <- 3005
      query <- query %>% filter(bcdata::INTERSECTS(bbox_poly))
    }
    
    # Download the data
    result <- collect(query)
    
    # Apply additional attribute filters if specified
    if(!is.null(filter)) {
      result <- result %>% filter(!!rlang::parse_expr(filter))
    }
  }
  
  saveRDS(result, cache_file)
  return(result)
}
# Set a smaller chunk size before downloading the watersheds
options(bcdata.chunk_limit = 5000) # Default is 10000
# Download all data with caching
# Note: For the DEM, we pass roi_poly = nechako_poly_3005 to the function
if ("watersheds" %in% names(datasets_to_get)) {
  map_layers$watersheds <- download_with_cache(datasets_to_get$watersheds, bbox = nechako_bbox_3005)
}

if ("dem_tiles" %in% names(datasets_to_get)) {
  # For raster, pass the polygon for tile selection and cropping
  map_layers$dem <- download_with_cache("dem-tile-index", 
                                        bbox = nechako_bbox_3005, 
                                        is_raster = TRUE, 
                                        roi_poly = nechako_poly_3005)
}

if ("rivers" %in% names(datasets_to_get)) {
  map_layers$rivers <- download_with_cache(datasets_to_get$rivers, 
                                           bbox = nechako_bbox_3005, 
                                           filter = "STREAM_ORDER > 3")
}

if ("lakes" %in% names(datasets_to_get)) {
  map_layers$reservoir <- download_with_cache(datasets_to_get$lakes, 
                                              bbox = nechako_bbox_3005, 
                                              filter = 'GNIS_NAME == "Nechako Reservoir"')
}

if ("obstructions" %in% names(datasets_to_get)) {
  map_layers$obstructions <- download_with_cache(datasets_to_get$obstructions, 
                                                 bbox = nechako_bbox_3005)
}

message("All data downloaded and filtered.")

# 5. INSPECT WATERSHED DATA
if (!is.null(map_layers$watersheds)) {
  print("Sample watershed identifiers:")
  print(head(unique(map_layers$watersheds$WATERSHED_KEY)))
  print(head(unique(map_layers$watersheds$GNIS_NAME_1)))
}

# 6. CREATE AND CUSTOMIZE THE MAP
final_map <- ggplot()

# Hillshade/DEM Base Layer
if (!is.null(map_layers$dem)) {
  slope <- terrain(map_layers$dem, "slope", unit = "radians")
  aspect <- terrain(map_layers$dem, "aspect", unit = "radians")
  hillshade <- shade(slope, aspect, angle = 45, direction = 315)
  
  final_map <- final_map +
    geom_spatraster(data = hillshade) +
    scale_fill_distiller(palette = "Greys", direction = 1, na.value = NA)
}

# Watershed Sub-basins
if (!is.null(map_layers$watersheds)) {
  final_map <- final_map +
    geom_sf(data = map_layers$watersheds, fill = NA, color = "grey50", linewidth = 0.2)
}

# Main Watershed Boundary
if (!is.null(map_layers$watersheds)) {
  main_ws_key <- "NE-001" # <-- CHANGE THIS after inspecting data in Step 5
  main_watershed <- map_layers$watersheds %>% filter(WATERSHED_KEY == main_ws_key)
  if (nrow(main_watershed) > 0) {
    final_map <- final_map +
      geom_sf(data = main_watershed, fill = NA, color = "black", linewidth = 1.2)
  } else {
    warning("Main watershed key not found. Skipping that layer.")
  }
}

# Major Rivers
if (!is.null(map_layers$rivers)) {
  final_map <- final_map +
    geom_sf(data = map_layers$rivers, color = "steelblue", linewidth = 0.3)
}

# Reservoir
if (!is.null(map_layers$reservoir)) {
  final_map <- final_map +
    geom_sf(data = map_layers$reservoir, fill = "lightblue", color = "darkblue", linewidth = 0.4)
}

# Obstructions
if (!is.null(map_layers$obstructions)) {
  final_map <- final_map +
    geom_sf(data = map_layers$obstructions, color = "red", size = 2, shape = 17)
}

# 7. APPLY FINAL THEME AND DISPLAY
final_map <- final_map +
  theme_void() +
  labs(title = "Nechako Watershed")

print(final_map)