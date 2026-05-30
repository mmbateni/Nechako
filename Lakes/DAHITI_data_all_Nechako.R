# ============================================================
# DAHITI: Download all lake/reservoir data in the Nechako basin (API v2)
# Improved with strict spatial filtering using local shapefile
# ============================================================
# 1. wse (Water Surface Elevation)
# What it means: The absolute height (elevation) of the lake or river's water surface at the exact time the satellite passed over. It is measured relative to a global reference "geoid" (usually EGM2008, which is a highly accurate mathematical model of global mean sea level).
# Unit: Meters (m)
# 2. wse_u (Water Surface Elevation Uncertainty)
# What it means: The statistical error margin (or standard deviation) of the wse measurement. Because radar/laser altimeters are measuring from space, things like atmospheric interference, waves, ice, or sensor noise can introduce slight inaccuracies. This value tells you the confidence interval of that specific water level reading
# Unit: Meters (m)
# 3. wsc (Water Storage Change / Volume Variation)
# What it means: The change in the volume of water stored in the lake or reservoir relative to a specific baseline average. Satellites cannot easily measure the absolute total volume of a lake because they cannot see the shape of the lake bottom. Instead, DAHITI calculates the change in volume by combining the Water Surface Elevation (wse) with the lake's Surface Area. It tells you how much water the lake has gained or lost.
# Unit: Typically Cubic Kilometers or Millions of Cubic Meters 
# 4. wsc_u (Water Storage Change Uncertainty)
# What it means: The statistical error margin of the wsc calculation. Because the volume variation is calculated using two different satellite measurements (water height + water area), the uncertainties of both of those measurements are combined into this final error margin.
# Unit: Same as wsc 
# ---- Packages ----
libs <- c("httr", "jsonlite", "dplyr", "purrr", "readr", "tibble", "fs", "glue", "sf")
for (p in libs) if (!requireNamespace(p, quietly = TRUE)) stop("Install package: ", p)

library(httr)
library(jsonlite)
library(dplyr)
library(purrr)
library(readr)
library(tibble)
library(fs)
library(glue)
library(sf)

# ---- CONFIG: put your DAHITI API key here ----
API_KEY <- Sys.getenv("DAHITI_API_KEY")
if (API_KEY == "") {
  API_KEY <- "C10E9AAF2751021683B291A721B1974C2EB30E025D3232828B4C07453E05E4C8"  # 
}
stopifnot(nchar(API_KEY) >= 30)  # lightweight sanity check

# ---- CONFIG: output directory & rate limiting ----
out_dir <- "D:/Nechako_Drought/Nechako/Lakes/dahiti_DE"
dir_create(out_dir)
# Be kind to DAHITI (v2 has ~250 requests/24h). Adjust if needed.
REQUEST_PAUSE_SEC <- 0.5

# ---- CONFIG: Nechako basin spatial filter ----
# 1. Read the provided shapefile
shp_path <- "D:/Nechako_Drought/Nechako/Spatial/nechakoBound_dissolve.shp"
if (!file.exists(shp_path)) {
  stop("Shapefile not found. Please check the path: ", shp_path)
}
basin_poly <- sf::st_read(shp_path, quiet = TRUE)

# --- Check the current Coordinate Reference System (CRS) ---
current_crs <- sf::st_crs(basin_poly)

if (is.na(current_crs)) {
  stop("ERROR: The shapefile does not have a defined coordinate system (missing .prj file).")
} else {
  message("Original Shapefile CRS: ", current_crs$name)
}

# 2. Transform to WGS84 (EPSG:4326) because DAHITI API expects Longitude/Latitude
if (is.na(current_crs$epsg) || current_crs$epsg != 4326) {
  message("Transforming shapefile to WGS84 (EPSG:4326)...")
  basin_poly_wgs84 <- sf::st_transform(basin_poly, 4326)
} else {
  message("Shapefile is already in WGS84 (EPSG:4326).")
  basin_poly_wgs84 <- basin_poly
}

# 3. Extract exact bounding box for the API query
bbox <- sf::st_bbox(basin_poly_wgs84)
NECHAKO_BBOX <- c(
  min_lon = as.numeric(bbox["xmin"]),
  max_lon = as.numeric(bbox["xmax"]),
  min_lat = as.numeric(bbox["ymin"]),
  max_lat = as.numeric(bbox["ymax"])
)

message("Bounding box calculated from shapefile: ", 
        paste(round(NECHAKO_BBOX, 3), collapse = ", "))

# ---- DAHITI API endpoints (v2) ----
base_v2 <- "https://dahiti.dgfi.tum.de/api/v2"
ep_list_targets          <- paste0(base_v2, "/list-targets/")
ep_get_target_info       <- paste0(base_v2, "/get-target-info/")
ep_download_water_level  <- paste0(base_v2, "/download-water-level/")
ep_download_surface_area <- paste0(base_v2, "/download-surface-area/")
ep_download_volume_var   <- paste0(base_v2, "/download-volume-variation/")
ep_download_hypsometry   <- paste0(base_v2, "/download-hypsometry/")
ep_download_bathymetry   <- paste0(base_v2, "/download-bathymetry/")
ep_download_water_occ    <- paste0(base_v2, "/download-water-occurrence-mask/")
ep_download_land_water   <- paste0(base_v2, "/download-land-water-mask/")

# ---- 1) Query targets in the region ----
query_args <- list(
  api_key = API_KEY,
  min_lon = NECHAKO_BBOX["min_lon"],
  max_lon = NECHAKO_BBOX["max_lon"],
  min_lat = NECHAKO_BBOX["min_lat"],
  max_lat = NECHAKO_BBOX["max_lat"]
)

resp <- httr::GET(ep_list_targets, query = query_args, httr::timeout(60))
stop_for_status(resp)
lst <- jsonlite::fromJSON(content(resp, "text", encoding = "UTF-8"), simplifyVector = FALSE)

if (is.null(lst$code) || lst$code != 200) {
  stop("DAHITI list-targets failed. Message: ", lst$message %||% "no details")
}

targets_raw <- lst$data
if (length(targets_raw) == 0) {
  stop("No DAHITI targets found in the bounding box.")
}

# Convert to a tibble; keep 'data_access' as list-col
targets_tbl <- tibble::tibble(
  dahiti_id   = vapply(targets_raw, `[[`, numeric(1), "dahiti_id"),
  target_name = vapply(targets_raw, `[[`, character(1), "target_name"),
  type        = vapply(targets_raw, `[[`, character(1), "type"),
  continent   = vapply(targets_raw, `[[`, character(1), "continent"),
  country     = vapply(targets_raw, `[[`, character(1), "country"),
  longitude   = vapply(targets_raw, `[[`, numeric(1), "longitude"),
  latitude    = vapply(targets_raw, `[[`, numeric(1), "latitude"),
  data_access = lapply(targets_raw, `[[`, "data_access")
)

# Keep Lakes, Reservoirs, AND Rivers
targets_to_keep <- targets_tbl %>% filter(type %in% c("Lake", "Reservoir", "River"))

if (nrow(targets_to_keep) == 0) {
  stop("No lakes, reservoirs, or rivers found by DAHITI in the bounding box.")
}

# ---- STRICT SPATIAL FILTERING ----
targets_sf <- sf::st_as_sf(targets_to_keep, coords = c("longitude", "latitude"), crs = 4326)
# Suppress the harmless planar warning
targets_in_basin_sf <- suppressWarnings(sf::st_filter(targets_sf, basin_poly_wgs84))

# Convert back to a standard dataframe
targets_to_keep <- targets_in_basin_sf %>% sf::st_drop_geometry()

if (nrow(targets_to_keep) == 0) {
  stop("Targets were found in the bounding box, but none fell strictly inside the Nechako basin polygon.")
}

message("Found ", nrow(targets_to_keep), " DAHITI targets strictly inside the Nechako basin.")

# ---- Helper: safely GET/POST JSON and binaries ----
dahiti_get_json <- function(url, query) {
  Sys.sleep(REQUEST_PAUSE_SEC)
  r <- httr::GET(url, query = query, httr::timeout(120))
  txt <- content(r, "text", encoding = "UTF-8")
  
  if (httr::http_error(r)) {
    return(list(ok = FALSE, err = paste("HTTP", r$status_code, "-", txt)))
  }
  
  out <- tryCatch(
    jsonlite::fromJSON(txt, simplifyVector = FALSE),
    error = function(e) e
  )
  
  if (inherits(out, "error")) {
    snippet <- substr(txt, 1, 80)
    snippet <- gsub("\n", " ", snippet) 
    return(list(ok = FALSE, err = paste("Invalid JSON received. Snippet:", snippet)))
  }
  
  if (is.list(out) && !is.null(out$code) && out$code != 200) {
    err_msg <- out$message %||% out$error %||% paste("API returned code:", out$code)
    return(list(ok = FALSE, err = err_msg))
  }
  
  if (is.list(out) && is.null(out$code) && !is.null(out$error)) {
    return(list(ok = FALSE, err = out$error))
  }
  
  list(ok = TRUE, data = out)
}

dahiti_get_binary <- function(url, query, path_out) {
  Sys.sleep(REQUEST_PAUSE_SEC)
  r <- httr::GET(url, query = query, httr::timeout(300), httr::write_disk(path_out, overwrite = TRUE))
  if (httr::http_error(r)) {
    list(ok = FALSE, err = paste("HTTP", r$status_code))
  } else {
    if (file_exists(path_out) && file_info(path_out)$size > 100) list(ok = TRUE, path = path_out) else list(ok = FALSE, err = "Empty/small file")
  }
}

# ---- Map DAHITI product flags to endpoints and file handling ----
product_map <- list(
  water_level_altimetry = list(flag = "water_level_altimetry", ep = ep_download_water_level,  mode = "json_table",  filename = "{id}_water_level.json"),
  surface_area          = list(flag = "surface_area",          ep = ep_download_surface_area, mode = "json_table",  filename = "{id}_surface_area.json"),
  volume_variation      = list(flag = "volume_variation",      ep = ep_download_volume_var,   mode = "json_table",  filename = "{id}_volume_variation.json"),
  hypsometry            = list(flag = "hypsometry",            ep = ep_download_hypsometry,   mode = "json_table",  filename = "{id}_hypsometry.json"),
  bathymetry            = list(flag = "bathymetry",            ep = ep_download_bathymetry,   mode = "binary_tif",  filename = "{id}_bathymetry.tif"),
  water_occurrence_mask = list(flag = "water_occurrence_mask", ep = ep_download_water_occ,    mode = "binary_tif",  filename = "{id}_water_occurrence_mask.tif"),
  land_water_mask       = list(flag = "land_water_mask",       ep = ep_download_land_water,   mode = "binary_tar",  filename = "{id}_land_water_mask.tar.gz")
)

# ---- 2) Download everything that is 'public' for each target ----
all_index <- list()
i <- 0L

for (row in seq_len(nrow(targets_to_keep))) {
  id   <- targets_to_keep$dahiti_id[row]
  name <- targets_to_keep$target_name[row]
  lon  <- targets_to_keep$longitude[row]   # ---> NEW: Extract Longitude
  lat  <- targets_to_keep$latitude[row]    # ---> NEW: Extract Latitude
  acc  <- targets_to_keep$data_access[[row]]
  
  message("\n=== Target ", id, " — ", name, " ===")
  
  # Create descriptive folder names
  # Replace spaces, commas, and special characters with underscores
  safe_name <- gsub("[^A-Za-z0-9]", "_", name)
  safe_name <- gsub("_+", "_", safe_name)      # Remove double underscores
  safe_name <- sub("^_|_$", "", safe_name)     # Remove leading/trailing underscores
  
  folder_name <- paste0(id, "_", safe_name)
  dir_tgt <- file.path(out_dir, folder_name)
  dir_create(dir_tgt)
  
  # Save basic metadata (get-target-info)
  info_q <- list(api_key = API_KEY, dahiti_id = id)
  info_res <- dahiti_get_json(ep_get_target_info, info_q)
  if (isTRUE(info_res$ok)) {
    write_json(info_res$data, file.path(dir_tgt, glue("{id}_info.json")), pretty = TRUE, auto_unbox = TRUE)
  } else {
    message("  (warn) get-target-info failed: ", info_res$err)
  }
  
  # Loop over products
  for (prod in names(product_map)) {
    meta <- product_map[[prod]]
    # Is it public for this target?
    is_public <- tryCatch(isTRUE(acc[[meta$flag]] == "public"), error = function(e) FALSE)
    if (!is_public) next
    
    fn  <- glue(meta$filename, id = id)
    path_out <- file.path(dir_tgt, fn)
    q <- list(api_key = API_KEY, dahiti_id = id)
    
    if (meta$mode == "json_table") {
      q$format <- "json" # ONLY request JSON format for the data download endpoints
      res <- dahiti_get_json(meta$ep, q)
      
      if (isTRUE(res$ok)) {
        write_json(res$data, path_out, pretty = TRUE, auto_unbox = TRUE)
        message("  ✓ ", prod, " → ", fn, " (JSON)")
        
        csv_path <- NA_character_
        
        # Attempt to create a CSV
        data_to_save <- NULL
        if (is.list(res$data) && !is.null(res$data$data)) {
          data_to_save <- res$data$data
        } else if (is.list(res$data) && is.null(res$data$code)) {
          data_to_save <- res$data
        }
        
        if (!is.null(data_to_save) && length(data_to_save) > 0) {
          df <- tryCatch(
            jsonlite::fromJSON(jsonlite::toJSON(data_to_save, auto_unbox = TRUE), simplifyVector = TRUE),
            error = function(e) NULL
          )
          if (!is.null(df) && is.data.frame(df)) {
            csv_out <- sub("\\.json$", ".csv", path_out)
            readr::write_csv(df, csv_out, na = "")
            message("    ↳ CSV table saved as ", basename(csv_out))
            csv_path <- csv_out
          }
        }
        
        i <- i + 1L
        all_index[[i]] <- tibble::tibble(
          dahiti_id   = id,
          target_name = name,
          longitude   = lon,     # ---> NEW: Add Longitude to index
          latitude    = lat,     # ---> NEW: Add Latitude to index
          product     = prod,
          file_json   = path_out,
          file_csv    = csv_path,
          file_bin    = NA_character_
        )
        
      } else {
        message("  (warn) ", prod, " failed: ", res$err)
      }
      
    } else if (meta$mode %in% c("binary_tif", "binary_tar")) {
      res <- dahiti_get_binary(meta$ep, q, path_out)
      if (isTRUE(res$ok)) {
        message("  ✓ ", prod, " → ", fn, " (", ifelse(meta$mode == "binary_tif", "GeoTIFF", "tar.gz"), ")")
        i <- i + 1L
        all_index[[i]] <- tibble::tibble(
          dahiti_id   = id,
          target_name = name,
          longitude   = lon,     # ---> NEW: Add Longitude to index
          latitude    = lat,     # ---> NEW: Add Latitude to index
          product     = prod,
          file_json   = NA_character_,
          file_csv    = NA_character_,
          file_bin    = path_out
        )
      } else {
        message("  (warn) ", prod, " failed: ", res$err)
      }
    }
  }
}

# ---- 3) Write a master index ----
index_tbl <- dplyr::bind_rows(all_index)
if (nrow(index_tbl) > 0) {
  readr::write_csv(index_tbl, file.path(out_dir, "download_index.csv"), na = "")
  message("\nMaster index written: ", file.path(out_dir, "download_index.csv"))
} else {
  message("\nNo public products were downloaded (check data_access flags and bbox).")
}