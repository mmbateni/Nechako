# ============================================================================
# NECHAKO RESERVOIR — HYPSOMETRIC TABLE DERIVATION
# Independent triple-source pipeline
# PIPELINE A — 3D-LAKES  (Huang et al., Sci. Data 12:1625, 2025)
# PIPELINE B — DAHITI    (DGFI-TUM, Schwatke et al.)
# PIPELINE C — BC CDED DEM (Local derivation from 25m DEM + FWA Polygons)
# ============================================================================

rm(list = ls())
setwd("D:/Nechako_Drought/Nechako/Lakes")

# ── PACKAGES ──────────────────────────────────────────────────────────────────
packages <- c("tidyverse", "lubridate", "httr", "jsonlite", "ggplot2", "terra", "sf")
for (p in packages) {
  if (!require(p, character.only = TRUE, quietly = TRUE)) {
    install.packages(p, quiet = TRUE)
    library(p, character.only = TRUE)
  }
}

# ============================================================================
#   # ── SECTION 1 ── TARGET CONFIGURATION
# ============================================================================
TARGET_LEVELS_M <- c(842.0, 844.0, 846.0, 848.0, 849.0,
                     850.0, 851.0, 852.0, 853.0, 854.0, 855.0)
NECHAKO_LAT <- 54.2493
NECHAKO_LON <- -125.7896

# ============================================================================
#   # ── SECTION 2 ── 3D-LAKES SETTINGS
# ============================================================================
NECHAKO_HYLAK_ID <- 6525L
LAKES_ST_FILE <- "3D-LAKES-ST.csv"
target_csv_name <- sprintf("%d_L1.csv", NECHAKO_HYLAK_ID)
AE_EXTRACT_DIR  <- "3DLakes/L1 A-E product"

LAKES_AE_CSV <- NA_character_
for (d in c(AE_EXTRACT_DIR, "3DLakes", ".")) {
  p <- file.path(d, target_csv_name)
  if (file.exists(p)) { LAKES_AE_CSV <- p; break }
}

if (is.na(LAKES_AE_CSV)) {
  zip_files  <- list.files("3DLakes", pattern = "\\.zip$",
                           full.names = TRUE, recursive = FALSE)
  for (zf in zip_files) {
    cat(sprintf("  [3D-LAKES] Checking ZIP: %s\n", basename(zf)))
    zip_contents  <- tryCatch(unzip(zf, list = TRUE)$Name,
                              error = function(e) character(0))
    hit  <- grep(paste0("(^|/)", target_csv_name, "$"), zip_contents, value = TRUE)
    if (length(hit) > 0L) {
      cat(sprintf("  [3D-LAKES] Found '%s' in ZIP. Extracting...\n", target_csv_name))
      dir.create(AE_EXTRACT_DIR, showWarnings = FALSE, recursive = TRUE)
      tryCatch({
        unzip(zf, files = hit[1L], exdir = AE_EXTRACT_DIR, junkpaths = TRUE)
        out_path  <- file.path(AE_EXTRACT_DIR, target_csv_name)
        if (file.exists(out_path)) LAKES_AE_CSV  <- out_path
      }, error = function(e) message("  ZIP extraction failed: ", conditionMessage(e)))
    }
    if (!is.na(LAKES_AE_CSV)) break
  }
}

if (!is.na(LAKES_AE_CSV)) {
  cat(sprintf("  [3D-LAKES] A-E file ready: %s\n", LAKES_AE_CSV))
} else {
  cat(sprintf("  [3D-LAKES] '%s' not found.\n", target_csv_name))
}

LAKES_AE_JSON <- NA_character_
LAKES_AE_ZIP  <- NA_character_

# ============================================================================
#   # ── SECTION 3 ── DAHITI SETTINGS
# ============================================================================
DAHITI_API_KEY <- Sys.getenv("DAHITI_API_KEY")
if (nchar(DAHITI_API_KEY) < 30L) {
  DAHITI_API_KEY <- "C10E9AAF2751021683B291A721B1974C2EB30E025D3232828B4C07453E05E4C8"
}
DAHITI_BASE_URL <- "https://dahiti.dgfi.tum.de/api/v2"
DAHITI_BBOX_PAD <- 1.5
RESERVOIR_ELEV_THRESHOLD_M <- 840.0
DAHITI_PAUSE_SEC <- 0.5

# ============================================================================
#   # ── SECTION 4 ── DEM SETTINGS (PIPELINE C)
# ============================================================================
DEM_GNIS_ID <- 16603
FWA_POLY_PATH <- "./FWA_LAKES_POLY/FWLKSPL_polygon.shp"
DEM_DOWNLOAD_DIR <- "dem_data_cded"

# ============================================================================
#   # ── SECTION 5 ── HELPER FUNCTIONS
# ============================================================================
dahiti_get <- function(endpoint, args, max_tries = 3L) {
  args[["api_key"]] <- DAHITI_API_KEY
  url <- paste0(DAHITI_BASE_URL, "/", endpoint, "/")
  for (i in seq_len(max_tries)) {
    Sys.sleep(DAHITI_PAUSE_SEC)
    resp <- tryCatch(
      httr::GET(url, query = args, httr::timeout(60L)),
      error = function(e) {
        message("  [attempt ", i, "] Network error: ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(resp) && httr::status_code(resp) == 200L) {
      txt <- httr::content(resp, "text", encoding = "UTF-8")
      return(jsonlite::fromJSON(txt, simplifyVector = FALSE))
    }
    if (!is.null(resp)) {
      message("  [attempt ", i, "] HTTP ", httr::status_code(resp), " from ", endpoint)
    }
    if (i < max_tries) Sys.sleep(2)
  }
  stop("DAHITI GET '", endpoint, "' failed after ", max_tries, " attempts.")
}

ae_to_ve <- function(elevation_m, area_m2) {
  keep <- !is.na(elevation_m) & !is.na(area_m2) & (area_m2 >= 0)
  h <- elevation_m[keep]; a <- area_m2[keep]
  ord <- order(h); h <- h[ord]; a <- a[ord]; n <- length(h)
  cum_vol_m3 <- numeric(n)
  for (i in seq_len(n - 1L)) {
    cum_vol_m3[i + 1L] <- cum_vol_m3[i] +
      0.5 * (a[i] + a[i + 1L]) * (h[i + 1L] - h[i])
  }
  data.frame(elevation_m = h, area_m2 = a, storage_Mm3 = cum_vol_m3 / 1e6)
}

interp_storage <- function(ve_df, target_levels_m) {
  data.frame(
    level_m     = target_levels_m,
    storage_Mm3 = approx(x = ve_df$elevation_m, y = ve_df$storage_Mm3,
                         xout = target_levels_m, method = "linear",
                         rule = 1L)$y
  )
}

dahiti_area_km2 <- function(h, x_wp, y_min, y_scale, z) {
  x_wp * (pmax(h - y_min, 0) / y_scale)^z
}

dahiti_storage_Mm3 <- function(h, x_wp, y_min, y_scale, z) {
  x_wp * y_scale / (z + 1) * (pmax(h - y_min, 0) / y_scale)^(z + 1)
}

pluck_num <- function(x, key) {
  v <- x[[key]]
  if (is.null(v)) return(NA_real_)
  as.numeric(v)
}

print_header <- function(title) {
  cat("\n\n")
  cat(strrep("=", 68), "\n")
  cat(title, "\n")
  cat(strrep("=", 68), "\n")
}

# ============================================================================
#   # ── SECTION 6 ── PIPELINE A: 3D-LAKES
# ============================================================================
print_header("PIPELINE A — 3D-LAKES  (Huang et al. 2025)")

hylak_id_used <- NECHAKO_HYLAK_ID
validated <- FALSE

if (file.exists(LAKES_ST_FILE)) {
  st <- tryCatch(
    read_csv(LAKES_ST_FILE,
             col_select     = c(Hylak_id, Pour_long, Pour_lat, Grand_id),
             show_col_types = FALSE),
    error = function(e) { warning("Could not read ", LAKES_ST_FILE); NULL }
  )
  
  if (!is.null(st)) {
    candidates <- st %>%
      mutate(dist_deg = sqrt((Pour_long - NECHAKO_LON)^2 + (Pour_lat - NECHAKO_LAT)^2)) %>%
      filter(dist_deg < 1.5, Grand_id > 0) %>%
      arrange(dist_deg)
    
    cat(sprintf("  Found %d reservoir candidates (Grand_id > 0) within 1.5 degrees.\n", nrow(candidates)))
    
    for (i in seq_len(min(10L, nrow(candidates)))) {
      cand_id <- candidates$Hylak_id[i]
      cand_csv_name <- sprintf("%d_L1.csv", cand_id)
      
      cand_csv_path <- NA_character_
      for (d in c(AE_EXTRACT_DIR, "3DLakes", ".")) {
        p <- file.path(d, cand_csv_name)
        if (file.exists(p)) { cand_csv_path <- p; break }
      }
      
      if (!is.na(cand_csv_path)) {
        ae_check <- tryCatch(
          read_csv(cand_csv_path, show_col_types = FALSE),
          error = function(e) NULL
        )
        
        if (!is.null(ae_check)) {
          nm <- names(ae_check)
          area_col <- nm[grepl("area", nm, ignore.case = TRUE)][1]
          elev_col <- nm[grepl("elev", nm, ignore.case = TRUE)][1]
          
          if (!is.na(area_col) && !is.na(elev_col)) {
            max_elev <- max(ae_check[[elev_col]], na.rm = TRUE)
            max_area_km2 <- max(ae_check[[area_col]], na.rm = TRUE) / 1e6
            
            if (max_elev > 830 && max_area_km2 > 100) {
              hylak_id_used <- cand_id
              LAKES_AE_CSV <- cand_csv_path
              validated <- TRUE
              cat(sprintf("  ✓ Validated Nechako Reservoir candidate: Hylak_id %d\n", hylak_id_used))
              break
            }
          }
        }
      }
    }
  }
}

ae_raw <- NULL
if (!is.na(LAKES_AE_CSV) && file.exists(LAKES_AE_CSV)) {
  cat(sprintf("\n  Reading A-E data from: %s\n", LAKES_AE_CSV))
  ae_raw <- tryCatch(
    read_csv(LAKES_AE_CSV, show_col_types = FALSE),
    error = function(e) { warning("Cannot read ", LAKES_AE_CSV); NULL }
  )
} else if (!is.na(LAKES_AE_JSON) && file.exists(LAKES_AE_JSON)) {
  cat(sprintf("\n  Reading A-E data from JSON (key '%d'): %s\n", hylak_id_used, LAKES_AE_JSON))
  json_all <- tryCatch(
    jsonlite::fromJSON(LAKES_AE_JSON, simplifyDataFrame = FALSE),
    error = function(e) { warning("Cannot parse JSON: ", conditionMessage(e)); NULL }
  )
  if (!is.null(json_all)) {
    key <- as.character(hylak_id_used)
    if (key %in% names(json_all)) {
      entry   <- json_all[[key]]
      ae_raw  <- data.frame(
        `Area (m2)`      = unlist(entry[["Area (m2)"]]),
        `Elevation (m)`  = unlist(entry[["Elevation (m)"]]),
        check.names      = FALSE
      )
    }
  }
}

if (!is.null(ae_raw)) {
  nm        <- names(ae_raw)
  area_col  <- nm[grepl("area", nm, ignore.case = TRUE)][1]
  elev_col  <- nm[grepl("elev", nm, ignore.case = TRUE)][1]
  if (is.na(area_col) || is.na(elev_col)) {
    warning("Unexpected column names in A-E file: ", paste(nm, collapse = ", "))
    ae_raw <- NULL
  } else {
    ae_raw <- ae_raw %>%
      rename(area_m2 = all_of(area_col), elevation_m = all_of(elev_col)) %>%
      filter(!is.na(area_m2), !is.na(elevation_m), area_m2 >= 0)
  }
}

results_3dlakes <- NULL
ve_3dlakes      <- NULL
if (!is.null(ae_raw) && nrow(ae_raw) >= 3L) {
  cat(sprintf("\n  A-E data loaded: %d points, elevation %.2f – %.2f m\n",
              nrow(ae_raw), min(ae_raw$elevation_m), max(ae_raw$elevation_m)))
  
  ve_3dlakes      <- ae_to_ve(ae_raw$elevation_m, ae_raw$area_m2)
  results_3dlakes <- interp_storage(ve_3dlakes, TARGET_LEVELS_M) %>%
    rename(storage_3dlakes_Mm3 = storage_Mm3)
  
  cat("\n  3D-LAKES — Storage at target levels:\n")
  print(results_3dlakes, row.names = FALSE, digits = 4)
}

# ============================================================================
#   # ── SECTION 7 ── PIPELINE B: DAHITI
# ============================================================================
print_header("PIPELINE B — DAHITI  (DGFI-TUM)")

results_dahiti <- NULL
if (nchar(DAHITI_API_KEY) < 30L) {
  cat("  ➜  DAHITI_API_KEY not set. Pipeline B skipped.\n")
} else {
  cat(sprintf(
    "\n  Querying list-targets in bbox (±%.1f° around reservoir centroid)...\n",
    DAHITI_BBOX_PAD))
  bbox_args <- list(
    min_lon = NECHAKO_LON - DAHITI_BBOX_PAD,
    max_lon = NECHAKO_LON + DAHITI_BBOX_PAD,
    min_lat = NECHAKO_LAT - DAHITI_BBOX_PAD,
    max_lat = NECHAKO_LAT + DAHITI_BBOX_PAD
  )
  tgt_resp <- tryCatch(
    dahiti_get("list-targets", bbox_args),
    error = function(e) {
      message("  list-targets failed: ", conditionMessage(e)); NULL
    }
  )
  
  if (is.null(tgt_resp) || length(tgt_resp$data) == 0L) {
    cat("  No DAHITI targets found in bounding box. Pipeline B skipped.\n")
  } else {
    tgt_list <- tgt_resp$data
    tgt_summary  <- data.frame(
      dahiti_id   = vapply(tgt_list, function(x) as.integer(x$dahiti_id), integer(1)),
      target_name = vapply(tgt_list, function(x) x$target_name %||% " ", character(1)),
      type        = vapply(tgt_list, function(x) x$type       %||% " ", character(1)),
      longitude   = vapply(tgt_list, function(x) as.numeric(x$longitude), numeric(1)),
      latitude    = vapply(tgt_list, function(x) as.numeric(x$latitude),  numeric(1)),
      has_hypso   = vapply(tgt_list, function(x) {
        isTRUE(tryCatch(x$data_access$hypsometry == "public",
                        error = function(e) FALSE))
      }, logical(1)),
      stringsAsFactors = FALSE
    )
    
    tgt_summary$dist_deg  <- sqrt(
      (tgt_summary$longitude - NECHAKO_LON)^2 +
        (tgt_summary$latitude  - NECHAKO_LAT)^2
    )
    tgt_summary  <- tgt_summary[order(tgt_summary$dist_deg), ]
    
    cat(sprintf("  Found %d targets in bbox; %d with public hypsometry:\n",
                nrow(tgt_summary), sum(tgt_summary$has_hypso)))
    
    hypso_candidates <- tgt_summary[tgt_summary$has_hypso, ]
    
    if (nrow(hypso_candidates) == 0L) {
      cat("\n  No public hypsometry products available. Pipeline B skipped.\n")
    } else {
      cat(sprintf("\n  Testing %d candidate station(s)...\n", nrow(hypso_candidates)))
      
      for (j in seq_len(nrow(hypso_candidates))) {
        dahiti_id  <- hypso_candidates$dahiti_id[j]
        stn_name   <- hypso_candidates$target_name[j]
        
        cat(sprintf("\n  ── Candidate %d/%d: ID %d (%s) ──\n",
                    j, nrow(hypso_candidates), dahiti_id, stn_name))
        
        hypso  <- tryCatch(
          dahiti_get("download-hypsometry",
                     list(dahiti_id = dahiti_id, format = "json")),
          error = function(e) { message("  download-hypsometry failed: ", conditionMessage(e)); NULL }
        )
        
        hp_data  <- NULL
        if (!is.null(hypso)) {
          if (!is.null(hypso$data$y_min)) hp_data  <- hypso$data
          else if (!is.null(hypso$data$data$y_min)) hp_data  <- hypso$data$data
          else if (!is.null(hypso$y_min)) hp_data  <- hypso
        }
        
        if (is.null(hp_data)) next
        
        x_wp    <- pluck_num(hp_data, "x_wp")
        y_min   <- pluck_num(hp_data, "y_min")
        y_scl   <- pluck_num(hp_data, "y_scale")
        z_exp   <- pluck_num(hp_data, "z")
        
        if (any(is.na(c(x_wp, y_min, y_scl, z_exp)))) next
        if (y_min < RESERVOIR_ELEV_THRESHOLD_M) {
          cat(sprintf("  ✗  y_min = %.1f m < threshold. Skipping river station.\n", y_min))
          next
        }
        
        cat(sprintf("  ✓  Reservoir-level station accepted (y_min = %.1f m).\n", y_min))
        
        area_km2     <- dahiti_area_km2(   TARGET_LEVELS_M, x_wp, y_min, y_scl, z_exp)
        storage_Mm3  <- dahiti_storage_Mm3(TARGET_LEVELS_M, x_wp, y_min, y_scl, z_exp)
        
        results_dahiti  <- data.frame(
          level_m            = TARGET_LEVELS_M,
          area_dahiti_km2    = round(area_km2,    3),
          storage_dahiti_Mm3 = round(storage_Mm3, 2)
        )
        
        cat("\n  DAHITI — Storage at target levels:\n")
        print(results_dahiti, row.names = FALSE, digits = 4)
        break
      }
    }
  }
}

# ============================================================================
#   # ── SECTION 8 ── PIPELINE C: BC CDED DEM (Local Derivation)
# ============================================================================
print_header("PIPELINE C — BC CDED DEM (Local Derivation)")
results_dem <- NULL; ve_dem <- NULL

if (!file.exists(FWA_POLY_PATH)) {
  cat(sprintf("  ➜  FWA polygon not found at %s. Pipeline C skipped.\n", FWA_POLY_PATH))
} else {
  # 8a. Download / Extract DEM tiles
  base_url <- "https://pub.data.gov.bc.ca/datasets/175624"
  url_dirs <- c("93e", "93f"); file_prefixes <- c("093e", "093f")
  dir.create(DEM_DOWNLOAD_DIR, showWarnings = FALSE, recursive = TRUE)
  
  urls <- c()
  for (i in seq_along(url_dirs)) {
    for (block in sprintf("%02d", 1:16)) {
      for (side in c("e", "w")) {
        urls <- c(urls, file.path(base_url, url_dirs[i], paste0(file_prefixes[i], block, "_", side, ".dem.zip")))
      }
    }
  }
  
  cat("  Checking / downloading CDED DEM tiles...\n")
  for (url in urls) {
    zip_file <- file.path(DEM_DOWNLOAD_DIR, basename(url))
    dem_file <- file.path(DEM_DOWNLOAD_DIR, sub("\\.zip$", "", basename(url)))
    if (!file.exists(dem_file)) {
      if (!file.exists(zip_file)) tryCatch(download.file(url, zip_file, mode = "wb", quiet = TRUE), error = function(e) NULL)
      if (file.exists(zip_file)) tryCatch(unzip(zip_file, exdir = DEM_DOWNLOAD_DIR, overwrite = FALSE), error = function(e) NULL)
    }
  }
  
  # 8b. Load and merge DEM
  dem_files <- list.files(DEM_DOWNLOAD_DIR, pattern = "\\.dem$", full.names = TRUE, ignore.case = TRUE)
  if (length(dem_files) == 0) {
    cat("  ➜  No .dem files found. Pipeline C skipped.\n")
  } else {
    cat(sprintf("  Loading and merging %d DEM tiles...\n", length(dem_files)))
    dem_list <- lapply(dem_files, terra::rast)
    dem <- do.call(terra::merge, dem_list)
    
    # 8c. Load reservoir polygon and buffer in metric CRS
    res_poly_sf <- sf::st_read(FWA_POLY_PATH, quiet = TRUE) |>
      dplyr::filter(dplyr::if_any(dplyr::starts_with("GNIS_ID"), ~ .x == DEM_GNIS_ID)) |>
      sf::st_union()
    
    if (length(res_poly_sf) == 0 || sf::st_is_empty(res_poly_sf)) {
      cat("  ➜  No polygon features matched GNIS_ID. Pipeline C skipped.\n")
    } else {
      res_poly_buf_sf <- sf::st_buffer(res_poly_sf, dist = 50) # 50m buffer
      res_poly_buf_terra <- terra::vect(res_poly_buf_sf)
      res_poly_terra <- terra::vect(res_poly_sf)
      
      # Project to DEM CRS (geographic) for cropping
      res_poly_buf_geo <- terra::project(res_poly_buf_terra, dem)
      res_poly_geo <- terra::project(res_poly_terra, dem)
      
      cat("  Clipping DEM to reservoir boundary...\n")
      dem_clipped_geo <- terra::mask(terra::crop(dem, res_poly_buf_geo), res_poly_geo)
      
      # Project clipped DEM to metric (BC Albers) for accurate area calculation
      cat("  Projecting clipped DEM to BC Albers (EPSG:3005)...\n")
      dem_clipped_m <- terra::project(dem_clipped_geo, "EPSG:3005", method = "near")
      
      # 8d. Extract hypsometry
      cat("  Extracting hypsometry (1m bins)...\n")
      hyps <- terra::freq(terra::round(dem_clipped_m)) |>
        as.data.frame() |>
        dplyr::mutate(
          elevation_m = value,
          area_m2 = count * prod(terra::res(dem_clipped_m)),
          .keep = "unused"
        ) |>
        dplyr::arrange(elevation_m) |>
        dplyr::filter(!is.na(elevation_m))
      
      if (nrow(hyps) > 2) {
        cat(sprintf("  DEM Elevation range: %.1f - %.1f m\n", min(hyps$elevation_m), max(hyps$elevation_m)))
        
        dz <- diff(hyps$elevation_m)
        mean_area <- (hyps$area_m2[-nrow(hyps)] + hyps$area_m2[-1]) / 2
        segment_volume <- dz * mean_area
        
        ve_dem <- hyps |>
          dplyr::mutate(storage_Mm3 = c(0, cumsum(segment_volume)) / 1e6) |>
          dplyr::select(elevation_m, area_m2, storage_Mm3)
        
        results_dem <- interp_storage(ve_dem, TARGET_LEVELS_M) |> dplyr::rename(storage_dem_Mm3 = storage_Mm3)
        area_dem_km2 <- approx(x = ve_dem$elevation_m, y = ve_dem$area_m2 / 1e6, xout = TARGET_LEVELS_M, rule = 1)$y
        results_dem$area_dem_km2 <- round(area_dem_km2, 3)
        
        cat("\n  DEM — Storage at target levels:\n")
        print(results_dem, row.names = FALSE, digits = 4)
      }
    }
  }
}

# ============================================================================
#   # ── SECTION 9 ── COMBINED OUTPUT TABLE
# ============================================================================
print_header("COMBINED RESULTS")
out <- data.frame(level_m = TARGET_LEVELS_M)
if (!is.null(results_3dlakes)) out <- left_join(out, results_3dlakes, by = "level_m")
if (!is.null(results_dahiti)) out <- left_join(out, results_dahiti[, c("level_m", "area_dahiti_km2", "storage_dahiti_Mm3")], by = "level_m")
if (!is.null(results_dem)) out <- left_join(out, results_dem[, c("level_m", "area_dem_km2", "storage_dem_Mm3")], by = "level_m")

if (!is.null(results_3dlakes) && !is.null(results_dahiti) &&
    all(c("storage_3dlakes_Mm3", "storage_dahiti_Mm3") %in% names(out))) {
  out <- out %>%
    mutate(
      diff_Mm3 = round(storage_dahiti_Mm3 - storage_3dlakes_Mm3, 2),
      diff_pct = round((storage_dahiti_Mm3 - storage_3dlakes_Mm3) /
                         pmax(storage_3dlakes_Mm3, 0.1) * 100, 1)
    )
}
cat("\n")
print(out, row.names = FALSE, digits = 5)
write_csv(out, "nechako_hypso_comparison.csv")
cat("\n  → CSV written to: nechako_hypso_comparison.csv\n")

# ============================================================================
#   # ── SECTION 10 ── DIAGNOSTIC PLOT
# ============================================================================
plot_list <- list()
if (!is.null(results_3dlakes) && "storage_3dlakes_Mm3" %in% names(out))
  plot_list[["3D-LAKES"]] <- data.frame(
    level_m     = out$level_m,
    storage_Mm3 = out$storage_3dlakes_Mm3,
    source      = "3D-LAKES (Huang et al. 2025)"
  )
if (!is.null(results_dahiti) && "storage_dahiti_Mm3" %in% names(out))
  plot_list[["DAHITI"]] <- data.frame(
    level_m     = out$level_m,
    storage_Mm3 = out$storage_dahiti_Mm3,
    source      = "DAHITI (DGFI-TUM)"
  )
if (!is.null(results_dem) && "storage_dem_Mm3" %in% names(out))
  plot_list[["DEM"]] <- data.frame(
    level_m     = out$level_m,
    storage_Mm3 = out$storage_dem_Mm3,
    source      = "BC CDED DEM (Local)"
  )

if (length(plot_list) > 0) {
  df_plot <- bind_rows(plot_list) %>% filter(!is.na(storage_Mm3))
  p  <- ggplot(df_plot, aes(x = storage_Mm3, y = level_m,
                            colour = source, shape = source)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 3) +
    scale_colour_manual(
      values = c("3D-LAKES (Huang et al. 2025)" = "#2c7bb6",
                 "DAHITI (DGFI-TUM)"            = "#d7191c",
                 "BC CDED DEM (Local)"          = "#1a9641")
    ) +
    scale_shape_manual(
      values = c("3D-LAKES (Huang et al. 2025)" = 16L,
                 "DAHITI (DGFI-TUM)"            = 17L,
                 "BC CDED DEM (Local)"          = 15L)
    ) +
    labs(
      title    = "Nechako Reservoir — Hypsometric Curve Comparison",
      subtitle = "Storage-elevation (S-E) relationship from three independent sources.",
      x = "Storage (Mm³)", y = "Water Surface Elevation (m ASL)",
      colour = NULL, shape = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold", hjust = 0.5),
      plot.subtitle   = element_text(hjust = 0.5, colour = "grey40", size = 9),
      legend.position = "bottom"
    )
  
  tryCatch({
    ggsave("nechako_hypso_comparison.pdf", p, width = 8, height = 6)
    cat("  → Plot written to: nechako_hypso_comparison.pdf\n")
  }, error = function(e) {
    cat("  ⚠ Failed to save PDF. Error:", conditionMessage(e), "\n")
  })
}

# ── Fine-resolution V-E table for HYPSO_TABLE replacement ────────────────────
ve_for_table <- if (!is.null(ve_3dlakes) && nrow(ve_3dlakes) > 2L) ve_3dlakes else if (!is.null(ve_dem) && nrow(ve_dem) > 2L) ve_dem else NULL
if (!is.null(ve_for_table)) {
  fine_levels  <- seq(842, 855, by = 0.5)
  fine_levels  <- fine_levels[fine_levels >= min(ve_for_table$elevation_m) &
                                fine_levels <= max(ve_for_table$elevation_m)]
  if (length(fine_levels) > 0L) {
    fine_ve  <- interp_storage(ve_for_table, fine_levels)
    source_name <- if (!is.null(ve_3dlakes)) "3D-LAKES" else "BC CDED DEM"
    cat(sprintf("\n  HYPSO_TABLE replacement (%s, 0.5 m resolution)\n", source_name))
    cat("  Paste into your NTSDI script:\n\n")
    cat("HYPSO_TABLE  <- data.frame(\n")
    cat("  level_m     = c(", paste(sprintf("%.1f", fine_ve$level_m),
                                    collapse = ", "), "),\n")
    cat("  storage_Mm3 = c(", paste(sprintf("%.1f", round(fine_ve$storage_Mm3, 1)),
                                    collapse = ", "), ")\n")
    cat(")\n")
  }
}
cat("\n\nScript complete.\n")