# ============================================================================
# 3D-LAKES BATHYMETRY — NECHAKO RIVER BASIN BATCH PROCESSING
# R Script (COMPLETE & OPTIMIZED - WITH H-V CURVE & VOLUME SUMMARY)
# Filename Format: {Hylak_id}_{LEVEL}.csv (e.g., "1_L1.csv", "1000032_L1.csv")
# ============================================================================

# ── INSTALL/LOAD REQUIRED PACKAGES ────────────────────────────────────────────
required_packages <- c("sf", "dplyr", "readr", "data.table", "lubridate",
                       "tidyr", "purrr", "stringr", "future.apply")

for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Disable S2 geometry for better compatibility with shapefiles
sf_use_s2(FALSE)

# ── USER SETTINGS ─────────────────────────────────────────────────────────────
WORKING_DIR <- "D:/Nechako_Drought/Nechako/Lakes/3DLakes"
THREED_LAKES_DIR <- "D:/Nechako_Drought/Nechako/Lakes/3DLakes/"
HYDROLAKES_PATH <- "D:/Nechako_Drought/Nechako/Lakes/HydroLAKES_NorthAmerica.shp"
NECHAKO_BOUNDARY <- "D:/Nechako_Drought/Nechako/Spatial/nechakoBound_dissolve.shp"
MIN_AREA_KM2 <- 1
PRODUCT_LEVEL <- "L1"  # "L1" or "L2"

# H-V Curve Settings
HV_MIN_ELEV_M <- 830    # Minimum elevation for H-V curve (meters)
HV_MAX_ELEV_M <- 860    # Maximum elevation for H-V curve (meters)
HV_STEP_M <- 0.25       # Elevation step (meters)

# Create output directory
if(!dir.exists(WORKING_DIR)) dir.create(WORKING_DIR, recursive = TRUE)
setwd(WORKING_DIR)

# ── HELPER: Safe CRS detection ───────────────────────────────────────────────
get_safe_crs <- function(filepath) {
  tryCatch({
    layers <- sf::st_layers(filepath)
    crs_cols <- c("crs", "wgs84_das", "coord_ref_sys", "srs")
    for(col in crs_cols) {
      if(col %in% names(layers) && length(layers[[col]]) > 0 &&
         !is.na(layers[[col]][1]) && layers[[col]][1] != "") {
        return(st_crs(layers[[col]][1]))
      }
    }
    return(st_crs(4326))
  }, error = function(e) {
    warning("CRS detection failed, assuming EPSG:4326")
    return(st_crs(4326))
  })
}

# ── HELPER: Generate H-V Curve for One Lake ──────────────────────────────────
generate_HV_curve <- function(hylak_id, lake_name, bathy_data, surface_elev, 
                              min_elev, max_elev, step) {
  tryCatch({
    # Check if bathymetry data has required columns
    if(!all(c("elevation", "area") %in% names(bathy_data))) {
      # Fallback: use HydroLAKES volume estimates with parabolic assumption
      cat(sprintf("  ⚠ Using estimated H-V curve for Hylak_id: %d\n", hylak_id))
      
      # Get lake properties from HydroLAKES
      lake_info <- nechako_lakes_df %>% filter(Hylak_id == hylak_id)
      if(nrow(lake_info) == 0) return(NULL)
      
      vol_total_m3 <- lake_info$Vol_total * 1e6  # Convert km³ to m³
      avg_depth <- lake_info$Depth_avg
      lake_area_m2 <- lake_info$Lake_area * 1e6   # Convert km² to m²
      
      # Generate elevation sequence
      elev_seq <- seq(min_elev, max_elev, by = step)
      min_bed_elev <- surface_elev - (avg_depth * 2)  # Estimate max depth
      
      # Parabolic volume estimation: V ∝ (h - h_min)^1.5
      hv_results <- data.frame(
        Hylak_id = hylak_id,
        Lake_name = lake_name,
        Elevation_m = elev_seq,
        Volume_m3 = sapply(elev_seq, function(e) {
          h <- max(0, e - min_bed_elev)
          h_max <- max(0, surface_elev - min_bed_elev)
          if(h_max > 0) {
            vol_total_m3 * (h / h_max)^1.5
          } else {
            0
          }
        }),
        Volume_Mm3 = sapply(elev_seq, function(e) {
          h <- max(0, e - min_bed_elev)
          h_max <- max(0, surface_elev - min_bed_elev)
          if(h_max > 0) {
            (vol_total_m3 * (h / h_max)^1.5) / 1e6
          } else {
            0
          }
        }),
        Method = "HydroLAKES_parabolic_estimate"
      )
      return(hv_results)
    }
    
    # Use actual bathymetry data if available
    cat(sprintf("  ✓ Using 3D-LAKES bathymetry for Hylak_id: %d\n", hylak_id))
    
    elev_seq <- seq(min_elev, max_elev, by = step)
    
    hv_results <- data.frame(
      Hylak_id = hylak_id,
      Lake_name = lake_name,
      Elevation_m = elev_seq,
      Volume_m3 = sapply(elev_seq, function(e) {
        # Calculate volume at each elevation level
        submerged <- bathy_data %>% filter(elevation <= e)
        if(nrow(submerged) == 0) return(0)
        sum(submerged$volume_m3, na.rm = TRUE)
      }),
      Volume_Mm3 = sapply(elev_seq, function(e) {
        submerged <- bathy_data %>% filter(elevation <= e)
        if(nrow(submerged) == 0) return(0)
        sum(submerged$volume_m3, na.rm = TRUE) / 1e6
      }),
      Method = "3D-LAKES_bathymetry"
    )
    
    return(hv_results)
    
  }, error = function(e) {
    cat(sprintf("  ⚠ Error generating H-V curve for Hylak_id %d: %s\n", 
                hylak_id, e$message))
    return(NULL)
  })
}

# ============================================================================
# STEP 1: LOAD NECHAKO BASIN BOUNDARY
# ============================================================================
cat("\n", rep("=", 70), "\n", sep = "")
cat("STEP 1: Loading Nechako Basin Boundary\n")
cat(rep("=", 70), "\n", sep = "")

nechako_basin <- st_read(NECHAKO_BOUNDARY, quiet = TRUE)
nechako_basin <- st_transform(nechako_basin, crs = 4326)
cat("✓ Basin CRS: EPSG:", st_crs(nechako_basin)$epsg, "\n")
cat("✓ Basin polygons:", nrow(nechako_basin), "\n")

# ============================================================================
# STEP 2: LOAD & FILTER HYDROLAKES TO NECHAKO BASIN
# ============================================================================
cat("\n", rep("=", 70), "\n", sep = "")
cat("STEP 2: Loading and Filtering HydroLAKES\n")
cat(rep("=", 70), "\n", sep = "")

if(!file.exists(HYDROLAKES_PATH)) {
  stop("❌ HydroLAKES file not found: ", HYDROLAKES_PATH)
}

hydrolakes_crs <- get_safe_crs(HYDROLAKES_PATH)
cat("✓ HydroLAKES CRS: EPSG:", hydrolakes_crs$epsg, "\n")

cat("Reading HydroLAKES (selected columns only)...\n")
start_time <- Sys.time()
hydrolakes <- tryCatch({
  st_read(HYDROLAKES_PATH, quiet = TRUE,
          select = c("Hylak_id", "Lake_name", "Lake_area", "Vol_total",
                     "Depth_avg", "Elevation"))
}, error = function(e) {
  cat("⚠ Fallback: Reading full dataset...\n")
  st_read(HYDROLAKES_PATH, quiet = FALSE) %>%
    select(Hylak_id, Lake_name, Lake_area, Vol_total, Depth_avg, Elevation)
})
load_time <- round(difftime(Sys.time(), start_time, units = "secs"), 2)
cat("✓ Load time:", load_time, "seconds | Loaded", nrow(hydrolakes), "lakes\n")

if(!is.na(st_crs(hydrolakes)$epsg) && st_crs(hydrolakes)$epsg != 4326) {
  cat("Transforming to WGS84...\n")
  hydrolakes <- st_transform(hydrolakes, crs = 4326)
}

# Spatial filter
cat("Applying spatial filter (st_intersects)...\n")
within_basin <- lengths(st_intersects(hydrolakes, nechako_basin)) > 0
hydrolakes_nechako <- hydrolakes[within_basin, ]
cat("✓ Lakes within basin boundary:", nrow(hydrolakes_nechako), "\n")

nechako_lakes <- hydrolakes_nechako %>%
  filter(Lake_area >= MIN_AREA_KM2) %>%
  arrange(desc(Lake_area))
cat("✓ Final: Lakes in Nechako Basin (≥", MIN_AREA_KM2, "km²):", nrow(nechako_lakes), "\n")

# Preview top lakes
cat("\nTop 10 lakes by area:\n")
print(nechako_lakes %>% select(Hylak_id, Lake_name, Lake_area) %>% head(10))

# ============================================================================
# STEP 3: EXPORT HYDROLAKES INVENTORY (CSV + SHAPEFILE)
# ============================================================================
cat("\n", rep("=", 70), "\n", sep = "")
cat("STEP 3: Exporting HydroLAKES Inventory (CSV + Shapefile)\n")
cat(rep("=", 70), "\n", sep = "")

# Export CSV (drop geometry)
nechako_lakes_df <- nechako_lakes %>% st_drop_geometry() %>% as.data.frame()
output_csv <- file.path(WORKING_DIR, "Nechako_HydroLAKES_all_lakes.csv")
write_csv(nechako_lakes_df, output_csv)
cat("✓ Exported CSV:", output_csv, "\n")

# Export Shapefile (with geometry)
output_shp <- file.path(WORKING_DIR, "Nechako_HydroLAKES_all_lakes.shp")
cat("Exporting Shapefile (with geometry)...\n")
if(is.na(st_crs(nechako_lakes)$epsg)) {
  st_crs(nechako_lakes) <- st_crs(4326)
}
st_write(nechako_lakes,
         dsn = output_shp,
         driver = "ESRI Shapefile",
         delete_dsn = TRUE,
         quiet = TRUE,
         check_exists = FALSE)
cat("✓ Exported Shapefile:", output_shp, "\n")

shp_files <- list.files(WORKING_DIR, pattern = "Nechako_HydroLAKES_all_lakes\\.",
                        full.names = TRUE)
cat("✓ Shapefile components created:", length(shp_files), "files\n")

# ============================================================================
# STEP 4: MATCH WITH 3D-LAKES BATHYMETRY DATA
# ============================================================================
cat("\n", rep("=", 70), "\n", sep = "")
cat("STEP 4: Matching with 3D-Lakes", PRODUCT_LEVEL, "Dataset\n")
cat(rep("=", 70), "\n", sep = "")

hylak_ids <- nechako_lakes_df$Hylak_id
cat("✓ Hylak_ids to search:", length(hylak_ids), "\n")

product_folder <- file.path(THREED_LAKES_DIR,
                            if(PRODUCT_LEVEL == "L2") "L2 A-E product" else "L1 A-E product")
if(!dir.exists(product_folder)) {
  stop("❌ Product folder not found: ", product_folder)
}

cat("Searching for 3D-Lakes files (direct lookup)...\n")
start_time <- Sys.time()
matched_files <- character(0)
matched_ids <- numeric(0)

for(hylak_id in hylak_ids) {
  test_file <- file.path(product_folder, paste0(hylak_id, "_", PRODUCT_LEVEL, ".csv"))
  if(file.exists(test_file)) {
    matched_files <- c(matched_files, test_file)
    matched_ids <- c(matched_ids, hylak_id)
  }
}

end_time <- Sys.time()
cat("✓ Lookup time:", round(difftime(end_time, start_time, units = "secs"), 2), "seconds\n")
cat("✓ Lakes with", PRODUCT_LEVEL, "data:", length(matched_files), "\n")

if(length(matched_files) > 0) {
  matched_rows <- data.frame(Hylak_id = matched_ids, file_path = matched_files)
  cat("\nMatched lakes:\n")
  for(i in seq_len(nrow(matched_rows))) {
    lake_name <- nechako_lakes_df %>%
      filter(Hylak_id == matched_ids[i]) %>%
      pull(Lake_name) %>%
      first()
    cat(sprintf(" [%d] Hylak_id: %d | Name: %s\n",
                i, matched_ids[i], ifelse(is.na(lake_name), "Unnamed", lake_name)))
  }
} else {
  cat("⚠️ WARNING: No 3D-Lakes", PRODUCT_LEVEL, "files found for Nechako lakes!\n")
  cat(" Possible reasons:\n")
  cat(" 1. Lakes are not in the 3D-Lakes dataset\n")
  cat(" 2. Try checking the other product level (L1 vs L2)\n")
  cat(" 3. Verify Hylak_id format matches between datasets\n")
}

# ============================================================================
# STEP 5: PROCESS & EXPORT BATHYMETRY DATA (PARALLEL)
# ============================================================================
cat("\n", rep("=", 70), "\n", sep = "")
cat("STEP 5: Processing Bathymetry Data (Parallel)\n")
cat(rep("=", 70), "\n", sep = "")

process_3dlake_file <- function(file_path, hylak_id) {
  tryCatch({
    lake_data <- fread(file_path, showProgress = FALSE)
    lake_data %>% mutate(Hylak_id = hylak_id, Source_File = basename(file_path))
  }, error = function(e) {
    cat("⚠ Error reading", basename(file_path), ":", e$message, "\n")
    return(NULL)
  })
}

bathymetry_list <- list()
if(exists("matched_files") && length(matched_files) > 0) {
  n_cores <- max(1, parallel::detectCores() - 1)
  cat(sprintf("Processing %d lakes using %d cores...\n", length(matched_files), n_cores))
  
  bathymetry_list <- future.apply::future_lapply(
    seq_along(matched_files),
    function(i) {
      if(i %% 10 == 0) cat(sprintf(" Progress: %d/%d lakes\n", i, length(matched_files)))
      
      hylak_id <- matched_ids[i]
      bathy_data <- process_3dlake_file(matched_files[i], hylak_id)
      
      if(!is.null(bathy_data)) {
        safe_name <- stringr::str_replace_all(
          nechako_lakes_df %>% filter(Hylak_id == hylak_id) %>% pull(Lake_name) %>% first(),
          "[^A-Za-z0-9]", "_")
        if(is.na(safe_name) || safe_name == "") safe_name <- "Unnamed"
        
        output_file <- file.path(WORKING_DIR,
                                 sprintf("Hylak_%d_%s_3DLakes_%s.csv",
                                         hylak_id, safe_name, PRODUCT_LEVEL))
        fwrite(bathy_data, output_file)
      }
      return(bathy_data)
    },
    future.seed = TRUE)
  
  bathymetry_list <- Filter(Negate(is.null), bathymetry_list)
  
  if(length(bathymetry_list) > 0) {
    cat("\nCombining bathymetry data...\n")
    all_bathymetry <- data.table::rbindlist(bathymetry_list, use.names = TRUE, fill = TRUE)
    combined_output <- file.path(WORKING_DIR,
                                 sprintf("Nechako_3DLakes_%s_combined.csv", PRODUCT_LEVEL))
    fwrite(all_bathymetry, combined_output)
    cat("✓ Combined bathymetry:", combined_output, "\n")
    cat("✓ Total rows:", nrow(all_bathymetry), "\n")
  }
} else {
  cat("⚠ No bathymetry data to process.\n")
}

# ============================================================================
# STEP 6: GENERATE H-V CURVES FOR ALL LAKES
# ============================================================================
cat("\n", rep("=", 70), "\n", sep = "")
cat("STEP 6: Generating Height-Volume (H-V) Curves\n")
cat(rep("=", 70), "\n", sep = "")

hv_curves_list <- list()

for(i in seq_len(nrow(nechako_lakes_df))) {
  hylak_id <- nechako_lakes_df$Hylak_id[i]
  lake_name <- nechako_lakes_df$Lake_name[i]
  surface_elev <- nechako_lakes_df$Elevation[i]
  
  if(is.na(lake_name) || lake_name == "") {
    lake_name <- paste0("Unnamed_", hylak_id)
  }
  
  cat(sprintf("[%d/%d] Processing H-V curve for: %s (Hylak_id: %d)\n",
              i, nrow(nechako_lakes_df), lake_name, hylak_id))
  
  # Check if bathymetry data exists for this lake
  bathy_data <- NULL
  if(exists("matched_ids") && hylak_id %in% matched_ids) {
    idx <- which(matched_ids == hylak_id)
    if(length(idx) > 0 && exists("bathymetry_list") && length(bathymetry_list) >= idx) {
      bathy_data <- bathymetry_list[[idx]]
    }
  }
  
  # Generate H-V curve
  hv_curve <- generate_HV_curve(
    hylak_id = hylak_id,
    lake_name = lake_name,
    bathy_data = bathy_data,
    surface_elev = surface_elev,
    min_elev = HV_MIN_ELEV_M,
    max_elev = HV_MAX_ELEV_M,
    step = HV_STEP_M
  )
  
  if(!is.null(hv_curve) && nrow(hv_curve) > 0) {
    hv_curves_list[[length(hv_curves_list) + 1]] <- hv_curve
    
    # Export individual H-V curve CSV
    safe_name <- stringr::str_replace_all(lake_name, "[^A-Za-z0-9]", "_")
    hv_output_file <- file.path(WORKING_DIR,
                                sprintf("HV_Curve_%d_%s.csv", hylak_id, safe_name))
    write_csv(hv_curve, hv_output_file)
  }
}

# Combine all H-V curves
if(length(hv_curves_list) > 0) {
  all_hv_curves <- data.table::rbindlist(hv_curves_list, use.names = TRUE, fill = TRUE)
  hv_combined_output <- file.path(WORKING_DIR, "Nechako_All_HV_Curves_combined.csv")
  fwrite(all_hv_curves, hv_combined_output)
  cat("✓ Combined H-V curves exported:", hv_combined_output, "\n")
  cat("✓ Total H-V curve records:", nrow(all_hv_curves), "\n")
} else {
  cat("⚠ No H-V curves generated.\n")
}

# ============================================================================
# STEP 7: CREATE SUMMARY REPORT WITH VOLUME STATISTICS
# ============================================================================
cat("\n", rep("=", 70), "\n", sep = "")
cat("STEP 7: Creating Summary Report with Volume Statistics\n")
cat(rep("=", 70), "\n", sep = "")

# Calculate volume statistics
total_lakes_in_basin <- nrow(nechako_lakes_df)
lakes_with_bathymetry <- if(exists("matched_ids")) length(matched_ids) else 0

# Sum of volumes from HydroLAKES Vol_total column (in km³)
total_volume_km3 <- sum(nechako_lakes_df$Vol_total, na.rm = TRUE)
total_volume_Mm3 <- total_volume_km3 * 1e6  # Convert to million cubic meters

# Volume of lakes with bathymetry data
if(exists("matched_ids") && length(matched_ids) > 0) {
  bathy_lakes_df <- nechako_lakes_df %>% filter(Hylak_id %in% matched_ids)
  volume_with_bathymetry_km3 <- sum(bathy_lakes_df$Vol_total, na.rm = TRUE)
  volume_with_bathymetry_Mm3 <- volume_with_bathymetry_km3 * 1e6
} else {
  volume_with_bathymetry_km3 <- 0
  volume_with_bathymetry_Mm3 <- 0
}

# Volume of lakes WITHOUT bathymetry data
volume_without_bathymetry_km3 <- total_volume_km3 - volume_with_bathymetry_km3
volume_without_bathymetry_Mm3 <- total_volume_Mm3 - volume_with_bathymetry_Mm3

# Coverage percentage
bathymetry_coverage_pct <- if(total_lakes_in_basin > 0) {
  (lakes_with_bathymetry / total_lakes_in_basin) * 100
} else {
  0
}

volume_coverage_pct <- if(total_volume_km3 > 0) {
  (volume_with_bathymetry_km3 / total_volume_km3) * 100
} else {
  0
}

summary_report <- list(
  Basin = "Nechako River Basin",
  Analysis_Date = as.character(Sys.Date()),
  Total_Lakes_In_Basin = total_lakes_in_basin,
  Lakes_With_3D_Data = lakes_with_bathymetry,
  Lakes_Without_3D_Data = total_lakes_in_basin - lakes_with_bathymetry,
  Bathymetry_Coverage_Pct = sprintf("%.2f%%", bathymetry_coverage_pct),
  Min_Area_KM2 = MIN_AREA_KM2,
  Total_Volume_HydroLAKES_km3 = sprintf("%.3f", total_volume_km3),
  Total_Volume_HydroLAKES_Mm3 = sprintf("%.0f", total_volume_Mm3),
  Volume_With_Bathymetry_km3 = sprintf("%.3f", volume_with_bathymetry_km3),
  Volume_With_Bathymetry_Mm3 = sprintf("%.0f", volume_with_bathymetry_Mm3),
  Volume_Without_Bathymetry_km3 = sprintf("%.3f", volume_without_bathymetry_km3),
  Volume_Without_Bathymetry_Mm3 = sprintf("%.0f", volume_without_bathymetry_Mm3),
  Volume_Coverage_Pct = sprintf("%.2f%%", volume_coverage_pct),
  Product_Level_Used = PRODUCT_LEVEL,
  HV_Curve_Elev_Range = sprintf("%.1f-%.1f m (step: %.2f m)", 
                                HV_MIN_ELEV_M, HV_MAX_ELEV_M, HV_STEP_M),
  Output_Directory = WORKING_DIR,
  Optimization_Method = "Direct Lookup + Parallel Processing + S2 Disabled + Shapefile Export + H-V Curves"
)

cat("\n══ SUMMARY ══\n")
for(name in names(summary_report)) {
  cat(sprintf("%-35s: %s\n", name, summary_report[[name]]))
}

# Export summary report
summary_file <- file.path(WORKING_DIR, "Analysis_Summary.txt")
sink(summary_file)
cat("Nechako Basin 3D-Lakes Analysis Summary\n")
cat(paste(rep("=", 50), collapse = ""), "\n\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
for(name in names(summary_report)) {
  cat(sprintf("%-35s: %s\n", name, summary_report[[name]]))
}
cat("\n\n══ VOLUME STATISTICS ══\n")
cat(sprintf("Total lakes in basin:              %d\n", total_lakes_in_basin))
cat(sprintf("Lakes with bathymetry data:        %d (%.2f%%)\n", 
            lakes_with_bathymetry, bathymetry_coverage_pct))
cat(sprintf("Lakes without bathymetry data:     %d\n", 
            total_lakes_in_basin - lakes_with_bathymetry))
cat("\n")
cat(sprintf("Total volume (all lakes):          %.3f km³ (%.0f Mm³)\n", 
            total_volume_km3, total_volume_Mm3))
cat(sprintf("Volume with bathymetry:            %.3f km³ (%.0f Mm³)\n", 
            volume_with_bathymetry_km3, volume_with_bathymetry_Mm3))
cat(sprintf("Volume without bathymetry:         %.3f km³ (%.0f Mm³)\n", 
            volume_without_bathymetry_km3, volume_without_bathymetry_Mm3))
cat(sprintf("Volume coverage:                   %.2f%%\n", volume_coverage_pct))
cat("\n\n══ H-V CURVE SETTINGS ══\n")
cat(sprintf("Elevation range:                   %.1f - %.1f m\n", 
            HV_MIN_ELEV_M, HV_MAX_ELEV_M))
cat(sprintf("Elevation step:                    %.2f m\n", HV_STEP_M))
cat(sprintf("Total elevation levels per lake:   %d\n", 
            length(seq(HV_MIN_ELEV_M, HV_MAX_ELEV_M, by = HV_STEP_M))))
sink()
cat("✓ Summary exported:", summary_file, "\n")

# Export volume statistics as CSV
volume_stats_df <- data.frame(
  Metric = c("Total_Lakes", "Lakes_With_Bathymetry", "Lakes_Without_Bathymetry",
             "Total_Volume_km3", "Total_Volume_Mm3",
             "Volume_With_Bathymetry_km3", "Volume_With_Bathymetry_Mm3",
             "Volume_Without_Bathymetry_km3", "Volume_Without_Bathymetry_Mm3",
             "Bathymetry_Coverage_Pct", "Volume_Coverage_Pct"),
  Value = c(total_lakes_in_basin, lakes_with_bathymetry, 
            total_lakes_in_basin - lakes_with_bathymetry,
            total_volume_km3, total_volume_Mm3,
            volume_with_bathymetry_km3, volume_with_bathymetry_Mm3,
            volume_without_bathymetry_km3, volume_without_bathymetry_Mm3,
            bathymetry_coverage_pct, volume_coverage_pct)
)
volume_stats_file <- file.path(WORKING_DIR, "Volume_Statistics.csv")
write_csv(volume_stats_df, volume_stats_file)
cat("✓ Volume statistics exported:", volume_stats_file, "\n")

# ============================================================================
# CLEANUP & FINAL MESSAGE
# ============================================================================
rm(list = ls()[!ls() %in% c("WORKING_DIR", "PRODUCT_LEVEL")])
gc()

cat("\n", rep("=", 70), "\n", sep = "")
cat("✅ PROCESSING COMPLETE!\n")
cat(rep("=", 70), "\n", sep = "")
cat("Output files are in:", WORKING_DIR, "\n\n")
cat("Key outputs:\n")
cat(" • Nechako_HydroLAKES_all_lakes.csv → Lake inventory (CSV)\n")
cat(" • Nechako_HydroLAKES_all_lakes.shp → Lake boundaries (Shapefile)\n")
cat(" • Nechako_3DLakes_", PRODUCT_LEVEL, "_combined.csv → All bathymetry data\n", sep = "")
cat(" • Hylak_*_3DLakes_", PRODUCT_LEVEL, ".csv → Individual lake files\n", sep = "")
cat(" • HV_Curve_*_*.csv → Individual H-V curves\n")
cat(" • Nechako_All_HV_Curves_combined.csv → Combined H-V curves\n")
cat(" • Volume_Statistics.csv → Volume summary statistics\n")
cat(" • Analysis_Summary.txt → Processing report\n")
cat("\nNext steps:\n")
cat(" 1. Open Nechako_HydroLAKES_all_lakes.shp in QGIS/ArcGIS for mapping\n")
cat(" 2. Review Volume_Statistics.csv for bathymetry coverage assessment\n")
cat(" 3. Use HV_Curve files for SWSI drought index calculations\n")
cat(" 4. Individual lake files can be used for detailed bathymetric analysis\n")