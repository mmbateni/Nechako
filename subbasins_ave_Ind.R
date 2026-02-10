###############################################################################
# Averages of Nechako's sub-basins  (CRS-aware, diagnostics)
# - FWA shapefile CRS: EPSG:3005 (NAD83 / BC Albers)
# - Index stacks stored in lat/lon (EPSG:4326 or similar)
# - Script reprojects basin polygon to raster CRS before extraction
# - Adds diagnostics: CRS, extent, resolution, nlayers, time range, NA counts
# - Computes averages for EACH sub-basin separately
#
# Edit the path variables below to match your environment.
###############################################################################

library(sf)
library(terra)
library(dplyr)
library(ggplot2)
library(lubridate)
library(openxlsx)
library(rlang)

# ---------------------------
# editable paths
# ---------------------------
base_dir <- "D:/Nechako_Drought"                 # project root
fwa_shp_path <- file.path(base_dir, "Spatial", "FWWTRSHDGR_polygon.shp") 
spi_dir  <- file.path(base_dir, "spi_results_seasonal")
spei_dir <- file.path(base_dir, "spei_results_seasonal")
output_dir <- file.path(base_dir, "basin_averaged_outputs")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# timescales to process (months)
timescales <- c(1, 3, 6, 9, 12)

# expected shapefile CRS (informational)
expected_shp_epsg <- 3005  # NAD83 / BC Albers


# ---------------------------
# Utility: print diagnostics for a SpatRaster
# ---------------------------
print_raster_diag <- function(r, name = "raster") {
  e <- ext(r)
  extent_str <- sprintf("xmin=%.4f, xmax=%.4f, ymin=%.4f, ymax=%.4f", 
                        e$xmin, e$xmax, e$ymin, e$ymax)
  
  cat("---- Raster diagnostics:", name, "----\n")
  cat("CRS:", crs(r), "\n")
  cat("Extent:", extent_str, "\n")
  cat("Resolution:", paste(round(res(r), 6), collapse = " x "), "\n")
  cat("Number of layers:", nlyr(r), "\n")
  
  # time metadata
  tvec <- time(r)
  if (!is.null(tvec) && !all(is.na(tvec))) {
    # Convert numeric time (days since origin) to Date for readable output
    if (is.numeric(tvec)) {
      tvec_dates <- as.Date(tvec, origin = "1970-01-01")
      cat("Time range:", min(tvec_dates, na.rm = TRUE), "to", max(tvec_dates, na.rm = TRUE), "\n")
    } else {
      cat("Time range:", min(tvec, na.rm = TRUE), "to", max(tvec, na.rm = TRUE), "\n")
    }
  } else {
    cat("Time metadata: none detected\n")
  }
  
  # CRITICAL FIXES:
  na_frac_val <- global(r[[1]], fun = function(x, ...) mean(is.na(x)), na.rm = FALSE)
  na_frac <- as.numeric(na_frac_val[1, 1])  # Extract scalar value
  
  cat("NA fraction (layer1):", round(na_frac, 4), "\n")
  cat("----------------------------------------\n")
}
# ---------------------------
# 1) Read FWA watershed groups and identify Nechako sub-basins
# ---------------------------
cat("Reading FWA watershed groups shapefile...\n")
fwa <- st_read(fwa_shp_path, quiet = TRUE)

cat("Shapefile CRS (reported):", st_crs(fwa)$input, "\n")
if (!is.na(st_crs(fwa)$epsg)) {
  cat("Shapefile EPSG:", st_crs(fwa)$epsg, "\n")
}

# Confirm expected CRS
if (!is.na(st_crs(fwa)$epsg) && st_crs(fwa)$epsg != expected_shp_epsg) {
  warning(sprintf("Shapefile EPSG is %s but expected %d. Proceeding with actual CRS.", st_crs(fwa)$epsg, expected_shp_epsg))
}

# Detect the boolean flag field for Nechako
flag_field_candidates <- c("Is_NRB", "Is NRB", "Is_NRB_", "IS_NRB", "IsNRB", "Is_Nrb", "Is_NRB", "IsNRB")
flag_field <- intersect(names(fwa), flag_field_candidates)
if (length(flag_field) == 0) {
  # fallback: search for fields containing 'NRB' or 'nechako'
  flag_field <- names(fwa)[grepl("NRB|nrb|Nechako|nechako", names(fwa), ignore.case = TRUE)]
  if (length(flag_field) == 0) {
    cat("Available fields:\n"); print(names(fwa))
    stop("Cannot find an 'Is NRB' field. Inspect names(fwa) and set the correct field name in the script.")
  } else {
    flag_field <- flag_field[1]
    cat("Auto-detected flag field:", flag_field, "\n")
  }
} else {
  flag_field <- flag_field[1]
  cat("Using flag field:", flag_field, "\n")
}

# Select polygons where flag is TRUE (handle logical, numeric, text)
fwa_selected <- fwa %>%
  filter((!!sym(flag_field)) %in% c(TRUE, 1, "1", "TRUE", "True", "true"))

if (nrow(fwa_selected) == 0) {
  stop("No polygons selected for Nechako (Is_NRB TRUE). Check attribute values and field name.")
}
cat(sprintf("Selected %d sub-basins flagged as Nechako.\n", nrow(fwa_selected)))

# MODIFIED: Keep individual sub-basins instead of dissolving
nechako_subbasins <- st_make_valid(fwa_selected)
cat("Keeping individual sub-basins for separate analysis.\n")

# Save the sub-basins for inspection
st_write(nechako_subbasins, file.path(output_dir, "nechako_subbasins.gpkg"), delete_dsn = TRUE, quiet = TRUE)
cat("Saved: nechako_subbasins.gpkg\n")

# ---------------------------
# Helper: load monthly raster stack for an index and timescale
# - robust filename matching and time parsing
# ---------------------------
load_index_stack <- function(index_type = c("spi", "spei"), scale) {
  index_type <- match.arg(index_type)
  folder <- if (index_type == "spi") spi_dir else spei_dir
  # candidate patterns
  patterns <- c(
    sprintf("(?i)%s.*scale[_-]?%02d.*\\.tif$", index_type, scale),
    sprintf("(?i)%s.*scale[_-]?%d.*\\.tif$", index_type, scale),
    sprintf("(?i)%s.*scale[_-]?%02d.*\\.nc$", index_type, scale),
    sprintf("(?i)%s.*scale[_-]?%d.*\\.nc$", index_type, scale)
  )
  files <- unlist(lapply(patterns, function(p) list.files(folder, pattern = p, full.names = TRUE)))
  files <- unique(files[file.exists(files)])
  if (length(files) == 0) {
    # try looser pattern: index + scale anywhere
    files <- list.files(folder, pattern = sprintf("(?i)%s.*%02d.*\\.(tif|nc)$", index_type, scale), full.names = TRUE)
  }
  if (length(files) == 0) {
    stop(sprintf("No raster stack found for %s scale %d in %s. Place a multi-layer GeoTIFF or NetCDF there.", index_type, scale, folder))
  }
  r <- rast(files[1])
  # diagnostics
  print_raster_diag(r, basename(files[1]))
  # ensure time metadata exists or attempt to parse layer names
  tvec <- time(r)
  if (is.null(tvec) || all(is.na(tvec))) {
    ln <- names(r)
    parsed <- sapply(ln, function(x) {
      m <- regmatches(x, regexpr("\\d{4}[_-]?\\d{2}", x))
      if (length(m) && nchar(m) > 0) {
        m2 <- gsub("[-_]", "", m)
        as.Date(paste0(substr(m2,1,4), "-", substr(m2,5,6), "-01"))
      } else NA
    })
    if (all(is.na(parsed))) {
      # fallback: ask user to confirm start date or assume 1950-01-01
      start <- as.Date("1950-01-01")
      parsed <- seq(start, by = "1 month", length.out = nlyr(r))
      warning(sprintf("No time metadata found for %s scale %d; assigning synthetic monthly dates starting %s.", index_type, scale, start))
    }
    time(r) <- parsed
  }
  return(r)
}

# ---------------------------
# MODIFIED: Compute sub-basin mean time series for each sub-basin
# - reprojects basin polygons to raster CRS before extraction
# ---------------------------
calculate_subbasin_averages <- function(index_type = c("spi", "spei"), scale) {
  index_type <- match.arg(index_type)
  cat(sprintf("Loading %s scale %d ...\n", index_type, scale))
  r <- load_index_stack(index_type, scale)
  # raster CRS
  raster_crs <- crs(r)
  cat("Raster CRS (terra):", raster_crs, "\n")
  # reproject sub-basins to raster CRS
  subbasins_proj <- st_transform(nechako_subbasins, raster_crs)
  # convert to terra vect
  subbasins_vect <- vect(subbasins_proj)
  
  # Identify name field for sub-basins
  name_field_candidates <- c("WTRSHDGRPN", "WTRSHD_NAM", "NAME", "Basin_Name", "GNIS_NAME")
  name_field <- intersect(names(subbasins_proj), name_field_candidates)
  if (length(name_field) == 0) {
    # Use first text column as name
    text_cols <- names(subbasins_proj)[sapply(subbasins_proj, is.character)]
    if (length(text_cols) > 0) {
      name_field <- text_cols[1]
    } else {
      name_field <- "ID"
      subbasins_proj$ID <- seq_len(nrow(subbasins_proj))
    }
  } else {
    name_field <- name_field[1]
  }
  cat(sprintf("Using '%s' field for sub-basin names\n", name_field))
  
  # Extract zonal mean for each sub-basin
  cat("Extracting zonal mean for each sub-basin (terra::extract) ...\n")
  ex <- terra::extract(r, subbasins_vect, fun = mean, na.rm = TRUE)
  
  if (is.null(ex) || nrow(ex) == 0) stop("terra::extract returned no values.")
  
  # Get time info
  dates <- time(r)
  
  # Combine results into long format dataframe
  results_list <- list()
  for (i in 1:nrow(ex)) {
    subbasin_name <- as.character(subbasins_proj[[name_field]][i])
    vals <- as.numeric(ex[i, -1])  # drop ID column
    
    df_sub <- data.frame(
      date = as.Date(dates),
      subbasin = subbasin_name,
      mean_value = vals,
      index = index_type,
      scale = scale
    )
    
    na_count <- sum(is.na(df_sub$mean_value))
    cat(sprintf("  Sub-basin '%s': %d months, %d NA months\n", subbasin_name, nrow(df_sub), na_count))
    
    results_list[[i]] <- df_sub
  }
  
  # Combine all sub-basins
  results_df <- do.call(rbind, results_list)
  
  return(results_df)
}

# ---------------------------
# Drought detection helper (modified for sub-basins)
# ---------------------------
detect_drought_events <- function(df, onset = -1.0, termination = 0.0, min_dur = 2) {
  if (!"mean_value" %in% names(df)) stop("df must contain mean_value and date")
  n <- nrow(df)
  events <- list()
  i <- 1
  while (i <= n) {
    if (!is.na(df$mean_value[i]) && df$mean_value[i] < onset) {
      start_i <- i
      j <- i + 1
      while (j <= n && (!is.na(df$mean_value[j]) && df$mean_value[j] < termination)) j <- j + 1
      end_i <- min(j - 1, n)
      dur <- end_i - start_i + 1
      if (dur >= min_dur) {
        min_val <- min(df$mean_value[start_i:end_i], na.rm = TRUE)
        events[[length(events) + 1]] <- data.frame(
          start_date = df$date[start_i],
          end_date = df$date[end_i],
          duration_months = dur,
          min_value = min_val
        )
      }
      i <- end_i + 1
    } else {
      i <- i + 1
    }
  }
  if (length(events) == 0) return(data.frame())
  return(do.call(rbind, events))
}

# ---------------------------
# Plot helper (modified for sub-basins)
# ---------------------------
create_timeseries_plot <- function(df, title = NULL) {
  p <- ggplot(df, aes(x = date, y = mean_value)) +
    geom_line(color = "steelblue") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    facet_wrap(~ subbasin, ncol = 3, scales = "free_y") +
    labs(x = "Date", y = "Mean value", title = title) +
    theme_minimal() +
    theme(strip.text = element_text(size = 8))
  return(p)
}

# ---------------------------
# MAIN: loop through indices and timescales, save CSVs and PNGs
# ---------------------------
for (index_type in c("spi", "spei")) {
  cat(sprintf("\nProcessing index: %s\n", toupper(index_type)))
  
  for (scale in timescales) {
    df <- calculate_subbasin_averages(index_type, scale)
    
    # Save combined CSV with all sub-basins
    csv_out <- file.path(output_dir, sprintf("%s_subbasin_means_scale%02d.csv", index_type, scale))
    write.csv(df, csv_out, row.names = FALSE)
    cat(sprintf("  + Saved CSV: %s\n", csv_out))
    
    # Save a PNG plot with all sub-basins in facets
    p <- create_timeseries_plot(df, sprintf("%s - Sub-basin means (scale %d months)", toupper(index_type), scale))
    png_out <- file.path(output_dir, sprintf("%s_subbasin_means_scale%02d.png", index_type, scale))
    ggsave(png_out, p, width = 16, height = 12)
    cat(sprintf("  + Saved plot: %s\n", png_out))
    
    # Also save individual plots for each sub-basin for clarity
    individual_plot_dir <- file.path(output_dir, sprintf("%s_scale%02d_individual_plots", index_type, scale))
    dir.create(individual_plot_dir, showWarnings = FALSE, recursive = TRUE)
    
    for (subbasin_name in unique(df$subbasin)) {
      df_sub <- df %>% filter(subbasin == subbasin_name)
      p_individual <- ggplot(df_sub, aes(x = date, y = mean_value)) +
        geom_line(color = "steelblue", linewidth = 0.7) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
        labs(x = "Date", y = "Mean value", 
             title = sprintf("%s - %s (scale %d months)", toupper(index_type), subbasin_name, scale)) +
        theme_minimal()
      
      # Clean filename
      safe_name <- gsub("[:/\\\\?*\\[\\] ]", "_", subbasin_name)
      png_individual <- file.path(individual_plot_dir, sprintf("%s.png", safe_name))
      ggsave(png_individual, p_individual, width = 10, height = 4)
    }
    cat(sprintf("  + Saved %d individual plots in: %s\n", length(unique(df$subbasin)), basename(individual_plot_dir)))
    
    # Detect droughts for recent period 2020-2025 for each sub-basin
    df_recent <- df %>% filter(date >= as.Date("2020-01-01") & date <= as.Date("2025-12-31"))
    
    all_events <- list()
    for (subbasin_name in unique(df_recent$subbasin)) {
      df_sub <- df_recent %>% filter(subbasin == subbasin_name)
      events <- detect_drought_events(df_sub, onset = -1.0, termination = 0.0, min_dur = 1)
      if (nrow(events) > 0) {
        events$subbasin <- subbasin_name
        all_events[[length(all_events) + 1]] <- events
        cat(sprintf("  + Sub-basin '%s': %d drought events (2020-2025)\n", subbasin_name, nrow(events)))
      }
    }
    
    if (length(all_events) > 0) {
      events_combined <- do.call(rbind, all_events)
      ev_out <- file.path(output_dir, sprintf("%s_drought_events_scale%02d_2020-2025.csv", index_type, scale))
      write.csv(events_combined, ev_out, row.names = FALSE)
      cat(sprintf("  + Saved drought events: %s\n", ev_out))
    } else {
      cat("  + No drought events detected (2020-2025) for any sub-basin.\n")
    }
    
    # Create Excel workbook with separate sheets for each sub-basin
    wb <- createWorkbook()
    for (subbasin_name in unique(df$subbasin)) {
      # Clean sheet name (Excel has 31 char limit and doesn't allow special chars)
      sheet_name <- gsub("[:/\\\\?*\\[\\]]", "_", subbasin_name)
      sheet_name <- substr(sheet_name, 1, 31)
      
      df_sub <- df %>% filter(subbasin == subbasin_name)
      addWorksheet(wb, sheet_name)
      writeData(wb, sheet_name, df_sub)
    }
    xlsx_out <- file.path(output_dir, sprintf("%s_subbasin_timeseries_scale%02d.xlsx", index_type, scale))
    saveWorkbook(wb, xlsx_out, overwrite = TRUE)
    cat(sprintf("  + Saved Excel workbook: %s\n", xlsx_out))
  }
}

cat("\nAll processing complete. Outputs in:\n", normalizePath(output_dir), "\n")