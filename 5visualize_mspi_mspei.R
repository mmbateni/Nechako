##############################################
# BASIN-AVERAGED MSPI & MSPEI VISUALIZATION
# - MSPI/MSPEI are SINGLE multi-scalar indices (no timescales)
# - Generates SEPARATE PDFs for each index + comparison plot
# - Drought definition: onset < -1.0, termination >= 0.0
# - Shaded drought periods with severity labels
# - Titles EMBEDDED in plots (PDF-safe)
# - Excel summary with detailed drought event tracking
# - NEW: Spatial maps for 20 months (Jan 2022-Dec 2025)
##############################################

library(openxlsx)
library(ggplot2)
library(scales)     # Required for oob = scales::squish
library(gridExtra)
library(grid)
library(lubridate)
library(terra)
library(sf)

#---- Helper: Safe PDF creation with automatic device cleanup ----
safe_pdf <- function(filename, width = 12, height = 8) {
  while (dev.cur() > 1) try(dev.off(), silent = TRUE)
  
  tryCatch({
    pdf(filename, width = width, height = height, family = "Helvetica")
    TRUE
  }, error = function(e) {
    cat(sprintf("ERROR: Cannot create PDF '%s': %s\n", filename, e$message))
    FALSE
  })
}

#---- Paths ----
# Set this to your project root. All other paths are relative to it.
# Adjust to your environment; do NOT use a hardcoded absolute path.
project_dir <- Sys.getenv("NECHAKO_PROJECT_DIR",
                          unset = getwd())  # fallback: current working directory
setwd(project_dir)
cat(sprintf("Working directory: %s\n", getwd()))

# CHANGED: Output directories match MSPI_calculation.R
mspi_output_dir <- "mspi_results"
mspei_output_dir <- "mspei_results"

# Create directories if they don't exist
dir.create(mspi_output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(mspei_output_dir, showWarnings = FALSE, recursive = TRUE)

mspi_file <- "mspi_results/mspi_basin_timeseries_1950-2025.xlsx"
mspei_file <- "mspei_results/mspei_basin_timeseries_1950-2025.xlsx"

# NEW: NetCDF files for spatial maps
mspi_nc_file <- "mspi_results/mspi_monthly_1950-2025.nc"
mspei_nc_file <- "mspei_results/mspei_monthly_1950-2025.nc"
basin_path <- "Spatial/nechakoBound_dissolve.shp"

#---- Configuration ----
drought_threshold_onset <- -1.0
drought_threshold_end <- 0.0
severe_threshold <- -1.3
extreme_threshold <- -1.6
exceptional_threshold <- -2.0
min_duration <- 2  # months

#---- DIAGNOSTIC: Check file availability ----
cat("\n--- Checking MSPI/MSPEI data files ---\n")
if (!file.exists(mspi_file)) {
  stop(sprintf("MSPI file not found: %s\n  Please run MSPI calculation script first.", mspi_file))
}
if (!file.exists(mspei_file)) {
  stop(sprintf("MSPEI file not found: %s\n  Please run MSPEI calculation script first.", mspei_file))
}
cat(sprintf("✓ MSPI file found: %s\n", basename(mspi_file)))
cat(sprintf("✓ MSPEI file found: %s\n", basename(mspei_file)))

# NEW: Check NetCDF files for spatial maps
cat("\n--- Checking NetCDF files for spatial maps ---\n")
if (file.exists(mspi_nc_file)) {
  cat(sprintf("✓ MSPI NetCDF found: %s\n", basename(mspi_nc_file)))
} else {
  cat(sprintf("⚠ MSPI NetCDF not found: %s (spatial maps will be skipped)\n", basename(mspi_nc_file)))
}
if (file.exists(mspei_nc_file)) {
  cat(sprintf("✓ MSPEI NetCDF found: %s\n", basename(mspei_nc_file)))
} else {
  cat(sprintf("⚠ MSPEI NetCDF not found: %s (spatial maps will be skipped)\n", basename(mspei_nc_file)))
}

#---- Load and prepare data ----
cat("\n--- Loading data ---\n")
mspi <- read.xlsx(mspi_file, sheet = "Basin_Timeseries")
mspei <- read.xlsx(mspei_file, sheet = "Basin_Timeseries")

# Convert date strings to Date objects
mspi$Date_formatted <- as.Date(paste0(mspi$Date, "-01"))
mspei$Date_formatted <- as.Date(paste0(mspei$Date, "-01"))

cat(sprintf("✓ MSPI data loaded: %d records from %s to %s\n", 
            nrow(mspi), 
            format(min(mspi$Date_formatted), "%Y-%m"), 
            format(max(mspi$Date_formatted), "%Y-%m")))
cat(sprintf("✓ MSPEI data loaded: %d records from %s to %s\n", 
            nrow(mspei), 
            format(min(mspei$Date_formatted), "%Y-%m"), 
            format(max(mspei$Date_formatted), "%Y-%m")))

#---- Detect drought events (kept here for plot shading; authoritative version
#     with CSV/Excel output lives in MSPI_calculation.R) ----
detect_drought_events <- function(df, onset_threshold = -1.0, termination_threshold = 0.0, min_duration = 2) {
  df <- df[order(df$Date_formatted), ]
  values <- df$Basin_Mean
  dates <- df$Date_formatted
  
  events <- data.frame(
    event_id = integer(),
    start_date = as.Date(character()),
    end_date = as.Date(character()),
    duration_months = integer(),
    min_value = numeric(),
    severity = character(),
    stringsAsFactors = FALSE
  )
  
  event_id <- 0
  i <- 1
  
  while (i <= length(values)) {
    if (!is.na(values[i]) && values[i] < onset_threshold) {
      event_id <- event_id + 1
      start_idx <- i
      min_val <- values[i]
      
      j <- i + 1
      while (j <= length(values) && !is.na(values[j]) && values[j] < termination_threshold) {
        min_val <- min(min_val, values[j], na.rm = TRUE)
        j <- j + 1
      }
      
      end_idx <- j - 1
      duration <- end_idx - start_idx + 1
      
      if (duration >= min_duration) {
        if (min_val < exceptional_threshold) {
          severity <- "Exceptional"
        } else if (min_val < extreme_threshold) {
          severity <- "Extreme"
        } else if (min_val < severe_threshold) {
          severity <- "Severe"
        } else {
          severity <- "Moderate"
        }
        
        events <- rbind(events, data.frame(
          event_id = event_id,
          start_date = dates[start_idx],
          end_date = dates[end_idx],
          duration_months = duration,
          min_value = round(min_val, 3),
          severity = severity,
          stringsAsFactors = FALSE
        ))
      }
      
      i <- end_idx + 1
    } else {
      i <- i + 1
    }
  }
  
  return(events)
}

##############################################
# PLOTTING FUNCTIONS (matching SPI/SPEI style)
##############################################

#---- Create full-period plot with drought shading (1950-2025) ----
create_timeseries_plot_full <- function(df, index_type) {
  events <- detect_drought_events(df, drought_threshold_onset, drought_threshold_end, min_duration)
  
  index_label <- toupper(index_type)
  title_text <- sprintf("Nechako Basin %s Basin-Wide Mean (1950-2025)", index_label)
  
  p <- ggplot(df, aes(x = Date_formatted, y = Basin_Mean)) +
    geom_line(color = if (index_type == "mspi") "#1f78b4" else "#e31a1c", linewidth = 0.9) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.7)
  
  if (nrow(events) > 0) {
    for (i in seq_len(nrow(events))) {
      shade_color <- switch(events$severity[i],
                            "Exceptional" = "purple",
                            "Extreme" = "darkred",
                            "Severe" = "red",
                            "Moderate" = "orange")
      
      p <- p + annotate("rect",
                        xmin = events$start_date[i], 
                        xmax = events$end_date[i],
                        ymin = -Inf, ymax = Inf,
                        fill = shade_color, alpha = 0.12)
    }
  }
  
  p <- p +
    geom_hline(yintercept = drought_threshold_onset, linetype = "dashed", color = "orange", linewidth = 0.8) +
    geom_hline(yintercept = severe_threshold, linetype = "dashed", color = "red", linewidth = 0.8) +
    geom_hline(yintercept = extreme_threshold, linetype = "dashed", color = "darkred", linewidth = 0.8) +
    annotate("text", x = as.Date("1955-01-01"), y = drought_threshold_onset, 
             label = "Moderate", hjust = 0, vjust = -0.5, size = 3.5, color = "orange", fontface = "bold") +
    annotate("text", x = as.Date("1955-01-01"), y = severe_threshold, 
             label = "Severe", hjust = 0, vjust = -0.5, size = 3.5, color = "red", fontface = "bold") +
    annotate("text", x = as.Date("1955-01-01"), y = extreme_threshold, 
             label = "Extreme", hjust = 0, vjust = -0.5, size = 3.5, color = "darkred", fontface = "bold") +
    labs(title = title_text,
         subtitle = sprintf("Drought thresholds: Moderate (< %.1f) | Severe (< %.1f) | Extreme (< %.1f) | %d drought events detected (≥%d months)", 
                            drought_threshold_onset, severe_threshold, extreme_threshold, 
                            nrow(events), min_duration),
         x = "Year", y = "Index Value") +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    ) +
    scale_x_date(date_breaks = "5 years", date_labels = "%Y", expand = expansion(mult = c(0.01, 0.01))) +
    coord_cartesian(ylim = c(-3.5, 3.5))
  
  return(p)
}

#---- Create recent-period plot with drought shading (2020-2025) ----
create_timeseries_plot_recent <- function(df, index_type) {
  df_recent <- df[df$Date_formatted >= as.Date("2020-01-01") & df$Date_formatted <= as.Date("2025-12-31"), ]
  events <- detect_drought_events(df_recent, drought_threshold_onset, drought_threshold_end, min_duration)
  
  index_label <- toupper(index_type)
  title_text <- sprintf("Nechako Basin %s Drought Events (2020-2025)", index_label)
  
  p <- ggplot(df_recent, aes(x = Date_formatted, y = Basin_Mean)) +
    geom_line(color = if (index_type == "mspi") "#1f78b4" else "#e31a1c", linewidth = 1.1) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.7) +
    geom_hline(yintercept = drought_threshold_onset, linetype = "dashed", color = "orange", linewidth = 0.9) +
    geom_hline(yintercept = severe_threshold, linetype = "dashed", color = "red", linewidth = 0.9) +
    geom_hline(yintercept = extreme_threshold, linetype = "dashed", color = "darkred", linewidth = 0.9) +
    labs(title = title_text,
         subtitle = sprintf("Drought definition: Onset < %.1f | Termination ≥ %.1f", 
                            drought_threshold_onset, drought_threshold_end),
         x = "Year", y = "Index Value",
         caption = "Shaded periods indicate drought events") +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30"),
      plot.caption = element_text(hjust = 0.5, size = 10, face = "italic", color = "gray40"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    ) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = expansion(mult = c(0.02, 0.02))) +
    coord_cartesian(ylim = c(-3.5, 3.5))
  
  if (nrow(events) > 0) {
    for (i in seq_len(nrow(events))) {
      shade_color <- switch(events$severity[i],
                            "Exceptional" = "purple",
                            "Extreme" = "darkred",
                            "Severe" = "red",
                            "Moderate" = "orange")
      
      p <- p + annotate("rect",
                        xmin = events$start_date[i], xmax = events$end_date[i],
                        ymin = -Inf, ymax = Inf,
                        fill = shade_color, alpha = 0.18)
      
      mid_date <- events$start_date[i] + (events$end_date[i] - events$start_date[i]) / 2
      label_y <- 3.2
      
      p <- p + annotate("text", 
                        x = mid_date, y = label_y,
                        label = sprintf("%s\n(%d mo)", events$severity[i], events$duration_months[i]),
                        size = 3.3,
                        color = shade_color,
                        fontface = "bold",
                        hjust = 0.5, vjust = 0.5)
    }
  }
  
  return(p)
}

#---- Create MSPI vs MSPEI comparison plot ----
create_comparison_plot <- function(mspi_df, mspei_df) {
  merged <- merge(mspi_df[, c("Date_formatted", "Basin_Mean")], 
                  mspei_df[, c("Date_formatted", "Basin_Mean")],
                  by = "Date_formatted", suffixes = c("_MSPI", "_MSPEI"))
  
  merged_long <- data.frame(
    Date = rep(merged$Date_formatted, 2),
    Index = rep(c("MSPI", "MSPEI"), each = nrow(merged)),
    Value = c(merged$Basin_Mean_MSPI, merged$Basin_Mean_MSPEI)
  )
  
  corr_val <- cor(merged$Basin_Mean_MSPI, merged$Basin_Mean_MSPEI, use = "complete.obs")
  
  mspi_events <- detect_drought_events(mspi_df, drought_threshold_onset, drought_threshold_end, min_duration)
  mspei_events <- detect_drought_events(mspei_df, drought_threshold_onset, drought_threshold_end, min_duration)
  
  title_text <- "Nechako Basin MSPI vs MSPEI Comparison (1950-2025)"
  
  p <- ggplot(merged_long, aes(x = Date, y = Value, color = Index)) +
    geom_line(linewidth = 0.75, alpha = 0.85) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.6) +
    geom_hline(yintercept = drought_threshold_onset, linetype = "dashed", color = "gray50", linewidth = 0.5, alpha = 0.7) +
    geom_hline(yintercept = severe_threshold, linetype = "dashed", color = "gray50", linewidth = 0.5, alpha = 0.7) +
    geom_hline(yintercept = extreme_threshold, linetype = "dashed", color = "gray50", linewidth = 0.5, alpha = 0.7) +
    scale_color_manual(values = c("MSPI" = "#1f78b4", "MSPEI" = "#e31a1c")) +
    theme_minimal(base_size = 13) +
    labs(
      title = title_text,
      subtitle = sprintf("Correlation: r = %.3f | MSPI events: %d | MSPEI events: %d", 
                         corr_val, nrow(mspi_events), nrow(mspei_events)),
      x = "Year",
      y = "Index Value",
      color = "Index Type",
      caption = "MSPI = precipitation only | MSPEI = precipitation minus evapotranspiration"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30"),
      axis.title = element_text(face = "bold"),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 11),
      panel.grid.minor = element_blank(),
      plot.caption = element_text(size = 9, hjust = 0.5, face = "italic", color = "gray40")
    ) +
    scale_x_date(date_breaks = "5 years", date_labels = "%Y", expand = expansion(mult = c(0.01, 0.01))) +
    coord_cartesian(ylim = c(-3.5, 3.5))
  
  return(p)
}

##############################################
# NEW: SPATIAL MAP FUNCTIONS
##############################################
#---- Create single spatial map ----
create_spatial_map <- function(rast_obj, layer_idx, date_val, basin_vect, index_type) {
  
  layer_rast <- rast_obj[[layer_idx]]
  
  # Crop to basin extent first
  layer_rast <- crop(layer_rast, basin_vect)
  
  # Mask with touches=TRUE (includes edge pixels)
  layer_rast <- mask(layer_rast, basin_vect, touches = TRUE)
  
  # Convert to data frame, remove NAs
  df <- as.data.frame(layer_rast, xy = TRUE, na.rm = TRUE)
  colnames(df) <- c("x", "y", "value")
  
  # Convert basin to sf for overlay
  basin_sf <- sf::st_as_sf(basin_vect)
  basin_sf_crop <- sf::st_crop(basin_sf, sf::st_bbox(layer_rast))
  
  # Calculate statistics on raw values (before any clamping)
  valid_vals <- df$value
  n_valid    <- nrow(df)
  val_mean   <- mean(valid_vals, na.rm = TRUE)
  val_sd     <- sd(valid_vals, na.rm = TRUE)
  val_min    <- min(valid_vals, na.rm = TRUE)
  val_max    <- max(valid_vals, na.rm = TRUE)
  
  # --- FIX: Dynamic color limits ---
  # Always show at least [-3, 3]; expand symmetrically if data goes beyond.
  # Using oob = scales::squish ensures out-of-range values are shown at the
  # extreme colour rather than turning gray (the default censoring behaviour).
  abs_max   <- max(3, ceiling(max(abs(val_min), abs(val_max))))
  clim      <- c(-abs_max, abs_max)
  
  # Count pixels clipped at the limits so the user is informed
  n_above <- sum(valid_vals > 3,  na.rm = TRUE)
  n_below <- sum(valid_vals < -3, na.rm = TRUE)
  clip_note <- if (n_above + n_below > 0) {
    sprintf(" | %d px > 3, %d px < -3 (squished to limit)", n_above, n_below)
  } else ""
  
  # Create plot
  p <- ggplot(df, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    geom_sf(data = basin_sf_crop, fill = NA, color = "black", linewidth = 0.8, 
            inherit.aes = FALSE) +
    # KEY FIX: oob = scales::squish prevents out-of-range values going gray
    scale_fill_gradient2(
      low      = "#d73027",
      mid      = "#ffffbf",
      high     = "#1a9850",
      midpoint = 0,
      limits   = clim,
      oob      = scales::squish,   # squish to colour limits instead of NA
      name     = "Index Value",
      breaks   = seq(-abs_max, abs_max, by = ifelse(abs_max > 3, 1, 1)),
      labels   = function(x) sprintf("%.0f", x)
    ) +
    labs(
      title   = sprintf("%s \u2013 %s", toupper(index_type), format(date_val, "%B %Y")),
      x       = "Longitude",
      y       = "Latitude",
      caption = sprintf("Valid pixels: %d | Mean: %.2f | SD: %.2f%s",
                        n_valid, val_mean, val_sd, clip_note)
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title   = element_text(face = "bold", size = 14, hjust = 0.5),
      panel.grid   = element_blank(),
      axis.text    = element_text(size = 10),
      axis.title   = element_text(face = "bold", size = 11),
      legend.position  = "right",
      legend.title     = element_text(face = "bold", size = 11),
      plot.caption     = element_text(size = 9, hjust = 0.5, face = "italic", 
                                      color = "gray40")
    ) +
    coord_sf()
  
  return(p)
}
#---- Select 20 months from Jan 2022-Dec 2025 ----
select_mapping_months <- function(dates_vector) {
  start_date <- as.Date("2022-01-01")
  end_date <- as.Date("2025-12-01")
  period_indices <- which(dates_vector >= start_date & dates_vector <= end_date)
  
  if (length(period_indices) < 20) {
    warning(sprintf("Not enough data in period Jan 2022-Dec 2025. Found: %d months", length(period_indices)))
    return(list(indices = period_indices, dates = dates_vector[period_indices]))
  }
  
  jan_2022_idx <- which(dates_vector == as.Date("2022-01-01"))
  dec_2025_idx <- which(dates_vector == as.Date("2025-12-01"))
  
  if (length(jan_2022_idx) == 0) jan_2022_idx <- period_indices[1]
  if (length(dec_2025_idx) == 0) dec_2025_idx <- period_indices[length(period_indices)]
  
  remaining_indices <- setdiff(period_indices, c(jan_2022_idx, dec_2025_idx))
  set.seed(123)
  random_indices <- sample(remaining_indices, size = min(18, length(remaining_indices)))
  
  selected_indices <- sort(c(jan_2022_idx, dec_2025_idx, random_indices))
  selected_dates <- dates_vector[selected_indices]
  
  return(list(indices = selected_indices, dates = selected_dates))
}

##############################################
# MAIN PROCESSING
##############################################

cat("\n============================================================\n")
cat("BASIN-AVERAGED MSPI & MSPEI VISUALIZATION\n")
cat("DROUGHT DEFINITION: Onset < -1.0, Termination >= 0.0\n")
cat(sprintf("MINIMUM DURATION: %d months\n", min_duration))
cat("============================================================\n")

while (dev.cur() > 1) try(dev.off(), silent = TRUE)

for (index_type in c("mspi", "mspei")) {
  index_label <- toupper(index_type)
  cat(sprintf("\n>>>>> Processing %s <<<<<\n", index_label))
  
  df <- if (index_type == "mspi") mspi else mspei
  
  cat(sprintf("  - Data range: %s to %s (n=%d)\n", 
              format(min(df$Date_formatted), "%Y-%m"), 
              format(max(df$Date_formatted), "%Y-%m"), 
              nrow(df)))
  
  write.csv(df, file.path(mspi_output_dir, sprintf("%s_basin_average_1950_2025.csv", index_type)), 
            row.names = FALSE)
  cat(sprintf("  + Saved: %s_basin_average_1950_2025.csv\n", index_type))
  
  p_full <- create_timeseries_plot_full(df, index_type)
  pdf_file <- file.path(mspi_output_dir, sprintf("%s_full_period_1950_2025.pdf", toupper(index_type)))
  
  if (safe_pdf(pdf_file, 12, 8)) {
    print(p_full)
    dev.off()
    cat(sprintf("  + Saved: %s\n", basename(pdf_file)))
  }
  
  df_recent <- df[df$Date_formatted >= as.Date("2020-01-01") & df$Date_formatted <= as.Date("2025-12-31"), ]
  events <- detect_drought_events(df_recent, drought_threshold_onset, drought_threshold_end, min_duration)
  
  if (nrow(events) > 0) {
    write.csv(events, file.path(mspi_output_dir, sprintf("%s_drought_events_2020_2025.csv", index_type)), 
              row.names = FALSE)
    cat(sprintf("  + Detected %d drought events (min duration %d months)\n", nrow(events), min_duration))
  } else {
    cat(sprintf("  + No drought events detected in 2020-2025 period (min duration %d months)\n", min_duration))
  }
  
  p_recent <- create_timeseries_plot_recent(df, index_type)
  pdf_file <- file.path(mspi_output_dir, sprintf("%s_recent_period_2020_2025.pdf", toupper(index_type)))
  
  if (safe_pdf(pdf_file, 12, 8)) {
    print(p_recent)
    dev.off()
    cat(sprintf("  + Saved: %s\n", basename(pdf_file)))
  }
  
  while (dev.cur() > 1) try(dev.off(), silent = TRUE)
}

cat("\n>>>>> Creating MSPI vs MSPEI Comparison Plot <<<<<\n")
p_comparison <- create_comparison_plot(mspi, mspei)
pdf_file <- file.path(mspi_output_dir, "MSPI_MSPEI_comparison_1950_2025.pdf")

if (safe_pdf(pdf_file, 14, 9)) {
  print(p_comparison)
  dev.off()
  cat(sprintf("  + Saved: %s\n", basename(pdf_file)))
}

##############################################
# NEW: SPATIAL MAPS SECTION
##############################################

cat("\n============================================================\n")
cat("SPATIAL MAP VISUALIZATION (Jan 2022 - Dec 2025)\n")
cat("============================================================\n")

if (file.exists(mspi_nc_file) && file.exists(mspei_nc_file) && file.exists(basin_path)) {
  
  cat("\n--- Loading NetCDF data ---\n")
  mspi_rast <- rast(mspi_nc_file)
  mspei_rast <- rast(mspei_nc_file)
  
  cat(sprintf("✓ MSPI raster loaded: %d layers\n", nlyr(mspi_rast)))
  cat(sprintf("✓ MSPEI raster loaded: %d layers\n", nlyr(mspei_rast)))
  
  mspi_times <- time(mspi_rast)
  mspei_times <- time(mspei_rast)
  
  if (is.null(mspi_times)) {
    mspi_times <- seq(as.Date("1950-01-01"), by = "month", length.out = nlyr(mspi_rast))
    mspei_times <- seq(as.Date("1950-01-01"), by = "month", length.out = nlyr(mspei_rast))
  }
  
  basin <- vect(basin_path)
  if (!same.crs(mspi_rast, basin)) {
    basin <- project(basin, crs(mspi_rast))
  }
  
  cat("\n--- Selecting 20 months for mapping ---\n")
  month_selection <- select_mapping_months(mspi_times)
  selected_indices <- month_selection$indices
  selected_dates <- month_selection$dates
  
  cat(sprintf("✓ Selected %d months from Jan 2022-Dec 2025\n", length(selected_indices)))
  cat(sprintf("✓ Guaranteed: %s and %s\n", 
              format(as.Date("2022-01-01"), "%Y-%m"), 
              format(as.Date("2025-12-01"), "%Y-%m")))
  
  cat("\n--- Creating MSPI Spatial Maps ---\n")
  
  mspi_pdf_file <- file.path(mspi_output_dir, "mspi_spatial_maps_2022_2025.pdf")
  pdf(mspi_pdf_file, width = 10, height = 8)
  
  for (i in seq_along(selected_indices)) {
    idx <- selected_indices[i]
    date_val <- selected_dates[i]
    
    cat(sprintf("  Processing MSPI: %s (layer %d/%d)\n", 
                format(date_val, "%Y-%m"), i, length(selected_indices)))
    
    p <- create_spatial_map(mspi_rast, idx, date_val, basin, "MSPI")
    print(p)
    
    if (i < length(selected_indices)) {
      grid::grid.newpage()
    }
  }
  
  dev.off()
  cat(sprintf("✓ Saved MSPI maps: %s\n", basename(mspi_pdf_file)))
  
  cat("\n--- Creating MSPEI Spatial Maps ---\n")
  
  mspei_pdf_file <- file.path(mspei_output_dir, "mspei_spatial_maps_2022_2025.pdf")
  pdf(mspei_pdf_file, width = 10, height = 8)
  
  for (i in seq_along(selected_indices)) {
    idx <- selected_indices[i]
    date_val <- selected_dates[i]
    
    cat(sprintf("  Processing MSPEI: %s (layer %d/%d)\n", 
                format(date_val, "%Y-%m"), i, length(selected_indices)))
    
    p <- create_spatial_map(mspei_rast, idx, date_val, basin, "MSPEI")
    print(p)
    
    if (i < length(selected_indices)) {
      grid::grid.newpage()
    }
  }
  
  dev.off()
  cat(sprintf("✓ Saved MSPEI maps: %s\n", basename(mspei_pdf_file)))
  
  selected_months_df <- data.frame(
    Index = 1:length(selected_indices),
    Date = format(selected_dates, "%Y-%m"),
    Guaranteed = ifelse(selected_dates %in% c(as.Date("2022-01-01"), as.Date("2025-12-01")), "Yes", "No")
  )
  
  write.csv(selected_months_df, 
            file.path(mspi_output_dir, "selected_months_for_mapping.csv"), 
            row.names = FALSE)
  cat(sprintf("✓ Saved selected months list: selected_months_for_mapping.csv\n"))
  
} else {
  cat("\n⚠ Skipping spatial maps (NetCDF files or basin boundary not found)\n")
}

##############################################
# SUMMARY STATISTICS
# The authoritative Summary Statistics Excel (including full drought catalogue)
# is produced by MSPI_calculation.R as MSPI_MSPEI_Summary_Statistics_YYYY-YYYY.xlsx.
# Here we just print a quick console summary from the already-loaded data.
##############################################

cat("\n--- Quick Drought Summary (from loaded data) ---\n")

mspi_droughts_full   <- detect_drought_events(mspi,  drought_threshold_onset,
                                              drought_threshold_end, min_duration)
mspei_droughts_full  <- detect_drought_events(mspei, drought_threshold_onset,
                                              drought_threshold_end, min_duration)
mspi_recent  <- mspi [mspi $Date_formatted >= as.Date("2020-01-01"), ]
mspei_recent <- mspei[mspei$Date_formatted >= as.Date("2020-01-01"), ]
mspi_droughts_recent  <- detect_drought_events(mspi_recent,  drought_threshold_onset,
                                               drought_threshold_end, min_duration)
mspei_droughts_recent <- detect_drought_events(mspei_recent, drought_threshold_onset,
                                               drought_threshold_end, min_duration)

##############################################
# FINAL SUMMARY
##############################################

cat("\n*** ALL VISUALIZATIONS COMPLETE ***\n")
cat(sprintf("Output directory (MSPI): %s\n", normalizePath(mspi_output_dir)))
cat(sprintf("Output directory (MSPEI): %s\n", normalizePath(mspei_output_dir)))
cat(sprintf("Generated:\n"))
cat(sprintf("  - 5 PDF files (2 indices x 2 periods + 1 comparison)\n"))
cat(sprintf("  - 2 basin-average CSV files\n"))
cat(sprintf("  - 2 PDF spatial map files (20 months each)\n"))
cat(sprintf("  - 1 CSV selected months list\n"))
cat(sprintf("  NOTE: Drought catalogue + summary statistics Excel are produced\n"))
cat(sprintf("        by MSPI_calculation.R (MSPI_MSPEI_Summary_Statistics_*.xlsx)\n"))
cat(sprintf("\nDrought Event Summary:\n"))
cat(sprintf("  MSPI  - Full period: %d events | Recent (2020-present): %d events\n", 
            nrow(mspi_droughts_full), nrow(mspi_droughts_recent)))
cat(sprintf("  MSPEI - Full period: %d events | Recent (2020-present): %d events\n", 
            nrow(mspei_droughts_full), nrow(mspei_droughts_recent)))