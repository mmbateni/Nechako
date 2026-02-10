##############################################
# BASIN-AVERAGED SPI & SPEI VISUALIZATION
# - CORRECT TIME PARSING for flexible date formats (1970-1-1, 1950-01-01, etc.)
# - Generates SEPARATE PDFs per timescale (no multi-panel grids)
# - Drought definition: onset < -1.0, termination >= 0.0
# - Shaded drought periods with severity labels
# - Titles EMBEDDED in plots (PDF-safe)
# - Timescales: 1, 3, 6, 9, 12 months
# - **NEW**: SPI vs SPEI comparison plots for each timescale with correlations
# - **NEW**: Correlation statistics in Excel summary file
##############################################

library(terra)
library(ncdf4)
library(ggplot2)
library(grid)
library(lubridate)
library(openxlsx)

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
setwd("D:/Nechako_Drought/")
spi_dir <- "spi_results_seasonal"
spei_dir <- "spei_results_seasonal"
plot_dir <- "basin_averaged_plots"
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

basin_path <- "Spatial/nechakoBound_dissolve.shp"
basin <- vect(basin_path)
cat("+ Basin boundary loaded\n")

# Reproject basin to BC Albers Equal Area
target_crs <- "EPSG:3005"
if (!same.crs(basin, target_crs)) {
  basin <- project(basin, target_crs)
  cat("+ Basin reprojected to BC Albers (EPSG:3005)\n")
}

#---- Configuration ----
timescales <- c(1, 3, 6, 9, 12)
drought_threshold_onset <- -1.0
drought_threshold_end <- 0.0
severe_threshold <- -1.3
extreme_threshold <- -1.6

#---- DIAGNOSTIC: Scan directories ----
cat("\n--- Scanning NetCDF file structure ---\n")
for (dir_path in c(spi_dir, spei_dir)) {
  if (!dir.exists(dir_path)) {
    cat(sprintf("WARNING: Directory '%s' does not exist!\n", dir_path))
    next
  }
  all_nc <- list.files(dir_path, pattern = "\\.nc$", full.names = TRUE, recursive = TRUE)
  cat(sprintf("\n%s (%d files):\n", dir_path, length(all_nc)))
  if (length(all_nc) > 0) {
    samples <- basename(all_nc[1:min(5, length(all_nc))])
    for (s in samples) cat(sprintf("  - %s\n", s))
  }
}

#---- CRITICAL FIX: Robust time parsing for ALL CF-convention formats ----
parse_time_from_nc <- function(nc_file) {
  nc <- nc_open(nc_file)
  time_vals <- ncvar_get(nc, "time")
  time_units <- ncatt_get(nc, "time", "units")$value
  nc_close(nc)
  
  cat(sprintf("    Time units: '%s'\n", time_units))
  
  # EXTRACT ORIGIN DATE WITH FLEXIBLE REGEX (handles 1970-1-1, 1950-01-01, etc.)
  origin_match <- regmatches(time_units, regexpr("\\d{4}-\\d{1,2}-\\d{1,2}", time_units))
  if (length(origin_match) == 0) {
    # SECONDARY ATTEMPT: Try space-separated format "1970 1 1"
    origin_match <- regmatches(time_units, regexpr("\\d{4} \\d{1,2} \\d{1,2}", time_units))
    if (length(origin_match) > 0) {
      # Convert space to dash for as.Date()
      origin_match <- gsub(" ", "-", origin_match)
    } else {
      stop(sprintf("Cannot parse origin date from time units: '%s'", time_units))
    }
  }
  
  # Convert to Date (handles both single and double-digit months/days)
  origin_date <- as.Date(origin_match, format = "%Y-%m-%d")
  
  # Handle time units conversion
  if (grepl("hours", time_units, ignore.case = TRUE)) {
    time_vals <- time_vals / 24
    cat(sprintf("    ⚠ Converted hours to days (origin: %s)\n", origin_date))
  } else if (!grepl("days", time_units, ignore.case = TRUE)) {
    cat(sprintf("    ⚠ Warning: Unknown time unit in '%s' - assuming days\n", time_units))
  }
  
  dates <- origin_date + time_vals
  
  # CRITICAL: Adjust for common NetCDF metadata error (1970 origin instead of 1950)
  # If dates start near 1970 but should be 1950, shift by 20 years
  if (year(min(dates)) >= 1968 && year(min(dates)) <= 1972) {
    # Likely metadata error - shift to 1950 start
    shift_years <- 1950 - year(min(dates))
    if (abs(shift_years) >= 15) {  # Only shift if significant difference
      dates <- dates + years(shift_years)
      cat(sprintf("    ⚠ Adjusted date origin: shifted %d years (metadata correction)\n", shift_years))
      cat(sprintf("    Adjusted time range: %s to %s\n", 
                  format(min(dates), "%Y-%m"), format(max(dates), "%Y-%m")))
    }
  }
  
  cat(sprintf("    Time range: %s to %s (%d time steps)\n", 
              format(min(dates), "%Y-%m"), format(max(dates), "%Y-%m"), length(dates)))
  
  list(dates = dates, origin = origin_date, units = time_units)
}

#---- REVISED: Basin average calculation with CORRECT time parsing ----
calculate_basin_average <- function(index_type, scale) {
  dir_path <- if (index_type == "spi") spi_dir else spei_dir
  
  # Match your actual file pattern: spi_01_month01_Jan.nc, etc.
  pattern <- sprintf("%s_%02d_month\\d{2}_.*\\.nc$", index_type, scale)
  files <- list.files(dir_path, pattern = pattern, full.names = TRUE, 
                      recursive = TRUE, ignore.case = TRUE)
  
  if (length(files) == 0) {
    stop(sprintf("No %s-%02d files found matching pattern '%s' in %s", 
                 toupper(index_type), scale, pattern, dir_path))
  }
  
  cat(sprintf("    ✓ Found %d file(s) for %s-%02d\n", length(files), toupper(index_type), scale))
  
  # Parse time from FIRST file (all files should share same time dimension)
  time_info <- parse_time_from_nc(files[1])
  dates <- time_info$dates
  
  # Process each time step - average across all monthly files for this timescale
  avg_values <- numeric(length(dates))
  for (t in seq_along(dates)) {
    r_list <- list()
    for (f in files) {
      tryCatch({
        r <- rast(f, lyr = t)
        if (!is.na(crs(r)) && !same.crs(r, target_crs)) {
          r <- project(r, target_crs)
        }
        r_list[[length(r_list) + 1]] <- r
      }, error = function(e) {
        cat(sprintf("    ⚠ Warning: Skipping layer %d in %s: %s\n", t, basename(f), e$message))
      })
    }
    
    if (length(r_list) == 0) {
      avg_values[t] <- NA
    } else {
      stack_r <- rast(r_list)
      avg_r <- mean(stack_r, na.rm = TRUE)
      extracted <- extract(avg_r, basin, fun = mean, na.rm = TRUE)
      avg_values[t] <- if (!is.null(extracted) && nrow(extracted) > 0) extracted[1, 1] else NA
    }
  }
  
  data.frame(date = dates, value = avg_values, stringsAsFactors = FALSE)
}

#---- Detect drought events ----
detect_drought_events <- function(df, onset_threshold = -1.0, termination_threshold = 0.0, min_duration = 1) {
  df <- df[order(df$date), ]
  values <- df$value
  dates <- df$date
  
  events <- data.frame(
    event_id = integer(),
    start_date = as.Date(character()),
    end_date = as.Date(character()),
    duration_months = integer(),
    min_value = numeric(),
    severity = character(),
    stringsAsFactors = FALSE
  )
  
  in_drought <- FALSE
  start_idx <- NA
  
  for (i in seq_along(values)) {
    if (!in_drought && !is.na(values[i]) && values[i] < onset_threshold) {
      in_drought <- TRUE
      start_idx <- i
    } else if (in_drought && !is.na(values[i]) && values[i] >= termination_threshold) {
      end_idx <- i - 1
      duration <- end_idx - start_idx + 1
      
      if (duration >= min_duration) {
        min_val <- min(values[start_idx:end_idx], na.rm = TRUE)
        severity <- if (min_val < extreme_threshold) "Extreme" else 
          if (min_val < severe_threshold) "Severe" else "Moderate"
        
        events <- rbind(events, data.frame(
          event_id = nrow(events) + 1,
          start_date = dates[start_idx],
          end_date = dates[end_idx],
          duration_months = duration,
          min_value = min_val,
          severity = severity
        ))
      }
      
      in_drought <- FALSE
      start_idx <- NA
    }
  }
  
  # Handle drought continuing to end of record
  if (in_drought && !is.na(start_idx)) {
    end_idx <- length(values)
    duration <- end_idx - start_idx + 1
    
    if (duration >= min_duration) {
      min_val <- min(values[start_idx:end_idx], na.rm = TRUE)
      severity <- if (min_val < extreme_threshold) "Extreme" else 
        if (min_val < severe_threshold) "Severe" else "Moderate"
      
      events <- rbind(events, data.frame(
        event_id = nrow(events) + 1,
        start_date = dates[start_idx],
        end_date = dates[end_idx],
        duration_months = duration,
        min_value = min_val,
        severity = severity
      ))
    }
  }
  
  rownames(events) <- NULL
  events
}

#---- Create full-period timeseries plot (1950-2025) ----
create_timeseries_plot_full <- function(df, index_type, scale) {
  index_label <- if (index_type == "spi") "SPI" else "SPEI"
  title_text <- sprintf("Nechako Basin %s-%02d (1950-2025)", index_label, scale)
  
  ggplot(df, aes(x = date, y = value)) +
    geom_line(color = if (index_type == "spi") "#1f78b4" else "#e31a1c", size = 0.9) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.6) +
    geom_hline(yintercept = drought_threshold_onset, linetype = "dashed", color = "orange", size = 0.8) +
    geom_hline(yintercept = severe_threshold, linetype = "dashed", color = "red", size = 0.8) +
    geom_hline(yintercept = extreme_threshold, linetype = "dashed", color = "darkred", size = 0.8) +
    annotate("text", x = as.Date("1955-01-01"), y = drought_threshold_onset, 
             label = "Moderate", hjust = 0, vjust = -0.5, size = 3.5, color = "orange", fontface = "bold") +
    annotate("text", x = as.Date("1955-01-01"), y = severe_threshold, 
             label = "Severe", hjust = 0, vjust = -0.5, size = 3.5, color = "red", fontface = "bold") +
    annotate("text", x = as.Date("1955-01-01"), y = extreme_threshold, 
             label = "Extreme", hjust = 0, vjust = -0.5, size = 3.5, color = "darkred", fontface = "bold") +
    labs(title = title_text,
         subtitle = sprintf("Drought thresholds: Moderate (< %.1f) | Severe (< %.1f) | Extreme (< %.1f)", 
                            drought_threshold_onset, severe_threshold, extreme_threshold),
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
}

#---- Create recent-period plot with drought shading (2020-2025) ----
create_timeseries_plot_recent <- function(df, index_type, scale) {
  df_recent <- df[df$date >= as.Date("2020-01-01") & df$date <= as.Date("2025-12-31"), ]
  events <- detect_drought_events(df_recent, drought_threshold_onset, drought_threshold_end, 1)
  
  index_label <- if (index_type == "spi") "SPI" else "SPEI"
  title_text <- sprintf("Nechako Basin %s-%02d Drought Events (2020-2025)", index_label, scale)
  
  p <- ggplot(df_recent, aes(x = date, y = value)) +
    geom_line(color = if (index_type == "spi") "#1f78b4" else "#e31a1c", size = 1.1) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.7) +
    geom_hline(yintercept = drought_threshold_onset, linetype = "dashed", color = "orange", size = 0.9) +
    geom_hline(yintercept = severe_threshold, linetype = "dashed", color = "red", size = 0.9) +
    geom_hline(yintercept = extreme_threshold, linetype = "dashed", color = "darkred", size = 0.9) +
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
  
  # Add drought shading and severity labels
  if (nrow(events) > 0) {
    for (i in seq_len(nrow(events))) {
      p <- p + annotate("rect",
                        xmin = events$start_date[i], xmax = events$end_date[i],
                        ymin = -Inf, ymax = Inf,
                        fill = if (events$severity[i] == "Extreme") "darkred" else 
                          if (events$severity[i] == "Severe") "red" else "orange",
                        alpha = 0.18)
      
      mid_date <- events$start_date[i] + (events$end_date[i] - events$start_date[i]) / 2
      label_y <- 3.2
      
      p <- p + annotate("text", 
                        x = mid_date, y = label_y,
                        label = sprintf("%s\n(%d mo)", events$severity[i], events$duration_months[i]),
                        size = 3.3,
                        color = if (events$severity[i] == "Extreme") "darkred" else 
                          if (events$severity[i] == "Severe") "red" else "orange",
                        fontface = "bold",
                        hjust = 0.5, vjust = 0.5)
    }
  }
  
  p
}

#---- Create comparison plot between SPI and SPEI for same timescale ----
create_comparison_plot <- function(spi_df, spei_df, scale) {
  # Rename columns for merging
  spi_df_merge <- spi_df
  names(spi_df_merge)[names(spi_df_merge) == "value"] <- "value_spi"
  
  spei_df_merge <- spei_df
  names(spei_df_merge)[names(spei_df_merge) == "value"] <- "value_spei"
  
  # Merge on date
  merged <- merge(spi_df_merge, spei_df_merge, by = "date", all = TRUE)
  
  # Reshape for ggplot
  merged_long <- data.frame(
    Date = rep(merged$date, 2),
    Index = rep(c(sprintf("SPI-%02d", scale), sprintf("SPEI-%02d", scale)), each = nrow(merged)),
    Value = c(merged$value_spi, merged$value_spei)
  )
  
  # Calculate correlation
  corr_val <- cor(merged$value_spi, merged$value_spei, use = "complete.obs")
  
  # Detect drought events for both indices
  spi_events <- detect_drought_events(spi_df, drought_threshold_onset, drought_threshold_end, 1)
  spei_events <- detect_drought_events(spei_df, drought_threshold_onset, drought_threshold_end, 1)
  
  title_text <- sprintf("Nechako Basin SPI-%02d vs SPEI-%02d Comparison (1950-2025)", scale, scale)
  
  # Create color vector with dynamic names
  color_vec <- c("#1f78b4", "#e31a1c")
  names(color_vec) <- c(sprintf("SPI-%02d", scale), sprintf("SPEI-%02d", scale))
  
  # Create plot
  p <- ggplot(merged_long, aes(x = Date, y = Value, color = Index)) +
    geom_line(size = 0.75, alpha = 0.85) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.6) +
    geom_hline(yintercept = drought_threshold_onset, linetype = "dashed", color = "gray50", size = 0.5, alpha = 0.7) +
    geom_hline(yintercept = severe_threshold, linetype = "dashed", color = "gray50", size = 0.5, alpha = 0.7) +
    geom_hline(yintercept = extreme_threshold, linetype = "dashed", color = "gray50", size = 0.5, alpha = 0.7) +
    scale_color_manual(values = color_vec) +
    theme_minimal(base_size = 13) +
    labs(
      title = title_text,
      subtitle = sprintf("Correlation: r = %.3f | SPI-%02d events: %d | SPEI-%02d events: %d", 
                         corr_val, scale, nrow(spi_events), scale, nrow(spei_events)),
      x = "Year",
      y = "Index Value",
      color = "Index Type",
      caption = sprintf("SPI-%02d = precipitation only | SPEI-%02d = precipitation minus evapotranspiration", 
                        scale, scale)
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

#---- MAIN PROCESSING ----
cat("\n============================================================\n")
cat("BASIN-AVERAGED SPI & SPEI VISUALIZATION (TIMESCALES: 1,3,6,9,12)\n")
cat("DROUGHT DEFINITION: Onset < -1.0, Termination >= 0.0\n")
cat("============================================================\n")

while (dev.cur() > 1) try(dev.off(), silent = TRUE)

for (index_type in c("spi", "spei")) {
  index_label <- toupper(index_type)
  cat(sprintf("\n>>>>> Processing %s <<<<<\n", index_label))
  
  for (scale in timescales) {
    cat(sprintf("  - Processing %s-%02d...\n", index_label, scale))
    
    df <- calculate_basin_average(index_type, scale)
    
    # Verify date range is correct before proceeding
    cat(sprintf("    Data range: %s to %s (n=%d)\n", 
                format(min(df$date), "%Y-%m"), format(max(df$date), "%Y-%m"), nrow(df)))
    
    # Safety check: abort if dates are clearly wrong
    if (min(df$date) < as.Date("1940-01-01") || max(df$date) > as.Date("2035-12-31")) {
      cat("    ⚠ WARNING: Date range seems incorrect - possible metadata error!\n")
    }
    
    write.csv(df, file.path(plot_dir, sprintf("%s_%02d_basin_average_1950_2025.csv", index_type, scale)), 
              row.names = FALSE)
    
    # FULL PERIOD (1950-2025)
    p_full <- create_timeseries_plot_full(df, index_type, scale)
    pdf_file <- file.path(plot_dir, sprintf("%s_%02d_full_period_1950_2025.pdf", index_type, scale))
    
    if (safe_pdf(pdf_file, 12, 8)) {
      print(p_full)
      dev.off()
      cat(sprintf("    + Saved: %s\n", basename(pdf_file)))
    }
    
    # RECENT PERIOD (2020-2025) with drought events
    df_recent <- df[df$date >= as.Date("2020-01-01") & df$date <= as.Date("2025-12-31"), ]
    events <- detect_drought_events(df_recent, drought_threshold_onset, drought_threshold_end, 1)
    
    if (nrow(events) > 0) {
      write.csv(events, file.path(plot_dir, sprintf("%s_%02d_drought_events_2020_2025.csv", index_type, scale)), 
                row.names = FALSE)
      cat(sprintf("    + Detected %d drought events (min duration 1 month)\n", nrow(events)))
    } else {
      cat("    + No drought events detected in 2020-2025 period\n")
    }
    
    p_recent <- create_timeseries_plot_recent(df, index_type, scale)
    pdf_file <- file.path(plot_dir, sprintf("%s_%02d_recent_period_2020_2025.pdf", index_type, scale))
    
    if (safe_pdf(pdf_file, 12, 8)) {
      print(p_recent)
      dev.off()
      cat(sprintf("    + Saved: %s\n", basename(pdf_file)))
    }
    
    while (dev.cur() > 1) try(dev.off(), silent = TRUE)
  }
}

#---- COMPARISON PLOTS: SPI vs SPEI for each timescale ----
cat("\n>>>>> Creating SPI vs SPEI Comparison Plots <<<<<\n")

for (scale in timescales) {
  cat(sprintf("  - Creating comparison plot for timescale %02d...\n", scale))
  
  # Calculate basin averages for both indices
  spi_df <- calculate_basin_average("spi", scale)
  spei_df <- calculate_basin_average("spei", scale)
  
  # Create comparison plot
  p_comparison <- create_comparison_plot(spi_df, spei_df, scale)
  pdf_file <- file.path(plot_dir, sprintf("SPI_SPEI_comparison_%02d_month_1950_2025.pdf", scale))
  
  if (safe_pdf(pdf_file, 14, 9)) {
    print(p_comparison)
    dev.off()
    cat(sprintf("    + Saved: %s\n", basename(pdf_file)))
  }
  
  while (dev.cur() > 1) try(dev.off(), silent = TRUE)
}

#---- Generate summary statistics Excel file ----
cat("\n--- Generating Summary Statistics ---\n")

summary_wb <- createWorkbook()
stats_data <- data.frame(Timescale = character(), Index = character(), 
                         Mean = numeric(), Median = numeric(), 
                         StdDev = numeric(), Min = numeric(), Max = numeric(),
                         Drought_Events_2020_2025 = integer(), stringsAsFactors = FALSE)

for (index_type in c("spi", "spei")) {
  for (scale in timescales) {
    df <- calculate_basin_average(index_type, scale)
    df_recent <- df[df$date >= as.Date("2020-01-01") & df$date <= as.Date("2025-12-31"), ]
    events <- detect_drought_events(df_recent, drought_threshold_onset, drought_threshold_end, 1)
    
    stats_data <- rbind(stats_data, data.frame(
      Timescale = sprintf("%02d", scale),
      Index = toupper(index_type),
      Mean = mean(df$value, na.rm = TRUE),
      Median = median(df$value, na.rm = TRUE),
      StdDev = sd(df$value, na.rm = TRUE),
      Min = min(df$value, na.rm = TRUE),
      Max = max(df$value, na.rm = TRUE),
      Drought_Events_2020_2025 = nrow(events)
    ))
  }
}

addWorksheet(summary_wb, "Summary_Statistics")
writeData(summary_wb, "Summary_Statistics", stats_data, rowNames = FALSE)
addStyle(summary_wb, "Summary_Statistics", 
         style = createStyle(textDecoration = "bold", fgFill = "#D3D3D3"), 
         rows = 1, cols = 1:ncol(stats_data))

#---- Correlation Statistics Sheet ----
cat("\n--- Calculating SPI vs SPEI Correlations ---\n")
correlation_data <- data.frame(
  Timescale = character(),
  Correlation = numeric(),
  SPI_Events = integer(),
  SPEI_Events = integer(),
  stringsAsFactors = FALSE
)

for (scale in timescales) {
  spi_df <- calculate_basin_average("spi", scale)
  spei_df <- calculate_basin_average("spei", scale)
  
  # Merge for correlation calculation
  merged <- merge(spi_df, spei_df, by = "date", all = TRUE)
  corr_val <- cor(merged$value.x, merged$value.y, use = "complete.obs")
  
  # Count drought events for recent period
  spi_recent <- spi_df[spi_df$date >= as.Date("2020-01-01") & spi_df$date <= as.Date("2025-12-31"), ]
  spei_recent <- spei_df[spei_df$date >= as.Date("2020-01-01") & spei_df$date <= as.Date("2025-12-31"), ]
  
  spi_events <- detect_drought_events(spi_recent, drought_threshold_onset, drought_threshold_end, 1)
  spei_events <- detect_drought_events(spei_recent, drought_threshold_onset, drought_threshold_end, 1)
  
  correlation_data <- rbind(correlation_data, data.frame(
    Timescale = sprintf("%02d-month", scale),
    Correlation = round(corr_val, 4),
    SPI_Events = nrow(spi_events),
    SPEI_Events = nrow(spei_events),
    stringsAsFactors = FALSE
  ))
  
  cat(sprintf("  Timescale %02d: r = %.4f | SPI events: %d | SPEI events: %d\n", 
              scale, corr_val, nrow(spi_events), nrow(spei_events)))
}

addWorksheet(summary_wb, "SPI_SPEI_Correlations")
writeData(summary_wb, "SPI_SPEI_Correlations", correlation_data, rowNames = FALSE)
addStyle(summary_wb, "SPI_SPEI_Correlations", 
         style = createStyle(textDecoration = "bold", fgFill = "#D3D3D3"), 
         rows = 1, cols = 1:ncol(correlation_data))

excel_file <- file.path(plot_dir, "SPI_SPEI_Summary_Statistics.xlsx")
saveWorkbook(summary_wb, excel_file, overwrite = TRUE)
cat(sprintf("+ Saved summary statistics: %s\n", basename(excel_file)))

cat("\n*** ALL VISUALIZATIONS COMPLETE ***\n")
cat(sprintf("Output directory: %s\n", normalizePath(plot_dir)))
cat(sprintf("Generated:\n"))
cat(sprintf("  - %d PDF files for individual indices (5 timescales × 2 periods × 2 indices)\n", 
            length(timescales) * 2 * 2))
cat(sprintf("  - %d PDF comparison plots (SPI vs SPEI for each timescale)\n", 
            length(timescales)))
cat(sprintf("  - Total: %d PDF files\n", 
            length(timescales) * 2 * 2 + length(timescales)))
cat(sprintf("  - 1 Excel summary file with correlation statistics\n"))