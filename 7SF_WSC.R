# ============================================================================
# NECHAKO StreamFlow RIVER BASIN DROUGHT ANALYSIS 
# Data Loading, Preprocessing & Drought Identification
# Based on Raut & Ganguli (2024) methodology
# All folders inside sf_results
# ============================================================================

# Clear workspace
rm(list = ls())

# ============================================================================
# SECTION 1: PACKAGE INSTALLATION & LOADING
# ============================================================================

# Install required packages (run once)
packages_needed <- c(
  "tidyverse",
  "lubridate", 
  "zoo",
  "circular",
  "MASS",
  "trend",
  "Kendall",
  "boot",
  "copula",
  "kde1d"
)

for(pkg in packages_needed) {
  if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.rstudio.com/")
    library(pkg, character.only = TRUE)
  }
}

# ============================================================================
# SECTION 2: STATION METADATA & DIRECTORY SETUP
# ============================================================================

# Define all 7 Nechako River Basin stations
stations <- data.frame(
  StationID = c("08JA015", "08JB002", "08JC001", "08JC002", 
                "08JC005", "08JE001", "08JE004"),
  StartYear = c(1976, 1929, 1915, 1950, 1953, 1929, 1975),
  EndYear = c(2022, 2024, 2023, 2023, 2022, 2024, 2024),
  stringsAsFactors = FALSE
)

# ============================================================================
# DIRECTORY SETUP - ALL FOLDERS INSIDE sf_results
# ============================================================================

# Set working directory
setwd("D:/Nechako_Drought/Nechako")

# Input data directory 
INPUT_DATA_DIR <- "D:/Nechako_Drought/Nechako/Hydrology/data_retrievalWaterSurveyofCanada/data_downloads_geomet_api"

# Create MAIN output folder
MAIN_OUTPUT_DIR <- "sf_results"
dir.create(MAIN_OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Create ALL subfolders INSIDE sf_results
dir.create(file.path(MAIN_OUTPUT_DIR, "data"), showWarnings = FALSE)
dir.create(file.path(MAIN_OUTPUT_DIR, "processed_data"), showWarnings = FALSE)
dir.create(file.path(MAIN_OUTPUT_DIR, "drought_events"), showWarnings = FALSE)
dir.create(file.path(MAIN_OUTPUT_DIR, "thresholds"), showWarnings = FALSE)
dir.create(file.path(MAIN_OUTPUT_DIR, "figures"), showWarnings = FALSE)
dir.create(file.path(MAIN_OUTPUT_DIR, "reports"), showWarnings = FALSE)

# List available CSV files in the directory
available_files <- list.files(INPUT_DATA_DIR, pattern = "\\.csv$", full.names = TRUE)
cat(sprintf("Available CSV files in directory: %d\n", length(available_files)))
if(length(available_files) > 0) {
  cat("Sample files:\n")
  print(head(basename(available_files)))
}
cat("\n")

# ============================================================================
# SECTION 3: DATA LOADING FUNCTION
# ============================================================================

load_discharge_data <- function(station_id, data_path = INPUT_DATA_DIR) {
  
  # Try different possible file names (adjusted for your file naming convention)
  file_patterns <- c(
    file.path(data_path, paste0(station_id, "_WSC_DISCHARGE.csv")),
    file.path(data_path, paste0(station_id, "_discharge.csv")),
    file.path(data_path, paste0(station_id, ".csv")),
    file.path(data_path, paste0("WSC_", station_id, ".csv")),
    file.path(data_path, paste0(station_id, "_WaterSurveyofCanada.csv"))
  )
  
  file_found <- FALSE
  data <- NULL
  used_file <- NA
  
  for(file_name in file_patterns) {
    if(file.exists(file_name)) {
      cat(sprintf("  ✓ Loading: %s\n", basename(file_name)))
      data <- read.csv(file_name, stringsAsFactors = FALSE)
      file_found <- TRUE
      used_file <- file_name
      break
    }
  }
  
  if(!file_found) {
    warning(sprintf("  ⚠ No data file found for station %s in %s", station_id, data_path))
    cat(sprintf("  Searched patterns: %s\n", paste(basename(file_patterns), collapse = ", ")))
    return(NULL)
  }
  
  # Standardize column names
  colnames(data) <- tolower(colnames(data))
  
  # Find date column
  date_col <- grep("date", colnames(data), value = TRUE, ignore.case = TRUE)
  if(length(date_col) == 0) {
    stop(sprintf("No date column found for station %s", station_id))
  }
  data$date <- as.Date(data[[date_col[1]]])
  
  # Find discharge column
  discharge_col <- grep("discharge|flow|q_", colnames(data), 
                        value = TRUE, ignore.case = TRUE)
  if(length(discharge_col) == 0) {
    stop(sprintf("No discharge column found for station %s", station_id))
  }
  data$discharge <- as.numeric(data[[discharge_col[1]]])
  
  # Remove NA values
  data <- data[!is.na(data$date) & !is.na(data$discharge), ]
  
  # Add time columns
  data$year <- year(data$date)
  data$month <- month(data$date)
  data$day <- day(data$date)
  data$doy <- yday(data$date)  # Day of year (1-366)
  
  # Add station ID
  data$station_id <- station_id
  
  cat(sprintf("  Records loaded: %d (%d-%d)\n", 
              nrow(data), min(data$year), max(data$year)))
  
  return(data)
}

# ============================================================================
# SECTION 4: DATA COMPLETENESS CHECK
# ============================================================================

check_data_completeness <- function(data, station_id) {
  
  if(is.null(data) || nrow(data) == 0) {
    return(list(
      station_id = station_id,
      start_year = NA,
      end_year = NA,
      total_days = 0,
      expected_days = 0,
      completeness = 0,
      meets_threshold = FALSE,
      message = "No data available"
    ))
  }
  
  # Get year range
  year_range <- range(data$year)
  total_years <- diff(year_range) + 1
  
  # Calculate expected days (accounting for leap years)
  expected_days <- sum(sapply(year_range[1]:year_range[2], function(y) {
    if(leap_year(y)) 366 else 365
  }))
  
  actual_days <- nrow(data)
  completeness <- (actual_days / expected_days) * 100
  
  cat(sprintf("\nStation: %s\n", station_id))
  cat(sprintf("  Period: %d-%d (%d years)\n", year_range[1], year_range[2], total_years))
  cat(sprintf("  Data: %d / %d days (%.2f%% complete)\n", 
              actual_days, expected_days, completeness))
  
  # Check if meets 70% threshold (per Raut & Ganguli 2024)
  meets_threshold <- completeness >= 70
  
  if(!meets_threshold) {
    warning(sprintf("  ⚠ Station %s has <70%% completeness (%.2f%%)", 
                    station_id, completeness))
  } else {
    cat(sprintf("  ✓ Meets 70%% completeness threshold\n"))
  }
  
  return(list(
    station_id = station_id,
    start_year = year_range[1],
    end_year = year_range[2],
    total_years = total_years,
    total_days = actual_days,
    expected_days = expected_days,
    completeness = completeness,
    meets_threshold = meets_threshold
  ))
}

# ============================================================================
# SECTION 5: MISSING DATA FILLING
# ============================================================================

fill_missing_data <- function(data, method = "zoo") {
  
  if(is.null(data) || nrow(data) == 0) {
    return(NULL)
  }
  
  # Create complete date sequence
  full_dates <- seq(min(data$date), max(data$date), by = "day")
  
  # Merge with existing data
  complete_data <- data.frame(date = full_dates)
  complete_data <- merge(complete_data, data, by = "date", all.x = TRUE)
  complete_data <- complete_data[order(complete_data$date), ]
  
  # Fill missing discharge values
  if(method == "zoo") {
    # Linear interpolation using zoo
    complete_data$discharge <- zoo::na.approx(
      complete_data$discharge, 
      na.rm = FALSE, 
      rule = 2  # Carry last observation forward for edges
    )
  } else if(method == "locf") {
    # Last observation carried forward
    complete_data$discharge <- zoo::na.locf(
      complete_data$discharge,
      na.rm = FALSE,
      rule = 2
    )
  }
  
  # Re-add time columns
  complete_data$year <- year(complete_data$date)
  complete_data$month <- month(complete_data$date)
  complete_data$day <- day(complete_data$date)
  complete_data$doy <- yday(complete_data$date)
  
  # Add station ID if not present
  if("station_id" %in% colnames(data)) {
    complete_data$station_id <- data$station_id[1]
  }
  
  # Count filled values
  original_na <- sum(is.na(data$discharge))
  filled_na <- sum(is.na(complete_data$discharge))
  
  cat(sprintf("  Missing data filled: %d values\n", original_na - filled_na))
  
  return(complete_data)
}

# ============================================================================
# SECTION 6: VARIABLE THRESHOLD CALCULATION
# ============================================================================

calculate_variable_threshold <- function(data, window_size = 31, 
                                         exceedance_prob = 0.80,
                                         station_id = NULL) {
  
  if(is.null(data) || nrow(data) == 0) {
    return(NULL)
  }
  
  # Initialize threshold dataframe
  thresholds <- data.frame(
    doy = 1:366,
    threshold = NA,
    threshold_smooth = NA
  )
  
  # Calculate threshold for each day of year
  for(d in 1:366) {
    # Get all discharge values for this day of year across all years
    doy_values <- data$discharge[data$doy == d]
    
    if(length(doy_values) >= 10) {  # Minimum 10 years of data
      # Calculate 80th percentile (20% exceedance = 80th percentile)
      # Per Raut & Ganguli 2024: 80% exceedance probability threshold
      thresholds$threshold[d] <- quantile(doy_values, probs = exceedance_prob, 
                                          na.rm = TRUE)
    } else {
      # Use neighboring days if insufficient data
      neighbors <- c()
      for(offset in 1:15) {
        d_prev <- ((d - offset - 1) %% 366) + 1
        d_next <- ((d + offset - 1) %% 366) + 1
        neighbors <- c(neighbors, 
                       data$discharge[data$doy == d_prev],
                       data$discharge[data$doy == d_next])
        if(length(neighbors) >= 10) break
      }
      thresholds$threshold[d] <- quantile(neighbors, probs = exceedance_prob, 
                                          na.rm = TRUE)
    }
  }
  
  # Apply centered moving average smoothing (31 days as per paper)
  thresholds$threshold_smooth <- zoo::rollapply(
    thresholds$threshold,
    width = window_size,
    FUN = mean,
    fill = NA,
    align = "center",
    partial = TRUE
  )
  
  # Handle leap year (day 366)
  thresholds$threshold_smooth[366] <- thresholds$threshold_smooth[365]
  
  # Add metadata
  thresholds$station_id <- station_id
  thresholds$window_size <- window_size
  thresholds$exceedance_prob <- exceedance_prob
  
  cat(sprintf("  Threshold calculated for DOY 1-366 (smoothed: %d-day window)\n", 
              window_size))
  
  return(thresholds)
}

# ============================================================================
# SECTION 7: DROUGHT EVENT IDENTIFICATION
# ============================================================================

identify_droughts <- function(data, thresholds, min_duration = 30, 
                              station_id = NULL) {
  
  if(is.null(data) || nrow(data) == 0 || is.null(thresholds)) {
    return(list(
      daily_data = data,
      drought_events = NULL,
      message = "Insufficient data for drought identification"
    ))
  }
  
  # Merge thresholds with daily data
  data$doy <- ifelse(leap_year(data$year) & data$doy > 366, 366, data$doy)
  data <- merge(data, thresholds[, c("doy", "threshold_smooth")], 
                by = "doy", all.x = TRUE)
  data <- data[order(data$date), ]
  
  # Identify days below threshold
  data$below_threshold <- data$discharge < data$threshold_smooth
  
  # Identify drought events (consecutive days below threshold)
  data$event_id <- with(rle(data$below_threshold), 
                        rep(seq_along(lengths), lengths))
  
  # Extract drought event information
  event_info <- data %>%
    group_by(event_id) %>%
    filter(below_threshold == TRUE) %>%
    summarise(
      start_date = min(date),
      end_date = max(date),
      duration = n(),
      start_doy = first(doy),
      end_doy = last(doy),
      discharge_min = min(discharge, na.rm = TRUE),
      discharge_mean = mean(discharge, na.rm = TRUE),
      threshold_mean = mean(threshold_smooth, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(duration >= min_duration) %>%  # Minimum 30 days as per paper
    mutate(
      event_id_new = row_number(),
      deficit_volume = sum(threshold_smooth - discharge),  # Total deficit
      deficit_mean = mean(threshold_smooth - discharge),   # Mean daily deficit
      year = year(start_date),
      month = month(start_date),
      season = case_when(
        month %in% c(12, 1, 2) ~ "Winter",
        month %in% c(3, 4, 5) ~ "Spring",
        month %in% c(6, 7, 8) ~ "Summer",
        month %in% c(9, 10, 11) ~ "Fall"
      )
    )
  
  # Add station ID
  if(!is.null(station_id)) {
    event_info$station_id <- station_id
  }
  
  # Merge event info back to daily data
  data <- data %>%
    left_join(event_info %>% select(event_id, event_id_new), by = "event_id")
  
  cat(sprintf("  Drought events identified: %d (min duration: %d days)\n", 
              nrow(event_info), min_duration))
  
  return(list(
    daily_data = data,
    drought_events = event_info,
    n_events = nrow(event_info),
    total_drought_days = sum(event_info$duration),
    message = "Success"
  ))
}

# ============================================================================
# SECTION 8: PROCESS ALL STATIONS
# ============================================================================

# Initialize storage lists
all_processed_data <- list()
all_completeness <- list()
all_thresholds <- list()
all_droughts <- list()

# Process each station
for(i in 1:nrow(stations)) {
  
  station_id <- stations$StationID[i]
  start_year <- stations$StartYear[i]
  end_year <- stations$EndYear[i]
  
  cat(sprintf("\n%s\n", paste(rep("=", 60), collapse = "")))
  cat(sprintf("PROCESSING STATION %d/%d: %s\n", i, nrow(stations), station_id))
  cat(sprintf("%s\n", paste(rep("=", 60), collapse = "")))
  
  # Step 1: Load data
  cat("\n[1/4] Loading discharge data...\n")
  data <- load_discharge_data(station_id)
  
  if(is.null(data)) {
    cat(sprintf("  ⚠ Skipping station %s (no data)\n", station_id))
    next
  }
  
  # Step 2: Check completeness
  cat("\n[2/4] Checking data completeness...\n")
  completeness <- check_data_completeness(data, station_id)
  all_completeness[[station_id]] <- completeness
  
  if(!completeness$meets_threshold) {
    cat(sprintf("  ⚠ Station %s does not meet 70%% threshold but will continue\n", 
                station_id))
  }
  
  # Step 3: Fill missing data
  cat("\n[3/4] Filling missing data...\n")
  data_filled <- fill_missing_data(data)
  all_processed_data[[station_id]] <- data_filled
  
  # Step 4: Calculate thresholds and identify droughts
  cat("\n[4/4] Calculating variable thresholds & identifying droughts...\n")
  
  # Calculate threshold for entire period
  thresholds <- calculate_variable_threshold(
    data_filled, 
    station_id = station_id
  )
  all_thresholds[[station_id]] <- thresholds
  
  # Identify droughts for entire period
  droughts <- identify_droughts(
    data_filled, 
    thresholds,
    station_id = station_id
  )
  all_droughts[[station_id]] <- droughts
  
  cat(sprintf("\n✓ Station %s processing complete\n", station_id))
}

# ============================================================================
# SECTION 9: SAVE ALL RESULTS
# ============================================================================

cat(sprintf("\n%s\n", paste(rep("=", 60), collapse = "")))
cat("SAVING ALL RESULTS...\n")
cat(sprintf("%s\n", paste(rep("=", 60), collapse = "")))

# Save station metadata
saveRDS(stations, file.path(MAIN_OUTPUT_DIR, "stations_metadata.rds"))
cat("✓ Saved: stations_metadata.rds\n")

# Save completeness reports
saveRDS(all_completeness, file.path(MAIN_OUTPUT_DIR, "completeness_reports.rds"))
cat("✓ Saved: completeness_reports.rds\n")

# Save processed data for each station
for(station_id in names(all_processed_data)) {
  saveRDS(all_processed_data[[station_id]], 
          file.path(MAIN_OUTPUT_DIR, "processed_data", 
                    paste0(station_id, "_processed_data.rds")))
}
cat("✓ Saved: processed_data/*.rds\n")

# Save thresholds
saveRDS(all_thresholds, file.path(MAIN_OUTPUT_DIR, "thresholds", "thresholds_all.rds"))
cat("✓ Saved: thresholds/thresholds_all.rds\n")

# Save drought events
saveRDS(all_droughts, file.path(MAIN_OUTPUT_DIR, "drought_events", "droughts_all.rds"))
cat("✓ Saved: drought_events/droughts_all.rds\n")

# Save all results as a single master file
all_results <- list(
  stations = stations,
  completeness = all_completeness,
  processed_data = all_processed_data,
  thresholds = all_thresholds,
  droughts = all_droughts,
  metadata = list(
    input_data_dir = INPUT_DATA_DIR,
    output_dir = MAIN_OUTPUT_DIR,
    analysis_type = "SINGLE_PERIOD",
    min_drought_duration = 30,
    threshold_exceedance_prob = 0.80,
    threshold_window_size = 31,
    processed_date = Sys.Date(),
    script_version = "2.0"
  )
)

saveRDS(all_results, file.path(MAIN_OUTPUT_DIR, "all_results_single_period.rds"))
cat("✓ Saved: all_results_single_period.rds (MASTER FILE)\n")

# ============================================================================
# SECTION 10: SUMMARY REPORT
# ============================================================================

cat(sprintf("\n%s\n", paste(rep("=", 60), collapse = "")))
cat("PROCESSING SUMMARY\n")
cat(sprintf("%s\n", paste(rep("=", 60), collapse = "")))

# Create summary table
summary_table <- data.frame(
  Station = character(),
  Period = character(),
  Years = character(),
  Completeness = numeric(),
  Drought_Events = integer(),
  Total_Drought_Days = integer(),
  Avg_Duration = numeric(),
  stringsAsFactors = FALSE
)

for(station_id in names(all_completeness)) {
  comp <- all_completeness[[station_id]]
  
  if(!is.null(all_droughts[[station_id]]$drought_events)) {
    events <- all_droughts[[station_id]]$drought_events
    summary_table <- rbind(summary_table, data.frame(
      Station = station_id,
      Period = "ENTIRE PERIOD",
      Years = sprintf("%d-%d", comp$start_year, comp$end_year),
      Completeness = comp$completeness,
      Drought_Events = nrow(events),
      Total_Drought_Days = sum(events$duration),
      Avg_Duration = ifelse(nrow(events) > 0, 
                            mean(events$duration), NA),
      stringsAsFactors = FALSE
    ))
  }
}

# Save summary
write.csv(summary_table, file.path(MAIN_OUTPUT_DIR, "reports", "station_summary.csv"), 
          row.names = FALSE)
cat("✓ Saved: reports/station_summary.csv\n")

# Print summary
cat("\n")
print(summary_table, row.names = FALSE)

cat(sprintf("\n%s\n", paste(rep("=", 60), collapse = "")))
cat("ANALYSIS COMPLETE - SINGLE PERIOD\n")
cat(sprintf("%s\n", paste(rep("=", 60), collapse = "")))

cat(sprintf("\nAll results saved to: %s/\n", MAIN_OUTPUT_DIR))
cat("Master file: sf_results/all_results_single_period.rds\n")

# Display folder structure
cat("\n=== OUTPUT FOLDER STRUCTURE ===\n")
cat("sf_results/\n")
cat("├── data/\n")
cat("├── processed_data/\n")
cat("├── drought_events/\n")
cat("├── thresholds/\n")
cat("├── figures/\n")
cat("├── reports/\n")
cat("├── stations_metadata.rds\n")
cat("├── completeness_reports.rds\n")
cat("└── all_results_single_period.rds (MASTER FILE)\n")