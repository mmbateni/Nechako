# Required packages
library(tidyverse)   # includes dplyr, readr, purrr, stringr
library(httr)        # HTTP requests
library(jsonlite)    # JSON parsing
library(purrr)       # map_dfr etc.

# --- 0 Setup output folder --------------------------------------------------
output_dir <- "data_downloads_geomet_api"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
message("WSC GeoMet API Data will be saved in: ", output_dir)
message("-------------------------------------------------------")

# --- 1 Read stations list and filter ACTIVE --------------------------------
stations_file <- if (file.exists("Stations.csv")) "Stations.csv" else "stations.csv"
if (!file.exists(stations_file)) stop("Could not find Stations.csv or stations.csv in the working directory.")
stations_df <- read_csv(stations_file, show_col_types = FALSE)

active_stations_data <- stations_df %>%
  filter(HYD_STATUS == "ACTIVE") %>%
  select(STATION_NUMBER, STATION_NAME)

message("Found ", nrow(active_stations_data), " ACTIVE stations to process.")

# --- 2 Base URL for hydrometric daily mean collection -----------------------
BASE_URL <- "https://api.weather.gc.ca/collections/hydrometric-daily-mean/items"

# --- 3 Helper to GET and parse JSON safely ----------------------------------
safe_get_json <- function(url, timeout_seconds = 60) {
  resp <- tryCatch(GET(url, timeout(timeout_seconds)), error = function(e) e)
  if (inherits(resp, "error")) {
    warning("HTTP request error: ", resp$message)
    return(NULL)
  }
  if (http_error(resp)) {
    warning("HTTP error ", status_code(resp), " for URL: ", url)
    return(NULL)
  }
  raw_text <- content(resp, "text", encoding = "UTF-8")
  parsed <- tryCatch(fromJSON(raw_text, simplifyVector = FALSE), error = function(e) {
    warning("Failed to parse JSON: ", e$message)
    return(NULL)
  })
  list(parsed = parsed, raw_text = raw_text)
}

# --- 4 Fetch function for DISCHARGE or LEVEL --------------------------------
fetch_wsc_parameter <- function(STATION_NUMBER, param_name, limit = 70000) {
  param_name <- toupper(param_name)
  if (!param_name %in% c("DISCHARGE", "LEVEL")) {
    stop("param_name must be 'DISCHARGE' or 'LEVEL'")
  }
  
  # Build URL: query by STATION_NUMBER and limit
  url <- paste0(BASE_URL,
                "?STATION_NUMBER=", URLencode(STATION_NUMBER, reserved = TRUE),
                "&limit=", limit)
  
  res <- safe_get_json(url)
  if (is.null(res) || is.null(res$parsed)) return(NULL)
  data_list <- res$parsed
  
  # Check features
  if (!("features" %in% names(data_list)) || length(data_list$features) == 0) {
    warning(sprintf("No features returned for station %s (param %s). Showing JSON snippet for debugging.", STATION_NUMBER, param_name))
    cat(substr(res$raw_text, 1, 2000), "\n")
    return(NULL)
  }
  
  # Convert features -> tibble of properties, replacing NULL with NA
  props_df <- map_dfr(data_list$features, function(f) {
    props <- f$properties
    # Replace NULLs with NA so as_tibble works
    props_clean <- lapply(props, function(x) if (is.null(x)) NA else x)
    # Convert to tibble row
    as_tibble(props_clean)
  })
  
  # Normalize column names to uppercase for robust matching
  names(props_df) <- toupper(names(props_df))
  
  # Ensure DATE exists (GeoMet uses DATE)
  if (!"DATE" %in% names(props_df)) {
    if ("DATETIME" %in% names(props_df)) {
      props_df <- props_df %>% rename(DATE = DATETIME)
    } else {
      warning("No DATE or DATETIME field found in properties for station ", STATION_NUMBER)
      return(NULL)
    }
  }
  
  # Ensure requested parameter exists
  if (!param_name %in% names(props_df)) {
    warning(sprintf("Requested parameter '%s' not present for station %s. Available properties: %s",
                    param_name, STATION_NUMBER, paste(names(props_df), collapse = ", ")))
    return(NULL)
  }
  
  # Select and clean
  df <- props_df %>%
    select(Date = DATE, Value = all_of(param_name)) %>%
    mutate(Date = as.Date(str_sub(as.character(Date), 1, 10))) %>%
    arrange(Date) %>%
    filter(!is.na(Value))
  
  if (nrow(df) == 0) return(NULL)
  return(df)
}

# --- 5 Main download function -----------------------------------------------
download_api_data <- function(STATION_NUMBER, STATION_NAME) {
  message("Processing station: ", STATION_NUMBER, " - ", STATION_NAME)
  
  # DISCHARGE
  flow_df <- fetch_wsc_parameter(STATION_NUMBER, "DISCHARGE")
  if (!is.null(flow_df) && nrow(flow_df) > 0) {
    flow_df <- flow_df %>% rename(Discharge_cms = Value)
    flow_file <- file.path(output_dir, paste0(STATION_NUMBER, "_WSC_DISCHARGE.csv"))
    write_csv(flow_df, flow_file)
    message("  -> Saved ", nrow(flow_df), " rows of DISCHARGE data to: ", flow_file)
  } else {
    message("  -> WARNING: No DISCHARGE data found for: ", STATION_NUMBER)
  }
  
  # LEVEL
  level_df <- fetch_wsc_parameter(STATION_NUMBER, "LEVEL")
  if (!is.null(level_df) && nrow(level_df) > 0) {
    level_df <- level_df %>% rename(Level_m = Value)
    level_file <- file.path(output_dir, paste0(STATION_NUMBER, "_WSC_LEVEL.csv"))
    write_csv(level_df, level_file)
    message("  -> Saved ", nrow(level_df), " rows of LEVEL data to: ", level_file)
  } else {
    message("  -> WARNING: No LEVEL data found for: ", STATION_NUMBER)
  }
}

# --- 6 Execute downloads for all active stations ----------------------------
pwalk(active_stations_data, download_api_data)

message("-------------------------------------------------------")
message("✅ Historical data download via WSC GeoMet API complete.")

# --- 7 Build basin-wide daily Excel matrices ---------------------------------
# These exports use the CSV files written above. Each row is one calendar day;
# each gauge is a separate column. Missing observations remain blank in Excel.
# Requires openxlsx: install.packages("openxlsx") if not already installed.
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop("Package 'openxlsx' is required for Excel export. Install it with install.packages('openxlsx').")
}

build_daily_wide_matrix <- function(parameter = c("DISCHARGE", "LEVEL")) {
  parameter <- match.arg(parameter)
  suffix <- if (parameter == "DISCHARGE") "DISCHARGE" else "LEVEL"
  value_column <- if (parameter == "DISCHARGE") "Discharge_cms" else "Level_m"
  files <- list.files(
    output_dir,
    pattern = paste0("_WSC_", suffix, "\\.csv$"),
    full.names = TRUE
  )

  if (length(files) == 0L) {
    warning("No downloaded ", parameter, " CSV files found in: ", output_dir)
    return(NULL)
  }

  station_lookup <- stations_df %>%
    dplyr::distinct(STATION_NUMBER, STATION_NAME) %>%
    dplyr::mutate(StationColumn = paste0(STATION_NUMBER, " - ", STATION_NAME))

  long_data <- purrr::map_dfr(files, function(file) {
    station_no <- sub(paste0("_WSC_", suffix, "\\.csv$"), "", basename(file))
    dat <- readr::read_csv(file, show_col_types = FALSE, col_types = readr::cols())

    if (!all(c("Date") %in% names(dat))) {
      warning("Skipping file without a Date column: ", file)
      return(tibble::tibble())
    }

    actual_value_column <- if (value_column %in% names(dat)) value_column else
      if ("Value" %in% names(dat)) "Value" else NA_character_
    if (is.na(actual_value_column)) {
      warning("Skipping file without a recognised value column: ", file)
      return(tibble::tibble())
    }

    station_name <- station_lookup$STATION_NAME[match(station_no, station_lookup$STATION_NUMBER)]
    if (length(station_name) == 0L || is.na(station_name)) station_name <- station_no

    dat %>%
      dplyr::transmute(
        Date = as.Date(Date),
        StationColumn = paste0(station_no, " - ", station_name),
        Value = suppressWarnings(as.numeric(.data[[actual_value_column]]))
      ) %>%
      dplyr::filter(!is.na(Date)) %>%
      dplyr::group_by(Date, StationColumn) %>%
      dplyr::summarise(Value = dplyr::last(Value[!is.na(Value)], default = NA_real_),
                       .groups = "drop")
  })

  if (nrow(long_data) == 0L) {
    warning("No valid daily ", parameter, " observations available for Excel export.")
    return(NULL)
  }

  # Continuous daily calendar spanning the actual earliest and latest observations.
  date_grid <- tibble::tibble(Date = seq(min(long_data$Date), max(long_data$Date), by = "day"))
  station_columns <- station_lookup$StationColumn[
    station_lookup$STATION_NUMBER %in% sub(" - .*", "", unique(long_data$StationColumn))
  ]
  station_columns <- unique(c(station_columns, unique(long_data$StationColumn)))

  wide <- date_grid %>%
    dplyr::left_join(
      long_data %>% tidyr::pivot_wider(names_from = StationColumn, values_from = Value),
      by = "Date"
    )

  missing_columns <- setdiff(station_columns, names(wide))
  for (nm in missing_columns) wide[[nm]] <- NA_real_
  wide <- wide %>% dplyr::select(Date, dplyr::all_of(station_columns)) %>% dplyr::arrange(Date)

  output_file <- file.path(
    output_dir,
    if (parameter == "DISCHARGE") "Nechako_Daily_Discharge.xlsx" else "Nechako_Daily_Heights.xlsx"
  )

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Daily data")
  openxlsx::writeData(wb, "Daily data", wide, withFilter = TRUE)
  openxlsx::freezePane(wb, "Daily data", firstRow = TRUE, firstCol = TRUE)
  openxlsx::addStyle(
    wb, "Daily data",
    style = openxlsx::createStyle(numFmt = "yyyy-mm-dd"),
    rows = 2:(nrow(wide) + 1L), cols = 1L, gridExpand = TRUE, stack = TRUE
  )
  openxlsx::setColWidths(wb, "Daily data", cols = 1L, widths = 13)
  if (ncol(wide) > 1L) openxlsx::setColWidths(wb, "Daily data", cols = 2:ncol(wide), widths = 22)
  openxlsx::saveWorkbook(wb, output_file, overwrite = TRUE)

  message(sprintf(
    "Excel matrix saved: %s | %s to %s | %s daily rows | %s gauges",
    output_file, min(wide$Date), max(wide$Date), nrow(wide), ncol(wide) - 1L
  ))
  invisible(wide)
}

message("\nCreating daily wide Excel workbooks...")
discharge_wide <- build_daily_wide_matrix("DISCHARGE")
heights_wide <- build_daily_wide_matrix("LEVEL")
message("Daily wide Excel export complete.")

# Complexity: for N downloaded observations and G gauges, reshaping is O(N)
# expected time plus O(D * G) output memory, where D is the number of calendar
# days from the global earliest to latest observation. Excel's worksheet limit
# is 1,048,576 rows and 16,384 columns; this daily layout is below those limits
# for ordinary hydrometric record spans, but exceptionally long/wide datasets
# should be checked before writing.
