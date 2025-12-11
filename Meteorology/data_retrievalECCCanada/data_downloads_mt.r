# --- 0. Libraries ---
library(tidyverse)
library(httr)
library(readr)
library(lubridate)
library(purrr)

# --- 1. Config ---
stations_file <- "Stations.csv"
output_dir    <- "climate_daily_downloads"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Optional override for date range (set both to NULL to auto-derive from station rows)
override_start_year <- NULL  # e.g., 1960
override_end_year   <- NULL  # e.g., 2025

# ECCC bulk daily CSV endpoint
# timeframe=2 => daily, format=csv
BULK_BASE <- "https://climate.weather.gc.ca/climate_data/bulk_data_e.html"

message("Starting ECCC bulk daily download (temperature & precipitation) via Stations.csv")
message("--------------------------------------------------------------------------")

# --- 2. Load stations and derive years ---
stations_df <- read_csv(stations_file, show_col_types = FALSE) %>%
  filter(!is.na(STATION_ID)) %>%
  distinct(STATION_ID, STATION_NAME, .keep_all = TRUE) %>%
  mutate(
    # Prefer daily first/last; fallback to station first/last if daily dates are missing
    start_date_raw = coalesce(DLY_FIRST_DATE, FIRST_DATE),
    end_date_raw   = coalesce(DLY_LAST_DATE, LAST_DATE),
    
    # Extract year (YYYY) safely
    start_year = suppressWarnings(year(as.Date(substr(as.character(start_date_raw), 1, 10)))),
    end_year   = suppressWarnings(year(as.Date(substr(as.character(end_date_raw), 1, 10)))),
    
    # Apply overrides if provided
    start_year = if (!is.null(override_start_year)) override_start_year else start_year,
    end_year   = if (!is.null(override_end_year))   override_end_year   else end_year
  )

message("Found ", nrow(stations_df), " stations to process")

# --- 3. Helpers ---

build_bulk_url <- function(station_id, start_year, end_year) {
  # Fallback years if missing
  if (is.na(start_year)) start_year <- 1840  # earliest possible
  if (is.na(end_year))   end_year   <- as.integer(format(Sys.Date(), "%Y"))
  
  paste0(
    BULK_BASE,
    "?format=csv",
    "&stationID=", URLencode(as.character(station_id), reserved = TRUE),
    "&timeframe=2",
    "&startYear=", start_year,
    "&endYear=",   end_year
  )
}

safe_download_csv <- function(url) {
  resp <- tryCatch(GET(url, timeout(120)), error = function(e) e)
  if (inherits(resp, "error")) {
    warning("HTTP request error: ", resp$message)
    return(NULL)
  }
  if (http_error(resp)) {
    warning("HTTP error ", status_code(resp), " for URL: ", url)
    return(NULL)
  }
  txt <- content(resp, "text", encoding = "UTF-8")
  # Some responses may be small if no data; guard against empty
  if (nchar(txt) < 50) return(NULL)
  df <- tryCatch(read_csv(txt, show_col_types = FALSE), error = function(e) NULL)
  df
}

normalize_daily <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  nm <- names(df)
  up <- toupper(nm)
  
  # Identify date column
  date_col <- case_when(
    "DATE/TIME" %in% up ~ nm[match("DATE/TIME", up)],
    "DATE"      %in% up ~ nm[match("DATE", up)],
    TRUE ~ NA_character_
  )
  
  # Candidates for temperature and precip fields (bulk CSV schema varies by era)
  tmax_candidates   <- c("MAX TEMP (°C)", "MAX TEMP (C)", "MAXIMUM TEMPERATURE (°C)", "TX")
  tmin_candidates   <- c("MIN TEMP (°C)", "MIN TEMP (C)", "MINIMUM TEMPERATURE (°C)", "TN")
  tmean_candidates  <- c("MEAN TEMP (°C)", "MEAN TEMP (C)", "MEAN TEMPERATURE (°C)", "TMEAN", "TAVG")
  precip_total_cand <- c("TOTAL PRECIP (MM)", "TOTAL PRECIPITATION (MM)", "TOTAL PRECIPITATION", "TOTAL PRECIP", "PRCP")
  rain_total_cand   <- c("TOTAL RAIN (MM)", "TOTAL RAINFALL (MM)", "TOTAL RAINFALL", "RAIN")
  snow_total_cand   <- c("TOTAL SNOW (CM)", "TOTAL SNOWFALL (CM)", "TOTAL SNOWFALL", "SNOW")
  
  find_col <- function(cands) {
    hits <- cands[cands %in% up]
    if (length(hits) > 0) nm[match(hits[1], up)] else NA_character_
  }
  
  tmax_col   <- find_col(tmax_candidates)
  tmin_col   <- find_col(tmin_candidates)
  tmean_col  <- find_col(tmean_candidates)
  precip_col <- find_col(precip_total_cand)
  rain_col   <- find_col(rain_total_cand)
  snow_col   <- find_col(snow_total_cand)
  
  sel <- c(date_col, tmax_col, tmin_col, tmean_col, precip_col, rain_col, snow_col)
  sel <- sel[!is.na(sel)]
  if (length(sel) == 0) return(NULL)
  
  out <- df %>%
    select(all_of(sel)) %>%
    rename(
      Date            = any_of(date_col),
      Tmax_C          = any_of(tmax_col),
      Tmin_C          = any_of(tmin_col),
      Tmean_C         = any_of(tmean_col),
      Total_Precip_mm = any_of(precip_col),
      Rain_mm         = any_of(rain_col),
      Snow_cm         = any_of(snow_col)
    ) %>%
    mutate(
      Date = suppressWarnings(as.Date(substr(as.character(Date), 1, 10))),
      # Coerce numeric columns safely
      across(c(Tmax_C, Tmin_C, Tmean_C, Total_Precip_mm, Rain_mm, Snow_cm),
             ~ suppressWarnings(as.numeric(.)))
    ) %>%
    arrange(Date)
  
  # Keep rows that have a date and at least one climate value
  val_cols <- intersect(c("Tmax_C", "Tmin_C", "Tmean_C", "Total_Precip_mm", "Rain_mm", "Snow_cm"), names(out))
  out <- out %>% filter(!is.na(Date), if_any(all_of(val_cols), ~ !is.na(.)))
  if (nrow(out) == 0) return(NULL)
  out
}

# --- 4. Loop with results tracking ---
results <- tibble(STATION_ID = character(), STATION_NAME = character(), rows_saved = integer())

download_station_daily <- function(station_id, station_name, start_year, end_year) {
  message("Processing: ", station_id, " - ", station_name,
          " (years: ", coalesce(start_year, NA_integer_), "–", coalesce(end_year, NA_integer_), ")")
  url <- build_bulk_url(station_id, start_year, end_year)
  raw_df <- safe_download_csv(url)
  
  if (is.null(raw_df) || nrow(raw_df) == 0) {
    message("  -> WARNING: No CSV rows returned.")
    results <<- bind_rows(results, tibble(STATION_ID = as.character(station_id),
                                          STATION_NAME = station_name,
                                          rows_saved = 0L))
    return(invisible(NULL))
  }
  
  daily <- normalize_daily(raw_df)
  if (is.null(daily) || nrow(daily) == 0) {
    message("  -> WARNING: CSV returned, but no recognized daily temp/precip fields.")
    results <<- bind_rows(results, tibble(STATION_ID = as.character(station_id),
                                          STATION_NAME = station_name,
                                          rows_saved = 0L))
    return(invisible(NULL))
  }
  
  out_file <- file.path(output_dir, paste0(station_id, "_CLIMATE_DAILY.csv"))
  write_csv(daily, out_file)
  message("  -> Saved ", nrow(daily), " rows to: ", out_file)
  
  results <<- bind_rows(results, tibble(STATION_ID = as.character(station_id),
                                        STATION_NAME = station_name,
                                        rows_saved = nrow(daily)))
  invisible(daily)
}

pwalk(
  stations_df %>% select(STATION_ID, STATION_NAME, start_year, end_year),
  function(STATION_ID, STATION_NAME, start_year, end_year) {
    download_station_daily(STATION_ID, STATION_NAME, start_year, end_year)
  }
)

# --- 5. Summary ---
# Ensure all stations represented
results <- results %>%
  right_join(stations_df %>% select(STATION_ID, STATION_NAME) %>%
               mutate(STATION_ID = as.character(STATION_ID)),
             by = c("STATION_ID", "STATION_NAME")) %>%
  mutate(rows_saved = replace_na(rows_saved, 0L))

total     <- nrow(stations_df)
with_data <- sum(results$rows_saved > 0)
no_data   <- sum(results$rows_saved == 0)

message("--------------------------------------------------------------------------")
message("Summary:")
message("  Total stations listed:   ", total)
message("  Stations with data:      ", with_data)
message("  Stations with no data:   ", no_data)

if (no_data > 0) {
  nd <- results %>% filter(rows_saved == 0) %>% arrange(STATION_NAME)
  message("Stations with no daily temp/precip data:")
  print(nd %>% select(STATION_ID, STATION_NAME))
}

if (with_data > 0) {
  wd <- results %>% filter(rows_saved > 0) %>% arrange(desc(rows_saved))
  message("Top stations by rows saved:")
  print(wd %>% select(STATION_ID, STATION_NAME, rows_saved) %>% head(10))
}

message("✅ Finished ECCC bulk daily downloads")
