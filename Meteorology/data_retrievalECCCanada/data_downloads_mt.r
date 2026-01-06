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

# --- 2. Load stations and derive year ranges ---
stations_df <- read_csv(stations_file, show_col_types = FALSE) %>%
  filter(!is.na(STATION_ID)) %>%
  distinct(STATION_ID, STATION_NAME, .keep_all = TRUE) %>%
  mutate(
    # Fallback logic to find the widest possible date range
    start_date_raw = coalesce(DLY_FIRST_DATE, FIRST_DATE),
    end_date_raw   = coalesce(DLY_LAST_DATE, LAST_DATE),
    
    start_year = year(as.Date(substr(as.character(start_date_raw), 1, 10))),
    end_year   = year(as.Date(substr(as.character(end_date_raw), 1, 10)))
  )

# --- 3. Helper Functions ---

# Function to clean and normalize the varying column names from ECCC
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
  
  # Standardized candidates for ECCC headers (which change over time)
  find_col <- function(cands) {
    hits <- cands[cands %in% up]
    if (length(hits) > 0) nm[match(hits[1], up)] else NA_character_
  }
  
  tmax_col   <- find_col(c( "MAX TEMP (°C)", "MAX TEMP (C)", "TX"))
  tmin_col   <- find_col(c( "MIN TEMP (°C)", "MIN TEMP (C)", "TN"))
  tmean_col  <- find_col(c( "MEAN TEMP (°C)", "MEAN TEMP (C)", "TMEAN"))
  precip_col <- find_col(c( "TOTAL PRECIP (MM)", "TOTAL PRECIPITATION (MM)", "PRCP"))
  rain_col   <- find_col(c( "TOTAL RAIN (MM)", "RAIN"))
  snow_col   <- find_col(c( "TOTAL SNOW (CM)", "SNOW"))
  
  sel <- c(date_col, tmax_col, tmin_col, tmean_col, precip_col, rain_col, snow_col)
  sel <- sel[!is.na(sel)]
  
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
      Date = as.Date(substr(as.character(Date), 1, 10)),
      across(where(is.character) & !Date, ~suppressWarnings(as.numeric(.)))
    ) %>%
    filter(!is.na(Date))
  
  return(out)
}

# --- 4. Main Download Logic (The Loop) ---

results <- tibble(STATION_ID = character(), STATION_NAME = character(), rows_saved = integer())

# This function loops through every year to ensure no data is missed
download_full_history <- function(stn_id, stn_name, s_year, e_year) {
  message("\nProcessing: ", stn_id, " - ", stn_name)
  
  if (is.na(s_year) | is.na(e_year)) {
    message("  -> Skip: No date range in metadata.")
    return(NULL)
  }
  
  all_years_list <- list()
  
  for (yr in s_year:e_year) {
    # Requesting 1 year at a time is the only reliable way to get long histories
    url <- paste0("https://climate.weather.gc.ca/climate_data/bulk_data_e.html?format=csv",
                  "&stationID=", stn_id,
                  "&Year=", yr,
                  "&Month=1&Day=1&timeframe=2")
    
    resp <- tryCatch(GET(url, timeout(60)), error = function(e) NULL)
    
    if (!is.null(resp) && status_code(resp) == 200) {
      txt <- content(resp, "text", encoding = "UTF-8")
      if (nchar(txt) > 200) { # Ignore empty CSVs
        raw_df <- read_csv(txt, show_col_types = FALSE, skip = 0)
        clean_df <- normalize_daily(raw_df)
        if (!is.null(clean_df) && nrow(clean_df) > 0) {
          all_years_list[[as.character(yr)]] <- clean_df
        }
      }
    }
  }
  
  if (length(all_years_list) > 0) {
    final_stn_df <- bind_rows(all_years_list) %>% arrange(Date) %>% distinct()
    out_path <- file.path(output_dir, paste0(stn_id, "_CLIMATE_DAILY.csv"))
    write_csv(final_stn_df, out_path)
    message("  -> SUCCESS: Saved ", nrow(final_stn_df), " rows.")
    results <<- bind_rows(results, tibble(STATION_ID=as.character(stn_id), STATION_NAME=stn_name, rows_saved=nrow(final_stn_df)))
  } else {
    message("  -> FAILED: No data found for any year in range.")
  }
}

# --- 5. Execution ---
pwalk(
  stations_df %>% select(STATION_ID, STATION_NAME, start_year, end_year),
  function(STATION_ID, STATION_NAME, start_year, end_year) {
    download_full_history(STATION_ID, STATION_NAME, start_year, end_year)
  }
)

message("\n✅ All downloads complete. Check the '", output_dir, "' folder.")
