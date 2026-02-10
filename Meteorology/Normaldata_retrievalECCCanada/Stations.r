library(tidyverse)
library(httr)
library(jsonlite)

# --- 1. Define Area and Output ---
MIN_LAT <- 53.0
MAX_LAT <- 55.0
MIN_LON <- -127.0
MAX_LON <- -122.0

output_file <- "Stations.csv"

message("Starting ECCC Meteorological Station download via GeoMet API...")
message("-------------------------------------------------------")

# --- 2. Source URL ---
BASE_URL <- "https://api.weather.gc.ca/collections/climate-stations/items"
url <- paste0(BASE_URL, "?limit=10000")
resp <- GET(url)
stop_for_status(resp)

stations_json <- content(resp, "text", encoding = "UTF-8")
parsed <- fromJSON(stations_json, simplifyVector = FALSE)
features <- parsed$features

stations_df <- map_dfr(features, function(f) {
  props <- f$properties
  props_clean <- lapply(props, function(x) if (is.null(x)) NA else x)
  as_tibble(props_clean)
})

# --- 3. Convert coordinates to decimal degrees ---
stations_df <- stations_df %>%
  mutate(
    LAT_DEC = LATITUDE / 1e7,
    LON_DEC = LONGITUDE / 1e7
  )

# --- 4. Filter Stations ---
nechako_met_stations <- stations_df %>%
  filter(
    PROV_STATE_TERR_CODE == "BC",
    LAT_DEC >= MIN_LAT,
    LAT_DEC <= MAX_LAT,
    LON_DEC >= MIN_LON,
    LON_DEC <= MAX_LON,
    !is.na(DLY_FIRST_DATE)
  ) %>%
  select(
    STATION_ID   = STN_ID,
    STATION_NAME,
    PROVINCE     = ENG_PROV_NAME,
    LATITUDE     = LAT_DEC,
    LONGITUDE    = LON_DEC,
    ELEVATION,
    FIRST_DATE,
    LAST_DATE,
    DLY_FIRST_DATE,
    DLY_LAST_DATE,
    HLY_FIRST_DATE,
    HLY_LAST_DATE,
    MLY_FIRST_DATE,
    MLY_LAST_DATE,
    HAS_HOURLY_DATA,
    HAS_MONTHLY_SUMMARY,
    HAS_NORMALS_DATA
  ) %>%
  arrange(desc(LATITUDE))

# --- 5. Save to CSV ---
write_csv(nechako_met_stations, output_file)

message("Found ", nrow(nechako_met_stations),
        " ECCC meteorological stations with daily data in the Nechako region.")
message("Station list successfully saved to: ", output_file)
message("-------------------------------------------------------")
print(head(nechako_met_stations))
