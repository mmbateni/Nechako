# =============================================================================
# Nechako Basin — Thiessen-Weighted Monthly & Annual Basin-Averaged Precipitation
# =============================================================================
#
# OVERVIEW
# --------
# This script:
#   1. Reads all NNN_CLIMATE_DAILY.csv station files from a folder
#   2. Builds Thiessen (Voronoi) polygons from station coordinates
#   3. Clips polygons to the basin boundary (shapefile OR auto convex hull)
#   4. Computes area-weighted fractions in an equal-area projection
#   5. Aggregates daily data → monthly totals with a completeness check
#   6. Applies Thiessen weights each month (re-normalised for missing stations)
#   7. Sums monthly → annual totals
#   8. Exports two CSV files: monthly and annual basin-averaged precip
#
# PROJECTION RATIONALE
# --------------------
# BC Albers (EPSG:3005) is used for ALL spatial operations.
# It is an Albers Equal-Area Conic projection designed specifically
# for British Columbia. It conserves AREA, which is the only metric
# that matters for Thiessen polygon weights.
#
# Projections deliberately AVOIDED:
#   • WGS84 / geographic (EPSG:4326) — lat/lon degrees; area wildly
#     distorted at ~54°N (cos distortion ~59%).
#   • UTM Zone 10N (EPSG:26910) — conformal (angle-preserving),
#     NOT area-preserving; area error grows quickly over a ~500 km basin.
#   • Any equidistant projection — preserves distances along specific
#     lines only; NOT suitable for area calculations.
#
# MISSING DATA HANDLING
# ---------------------
# A station-month is only included in the weighted average when the
# fraction of non-missing daily precip values >= completeness_threshold.
# If a station is excluded, its Thiessen weight is redistributed
# proportionally to the remaining valid stations (re-normalisation).
# Annual totals are the sum of 12 monthly basin averages.
# A 'months_available' column flags how many months contributed.
#
# REQUIRED PACKAGES
# -----------------
# install.packages(c("tidyverse", "lubridate", "sf"))
#
# =============================================================================

library(tidyverse)
library(lubridate)
library(sf)

# =============================================================================
# SECTION 1 — USER CONFIGURATION
# =============================================================================

# Folder containing the downloaded NNN_CLIMATE_DAILY.csv files
data_dir  <- "climate_daily_downloads"

# Station metadata file produced by Stations.r
stations_file <- "Stations.csv"

# Output folder (created if absent)
output_dir <- "thiessen_output"

# Path to a basin boundary polygon shapefile / GeoPackage.
# Set to NULL to auto-build from convex hull of station points + buffer.
#   Example: basin_shapefile <- "Nechako_watershed_boundary.shp"
basin_shapefile <- "D:/Nechako_Drought/Nechako/Spatial/nechakoBound_dissolve.shp"

# Buffer (metres) applied around the convex hull when basin_shapefile = NULL.
# Ensures edge stations are fully inside the boundary.
hull_buffer_m <- 15000   # 15 km — adjust if needed

# Minimum fraction of valid daily precip records required to keep
# a station-month as valid [0–1].  0.80 = at least 80% of days present.
completeness_threshold <- 0.80

# Equal-area CRS for all spatial operations
EQUAL_AREA_CRS <- 3005   # NAD83 / BC Albers (Albers Equal-Area Conic)

# =============================================================================
# SECTION 2 — SETUP
# =============================================================================

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

message("\n", strrep("=", 70))
message("  Nechako Basin — Thiessen-Weighted Precipitation Aggregation")
message(strrep("=", 70))

# =============================================================================
# SECTION 3 — LOAD STATION METADATA & MATCH TO AVAILABLE FILES
# =============================================================================

message("\n[1/6] Loading station metadata ...")

stations_meta <- read_csv(stations_file, show_col_types = FALSE) %>%
  filter(!is.na(STATION_ID), !is.na(LATITUDE), !is.na(LONGITUDE)) %>%
  # Some IDs come in as numeric — normalise to character for matching filenames
  mutate(STATION_ID = as.character(STATION_ID)) %>%
  distinct(STATION_ID, .keep_all = TRUE) %>%
  select(STATION_ID, STATION_NAME, LATITUDE, LONGITUDE)

# Discover all CSV data files in the data directory
csv_files <- tibble(
  filepath   = list.files(data_dir, pattern = "_CLIMATE_DAILY\\.csv$",
                          full.names = TRUE, ignore.case = TRUE),
  filename   = basename(filepath),
  STATION_ID = str_extract(filename, "^[0-9]+")
)

if (nrow(csv_files) == 0) {
  stop("No *_CLIMATE_DAILY.csv files found in: ", data_dir,
       "\nSet data_dir to the correct folder.")
}

# Inner-join: keep only stations with BOTH metadata AND a data file
stations_active <- inner_join(stations_meta, csv_files,
                              by = "STATION_ID")

message("  Stations with metadata: ", nrow(stations_meta))
message("  CSV files found:        ", nrow(csv_files))
message("  Matched (both):         ", nrow(stations_active))

if (nrow(stations_active) == 0) {
  stop("No station IDs matched between Stations.csv and the CSV files in ",
       data_dir, ".\nCheck that STATION_ID values match the filename prefixes.")
}

# =============================================================================
# SECTION 4 — BUILD SPATIAL OBJECTS & THIESSEN WEIGHTS
# =============================================================================

message("\n[2/6] Building Thiessen polygon weights (CRS: EPSG:", EQUAL_AREA_CRS, ") ...")

# --- 4a. Station points in equal-area CRS ---
stn_pts_geo <- st_as_sf(stations_active,
                        coords = c("LONGITUDE", "LATITUDE"),
                        crs    = 4326)                      # raw GPS coords
stn_pts     <- st_transform(stn_pts_geo, crs = EQUAL_AREA_CRS)

# --- 4b. Basin boundary ---
if (!is.null(basin_shapefile)) {
  message("  Loading basin boundary from: ", basin_shapefile)
  basin_boundary <- st_read(basin_shapefile, quiet = TRUE) %>%
    st_union() %>%
    st_transform(crs = EQUAL_AREA_CRS)
} else {
  message("  No shapefile supplied — building convex hull of stations + ",
          hull_buffer_m / 1000, " km buffer.")
  basin_boundary <- stn_pts %>%
    st_union() %>%
    st_convex_hull() %>%
    st_buffer(dist = hull_buffer_m)
}

# Ensure basin_boundary is a single-geometry sfc (for clipping)
basin_boundary <- st_union(basin_boundary)

# --- 4c. Voronoi (Thiessen) tessellation ---
# st_voronoi works on a MULTIPOINT geometry; envelope ensures polygons
# extend past the outermost stations before clipping.
envelope <- st_buffer(basin_boundary, dist = 50000)  # generous envelope

voronoi_raw <- stn_pts %>%
  st_union() %>%                              # combine into MULTIPOINT
  st_voronoi(envelope = envelope) %>%         # tessellate
  st_collection_extract("POLYGON") %>%        # pull individual polygons
  st_sf()                                     # convert to sf data frame

# --- 4d. Clip Voronoi cells to basin boundary ---
voronoi_clipped <- st_intersection(voronoi_raw, basin_boundary)

# --- 4e. Assign each Voronoi cell back to the nearest station point ---
# (st_voronoi does not preserve point order reliably, so we use nearest-join)
nearest_idx    <- st_nearest_feature(voronoi_clipped, stn_pts)

thiessen <- voronoi_clipped %>%
  mutate(
    STATION_ID   = stations_active$STATION_ID[nearest_idx],
    STATION_NAME = stations_active$STATION_NAME[nearest_idx],
    cell_area_m2 = as.numeric(st_area(geometry))
  )

# Check for duplicate assignments (shouldn't happen but guard anyway)
if (anyDuplicated(thiessen$STATION_ID)) {
  warning("Some stations were assigned more than one Voronoi cell — ",
          "merging them now.")
  thiessen <- thiessen %>%
    group_by(STATION_ID, STATION_NAME) %>%
    summarise(cell_area_m2 = sum(cell_area_m2), .groups = "drop")
}

# --- 4f. Compute weights ---
total_area <- sum(thiessen$cell_area_m2)

weights_df <- thiessen %>%
  st_drop_geometry() %>%
  mutate(thiessen_weight = cell_area_m2 / total_area) %>%
  select(STATION_ID, STATION_NAME, cell_area_m2, thiessen_weight)

# Sanity check: weights should sum to ~1
weight_sum <- sum(weights_df$thiessen_weight)
message(sprintf("  Basin area: %.0f km²   Weight sum: %.6f",
                total_area / 1e6, weight_sum))

if (abs(weight_sum - 1) > 0.001) {
  warning("Thiessen weights sum to ", round(weight_sum, 4),
          " (expected 1.000). Check basin boundary geometry.")
}

# Save weight table for reference
write_csv(weights_df %>%
            mutate(cell_area_km2 = cell_area_m2 / 1e6) %>%
            select(STATION_ID, STATION_NAME,
                   cell_area_km2, thiessen_weight),
          file.path(output_dir, "thiessen_weights.csv"))
message("  Weight table saved → thiessen_weights.csv")

# =============================================================================
# SECTION 5 — LOAD DAILY DATA & COMPUTE MONTHLY TOTALS
# =============================================================================

message("\n[3/6] Reading daily CSVs and computing monthly totals ...")
message("      (completeness threshold = ", completeness_threshold * 100, "% of days)")

# Helper: read one station file and return a tidy monthly summary
read_and_monthly <- function(filepath, station_id) {
  
  daily <- tryCatch(
    read_csv(filepath,
             col_types = cols(Date = col_date(format = "%Y-%m-%d"),
                              .default = col_double()),
             show_col_types = FALSE),
    error = function(e) {
      message("  WARNING: Could not read ", basename(filepath), " — skipping.")
      return(NULL)
    }
  )
  
  if (is.null(daily) || nrow(daily) == 0) return(NULL)
  
  # Ensure the Total_Precip_mm column exists
  if (!"Total_Precip_mm" %in% names(daily)) {
    # Try alternative names produced by different ECCC download vintages
    alt <- c("TOTAL_PRECIP_MM", "Total.Precip..mm.", "Precip..mm.")
    found <- intersect(alt, toupper(names(daily)))
    if (length(found) == 0) {
      message("  WARNING: No precipitation column in ", basename(filepath),
              " — skipping.")
      return(NULL)
    }
    daily <- daily %>%
      rename(Total_Precip_mm = names(daily)[match(found[1], toupper(names(daily)))])
  }
  
  daily %>%
    filter(!is.na(Date)) %>%
    mutate(
      Year  = year(Date),
      Month = month(Date),
      # Days in this calendar month (used for completeness denominator)
      days_in_month = days_in_month(Date)
    ) %>%
    group_by(Year, Month, days_in_month) %>%
    summarise(
      n_valid      = sum(!is.na(Total_Precip_mm)),   # days with a reading
      monthly_precip_mm = ifelse(
        n_valid / first(days_in_month) >= completeness_threshold,
        sum(Total_Precip_mm, na.rm = TRUE),           # monthly total
        NA_real_                                       # mark insufficient
      ),
      .groups = "drop"
    ) %>%
    mutate(STATION_ID = station_id) %>%
    select(STATION_ID, Year, Month, monthly_precip_mm, n_valid, days_in_month)
}

# Run for all matched stations
monthly_all <- map2_dfr(
  stations_active$filepath,
  stations_active$STATION_ID,
  read_and_monthly,
  .id = NULL
)

n_stns_read <- n_distinct(monthly_all$STATION_ID)
message("  Successfully read data from ", n_stns_read, " stations.")

# =============================================================================
# SECTION 6 — APPLY THIESSEN WEIGHTS — MONTHLY BASIN AVERAGE
# =============================================================================

message("\n[4/6] Applying Thiessen weights (re-normalised per month) ...")

# Join weights into the monthly panel
monthly_weighted <- monthly_all %>%
  inner_join(weights_df %>% select(STATION_ID, thiessen_weight),
             by = "STATION_ID")

# Determine full date range
all_years  <- sort(unique(monthly_weighted$Year))
all_months <- 1:12

# Build a complete Year × Month grid and fill with weighted basin average
basin_monthly <- expand_grid(
  Year  = seq(min(all_years), max(all_years)),
  Month = all_months
) %>%
  left_join(monthly_weighted, by = c("Year", "Month")) %>%
  group_by(Year, Month) %>%
  summarise(
    n_stations_total   = n_distinct(STATION_ID[!is.na(STATION_ID)]),
    n_stations_valid   = sum(!is.na(monthly_precip_mm)),
    sum_valid_weight   = sum(thiessen_weight[!is.na(monthly_precip_mm)],
                             na.rm = TRUE),
    basin_precip_mm    = ifelse(
      # Require at least 1 valid station AND combined weight > 0
      n_stations_valid > 0 & sum_valid_weight > 0,
      # Re-normalise: weighted sum / total valid weight
      sum(thiessen_weight[!is.na(monthly_precip_mm)] *
            monthly_precip_mm[!is.na(monthly_precip_mm)],
          na.rm = TRUE) / sum_valid_weight,
      NA_real_
    ),
    weight_coverage_pct = round(sum_valid_weight * 100, 1),  # % basin covered
    .groups = "drop"
  ) %>%
  mutate(
    YearMonth = sprintf("%04d-%02d", Year, Month)
  ) %>%
  select(Year, Month, YearMonth,
         basin_precip_mm, n_stations_valid,
         weight_coverage_pct)

message("  Monthly basin series: ", nrow(basin_monthly), " rows  (",
        min(all_years), "–", max(all_years), ")")

# =============================================================================
# SECTION 7 — ANNUAL TOTALS
# =============================================================================

message("\n[5/6] Computing annual basin totals ...")

basin_annual <- basin_monthly %>%
  group_by(Year) %>%
  summarise(
    months_available    = sum(!is.na(basin_precip_mm)),
    basin_annual_precip_mm = ifelse(
      months_available == 12,
      sum(basin_precip_mm, na.rm = TRUE),
      NA_real_                        # flag incomplete years
    ),
    # Partial sum even when incomplete (clearly labelled)
    basin_partial_sum_mm = sum(basin_precip_mm, na.rm = TRUE),
    min_weight_coverage_pct = min(weight_coverage_pct, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    complete_year = months_available == 12
  )

n_complete <- sum(basin_annual$complete_year, na.rm = TRUE)
message("  Annual records: ", nrow(basin_annual),
        "  (", n_complete, " complete years with all 12 months)")

# =============================================================================
# SECTION 8 — EXPORT RESULTS
# =============================================================================

message("\n[6/6] Exporting CSV files ...")

# --- Monthly output ---
out_monthly <- file.path(output_dir, "basin_precip_monthly.csv")

basin_monthly_out <- basin_monthly %>%
  rename(
    basin_precip_total_mm    = basin_precip_mm,
    n_valid_stations         = n_stations_valid,
    basin_weight_coverage_pct = weight_coverage_pct
  ) %>%
  arrange(Year, Month)

write_csv(basin_monthly_out, out_monthly)
message("  Monthly series  → ", out_monthly)

# --- Annual output ---
out_annual <- file.path(output_dir, "basin_precip_annual.csv")

basin_annual_out <- basin_annual %>%
  rename(
    annual_precip_mm          = basin_annual_precip_mm,
    partial_sum_mm            = basin_partial_sum_mm,
    min_basin_coverage_pct    = min_weight_coverage_pct
  ) %>%
  arrange(Year)

write_csv(basin_annual_out, out_annual)
message("  Annual series   → ", out_annual)

# =============================================================================
# SECTION 9 — QUICK CONSOLE SUMMARY
# =============================================================================

message("\n", strrep("=", 70))
message("  SUMMARY")
message(strrep("=", 70))

long_term_mean <- mean(basin_annual_out$annual_precip_mm, na.rm = TRUE)
long_term_sd   <- sd(basin_annual_out$annual_precip_mm,   na.rm = TRUE)

message(sprintf("  Basin area (km²)         : %.0f", total_area / 1e6))
message(sprintf("  Stations used            : %d", nrow(weights_df)))
message(sprintf("  Period                   : %d – %d",
                min(basin_annual_out$Year), max(basin_annual_out$Year)))
message(sprintf("  Complete years           : %d", n_complete))
message(sprintf("  Long-term mean precip    : %.1f mm/yr  (SD: %.1f mm)",
                long_term_mean, long_term_sd))
message(strrep("-", 70))
message("  Output files:")
message("    • ", out_monthly, "  (columns: Year, Month, YearMonth,")
message("         basin_precip_total_mm, n_valid_stations,")
message("         basin_weight_coverage_pct)")
message("    • ", out_annual,
        "  (columns: Year, annual_precip_mm [complete years only],")
message("         partial_sum_mm, months_available, complete_year)")
message("    • ", file.path(output_dir, "thiessen_weights.csv"),
        "  (station areas and weights)")
message(strrep("=", 70))
message("\n✅  Done.\n")