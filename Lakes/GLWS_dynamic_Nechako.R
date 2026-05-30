# =============================================================================
# Extract GLWS v1.1 Lake Data for Nechako River Basin
# =============================================================================

# Load required libraries
library(sf)        # For spatial data handling
library(dplyr)     # For data manipulation
library(readr)     # For reading CSV files
library(tidyr)     # For data tidying
library(ggplot2)   # For plotting

# =============================================================================
# ROOT PATH — update these two variables if the data moves
# =============================================================================

glws_root  <- "D:/Nechako_Drought/Nechako/Lakes/GLWS v1.1"
output_dir <- "D:/Nechako_Drought/Nechako/Lakes/GLWS v1.1/Results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# STEP 1: Load GLWS Shapefile Data
# =============================================================================

# Actual shapefile path (folder name contains spaces; .shp filename matches folder)
glws_shapefile_path <- file.path(
  glws_root,
  "Global database of lake water storage GLWS shapefile v1.1",
  "Global database of lake water storage GLWS shapefile v1.1.shp"
)

glws_lakes <- st_read(glws_shapefile_path)
print(st_drop_geometry(glws_lakes))

# =============================================================================
# STEP 2: Load Nechako River Basin Boundary
# =============================================================================

nechako_basin_path <- "D:/Nechako_Drought/Nechako/Spatial/nechakoBound_dissolve.shp"
nechako_basin <- st_read(nechako_basin_path)

# Option B: Load from HydroBASINS or other hydrological dataset
# hydrobasins <- st_read("path/to/HydroBASINS/shp/file.shp")
# nechako_basin <- hydrobasins %>% filter(PFAF_ID == "your_nechako_id")

nechako_basin <- st_transform(nechako_basin, st_crs(glws_lakes))

# =============================================================================
# STEP 3: Spatial Intersection - Find Lakes Within Nechako Basin
# =============================================================================

lakes_in_nechako <- st_intersection(glws_lakes, nechako_basin)
cat("Number of lakes in Nechako River basin:", nrow(lakes_in_nechako), "\n")

# =============================================================================
# STEP 4: Extract Lake Attributes
# =============================================================================

nechako_lake_data <- lakes_in_nechako %>%
  st_drop_geometry() %>%
  select(
    LakeID,
    LakeName,
    LakeArea,
    TypeName,
    Dryland,
    Trend,
    TrendPval,
    TrendErr,
    Notrend,
    Domidriver,
    IsRunoff,
    Modelagree,
    ModelR2,
    PrecSeason,   # truncated to 10 chars by ESRI Shapefile format (was PrepSeason)
    SediRate,
    SediLoss,
    SediLossEr,   # truncated to 10 chars by ESRI Shapefile format (was SediLosErr)
    NReservoir,
    TrendOri,
    TrendOriPv,
    Hydropower,
    InsituVali,
    ReguNatura
  )

print(nechako_lake_data)

# =============================================================================
# STEP 5: Save Shapefile
# Note: the attributes CSV is saved later in Step 6b, after time series
#       summary statistics are computed and joined to nechako_lake_data.
# =============================================================================

st_write(lakes_in_nechako, file.path(output_dir, "nechako_lakes_glws_v1.1.shp"),
         delete_dsn = TRUE)

# =============================================================================
# STEP 6: Load Time Series Data for Identified Lakes
# =============================================================================

timeseries_root <- file.path(
  glws_root,
  "Global database of lake water storage GLWS time series v1.1"
)

# -----------------------------------------------------------------------
# Method priority order by lake type:
#
#   Natural lakes : Poly -> SimpleLinear -> ConstantArea
#                   (LinearLi does NOT exist for natural lakes)
#                   (smaller lakes are often only in SimpleLinear)
#
#   Reservoirs    : Poly -> SimpleLinear -> ConstantArea
#
# The loader tries each method in order and uses the first one that has
# a matching file, recording which method was actually used.
# -----------------------------------------------------------------------

methods_natural_lake <- c("Poly", "SimpleLinear", "ConstantArea")
methods_reservoir    <- c("Poly", "SimpleLinear", "ConstantArea")

# -----------------------------------------------------------------------
# File naming convention in GLWS v1.1:
#   ID{LakeID}{LakeName}30m_GEE1992_2020_monthlyVolume_{Satellite}Levelrange.csv
# e.g. ID137Vygozero30m_GEE1992_2020_monthlyVolume_ICESat2Levelrange.csv
#
# Files are stored under:
#   Natural lakes  -> .../Natural lakes/{method}/
#   Reservoirs     -> .../Reservoirs/{method}/
# -----------------------------------------------------------------------

# Helper: search one method folder for a lake file (flat file or sub-directory)
find_lake_file <- function(lake_id, folder_path) {
  pattern   <- paste0("^ID", lake_id, "[^0-9]")
  all_files <- list.files(folder_path, pattern = pattern,
                          full.names = TRUE, recursive = FALSE)
  if (length(all_files) == 0) {
    subdirs      <- list.dirs(folder_path, recursive = FALSE)
    matching_dir <- subdirs[grepl(pattern, basename(subdirs))]
    if (length(matching_dir) > 0) {
      all_files <- list.files(matching_dir[1], pattern = "\\.csv$",
                              full.names = TRUE)
    }
  }
  return(all_files)
}

load_lake_timeseries <- function(lake_id, lake_type, ts_root) {
  
  is_reservoir <- grepl("reservoir", lake_type, ignore.case = TRUE)
  type_folder  <- if (is_reservoir) "Reservoirs" else "Natural lakes"
  methods      <- if (is_reservoir) methods_reservoir else methods_natural_lake
  
  for (method in methods) {
    folder_path <- file.path(ts_root, type_folder, method)
    all_files   <- find_lake_file(lake_id, folder_path)
    
    if (length(all_files) > 0) {
      data            <- read_csv(all_files[1], show_col_types = FALSE)
      data$LakeID     <- lake_id
      data$FileName   <- basename(all_files[1])
      data$MethodUsed <- method
      message(paste("Loaded LakeID:", lake_id, "| Type:", lake_type,
                    "| Method used:", method))
      return(data)
    }
  }
  
  message(paste("File not found for LakeID:", lake_id,
                "| Type:", lake_type,
                "| Tried:", paste(methods, collapse = " -> ")))
  return(NULL)
}

# Get list of LakeIDs and their types in Nechako basin
nechako_lake_ids   <- unique(nechako_lake_data$LakeID)
nechako_lake_types <- nechako_lake_data %>% select(LakeID, TypeName)

# Load time series for all lakes (auto-selects best available method)
timeseries_list <- lapply(nechako_lake_ids, function(id) {
  ltype <- nechako_lake_types$TypeName[nechako_lake_types$LakeID == id][1]
  load_lake_timeseries(id, ltype, timeseries_root)
})

# Combine all time series data
nechako_timeseries <- bind_rows(timeseries_list)

# Save combined time series
write_csv(nechako_timeseries,
          file.path(output_dir, "nechako_lakes_timeseries_glws_v1.1.csv"))

# =============================================================================
# STEP 6b: Summarise time series per lake and append to attributes table
# =============================================================================

# The Time column is in YYYYDDD format (year + 3-digit day of year).
# Extract the year as the first 4 digits.
#
# rws = water storage anomaly relative to the first record (Gt).
# Only non-NA rws rows are counted and averaged.

timeseries_summary <- nechako_timeseries %>%
  mutate(Year = as.integer(substr(as.character(Time), 1, 4))) %>%
  filter(!is.na(rws)) %>%
  group_by(LakeID) %>%
  summarise(
    Avg_rws_Gt   = mean(rws),                 # mean storage anomaly (Gt)
    N_obs        = n(),                        # total monthly records available
    Year_start   = min(Year),                  # first year with data
    Year_end     = max(Year),                  # last year with data
    Time_span_yr = Year_end - Year_start + 1,  # inclusive span in years
    MethodUsed   = first(MethodUsed),          # hypsometry method loaded
    .groups = "drop"
  )

# Join summary columns onto the lake attributes table
nechako_lake_data <- nechako_lake_data %>%
  left_join(timeseries_summary, by = "LakeID")

# Preview the enriched table
print(nechako_lake_data %>%
        select(LakeID, LakeName, TypeName, Avg_rws_Gt,
               N_obs, Year_start, Year_end, Time_span_yr, MethodUsed))

# Save the enriched attributes CSV (includes time series summary columns)
write_csv(nechako_lake_data,
          file.path(output_dir, "nechako_lakes_glws_v1.1.csv"))

cat("\nOutputs saved to:", output_dir, "\n")

# =============================================================================
# STEP 7: Alternative - Using Lake Coordinates CSV (If No Shapefile Access)
# =============================================================================

coords_path <- file.path(glws_root, "GLWS lake coordinates v1.1.csv")
lake_coords <- read_csv(coords_path, show_col_types = FALSE)

# Column names are lowercase: latitude / longitude
lake_coords_sf <- st_as_sf(lake_coords,
                           coords = c("longitude", "latitude"),
                           crs = 4326)

lake_coords_sf <- st_transform(lake_coords_sf, st_crs(nechako_basin))

lakes_in_nechako_coords <- st_join(lake_coords_sf, nechako_basin,
                                   join = st_within)

lakes_in_nechako_coords <- lakes_in_nechako_coords %>%
  filter(!is.na(LakeID))

# =============================================================================
# STEP 8: Basin Map
# =============================================================================

ggplot() +
  geom_sf(data = nechako_basin, fill = "lightblue", alpha = 0.5) +
  geom_sf(data = lakes_in_nechako, color = "darkblue", size = 2) +
  labs(title = "GLWS Lakes in Nechako River Basin",
       subtitle = paste("Number of lakes:", nrow(lakes_in_nechako))) +
  theme_minimal()

ggsave(file.path(output_dir, "nechako_lakes_map.png"),
       width = 10, height = 8, dpi = 300)

# =============================================================================
# STEP 9: Per-lake water storage anomaly time series — plots and CSVs
# =============================================================================

# GLWS provides water storage *anomaly* (rws) relative to the first record,
# not absolute storage. Absolute storage values are not available in GLWS v1.1.
#
# The Time column is in YYYYDDD format. We convert it to a proper Date by
# treating the DDD part as the day-of-year, giving an approximate monthly date.
#
# rws_wSediAdj is used where available (sedimentation-corrected anomaly).
# If that column is absent or all-NA for a lake, rws is used as fallback.
#
# Output per lake (saved to output_dir/TimeSeries_per_lake/):
#   - ID{LakeID}_{LakeName}_rws_timeseries.csv   (time series data)
#   - ID{LakeID}_{LakeName}_rws_timeseries.png   (plot)

ts_output_dir <- file.path(output_dir, "TimeSeries_per_lake")
dir.create(ts_output_dir, showWarnings = FALSE, recursive = TRUE)

# Helper: convert YYYYDDD integer to Date
yyyyddd_to_date <- function(yyyyddd) {
  yr  <- as.integer(substr(as.character(yyyyddd), 1, 4))
  doy <- as.integer(substr(as.character(yyyyddd), 5, 7))
  as.Date(paste(yr, doy), format = "%Y %j")
}

# Build a lookup of LakeID -> display name (use LakeName if available, else ID)
lake_names_lookup <- nechako_lake_data %>%
  mutate(DisplayName = ifelse(is.na(LakeName) | LakeName == "",
                              paste0("ID", LakeID),
                              paste0("ID", LakeID, "_", LakeName))) %>%
  select(LakeID, DisplayName, TypeName)

# Loop over each lake that has data loaded
lakes_with_data <- unique(nechako_timeseries$LakeID)

for (lid in lakes_with_data) {
  
  lake_ts   <- nechako_timeseries %>% filter(LakeID == lid)
  meta      <- lake_names_lookup %>% filter(LakeID == lid)
  dname     <- meta$DisplayName
  lake_type <- meta$TypeName
  method    <- unique(lake_ts$MethodUsed)[1]
  
  # Convert Time to Date
  lake_ts <- lake_ts %>%
    mutate(Date = yyyyddd_to_date(Time))
  
  # Choose storage variable:
  # Prefer rws_wSediAdj (sedimentation-corrected) for reservoirs if available,
  # otherwise fall back to rws (anomaly without sedimentation correction).
  use_sedi_adj <- "rws_wSediAdj" %in% names(lake_ts) &&
    any(!is.na(lake_ts$rws_wSediAdj))
  
  if (use_sedi_adj) {
    lake_ts  <- lake_ts %>% mutate(Storage_Gt = rws_wSediAdj)
    y_label  <- "Water Storage Anomaly (Gt)\n[sedimentation-corrected: rws_wSediAdj]"
    var_note <- "rws_wSediAdj"
  } else {
    lake_ts  <- lake_ts %>% mutate(Storage_Gt = rws)
    y_label  <- "Water Storage Anomaly (Gt)\n[rws, relative to first record]"
    var_note <- "rws"
  }
  
  # --- Save individual CSV ---------------------------------------------------
  csv_cols <- c("Date", "Time", "Mean_Levels", "Mean_Level_Errors",
                "Cleaned_Area", "rws", "rws_err", "Storage_Gt", "MethodUsed")
  csv_cols <- csv_cols[csv_cols %in% names(lake_ts)]
  
  write_csv(
    lake_ts %>% select(all_of(csv_cols)),
    file.path(ts_output_dir, paste0(dname, "_rws_timeseries.csv"))
  )
  
  # --- Plot ------------------------------------------------------------------
  ts_clean <- lake_ts %>% filter(!is.na(Storage_Gt))
  
  if (nrow(ts_clean) == 0) {
    message(paste("No valid rws data to plot for", dname))
    next
  }
  
  p <- ggplot(ts_clean, aes(x = Date, y = Storage_Gt)) +
    geom_line(colour = "steelblue", linewidth = 0.7) +
    geom_point(colour = "steelblue", size = 1.2, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    labs(
      title    = paste("Water Storage Anomaly —", gsub("_", " ", dname)),
      subtitle = paste0("Type: ", lake_type,
                        "  |  Method: ", method,
                        "  |  Variable: ", var_note,
                        "  |  N = ", nrow(ts_clean), " months"),
      x        = "Date",
      y        = y_label,
      caption  = "Source: GLWS v1.1 (Yao et al., 2023)"
    ) +
    scale_x_date(date_breaks = "3 years", date_labels = "%Y") +
    theme_bw() +
    theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, colour = "grey30"),
      axis.text.x   = element_text(angle = 45, hjust = 1),
      plot.caption  = element_text(size = 7, colour = "grey50")
    )
  
  ggsave(
    file.path(ts_output_dir, paste0(dname, "_rws_timeseries.png")),
    plot = p, width = 10, height = 5, dpi = 300
  )
  
  message(paste("Saved plot and CSV for", dname))
}

cat("\nPer-lake time series saved to:", ts_output_dir, "\n")

# =============================================================================
# NOTES & IMPORTANT INFORMATION
# =============================================================================

# 1. Two paths are set at the top of the script:
#    glws_root  <- input data root
#    output_dir <- where all outputs are saved (CSV, SHP, PNG)
#    output_dir is created automatically if it does not exist.

# 2. Coordinate columns in GLWS_lake_coordinates_v1_1.csv are LOWERCASE:
#    "latitude" and "longitude" (not capitalised).

# 3. Time series filenames follow the pattern:
#    ID{LakeID}{LakeName}30m_GEE1992_2020_monthlyVolume_{Satellite}Levelrange.csv
#    Files live under:  .../Natural lakes/{method}/  or  .../Reservoirs/{method}/
#    The loader finds the right file via ID-prefix matching.

# 4. Hypsometry method selection (Step 6):
#    Natural lakes : Poly -> SimpleLinear -> ConstantArea
#                    (LinearLi does NOT exist for natural lakes)
#    Reservoirs    : Poly -> SimpleLinear -> ConstantArea
#    The first method with a matching file is used. "MethodUsed" column
#    in the output records which method was actually loaded.
#    To change priority, edit methods_natural_lake / methods_reservoir.

# 5. Summary columns added to nechako_lakes_glws_v1.1.csv (Step 6b):
#    Avg_rws_Gt   - mean water storage anomaly across all available months (Gt)
#    N_obs        - total number of monthly observations available
#    Year_start   - first year with data in GLWS
#    Year_end     - last year with data in GLWS
#    Time_span_yr - inclusive number of years covered (Year_end - Year_start + 1)
#    MethodUsed   - hypsometry method the time series was loaded from
#    NA in these columns means no time series file was found for that lake.

# 6. Per-lake time series outputs (Step 9) are saved to:
#    output_dir/TimeSeries_per_lake/
#    GLWS v1.1 provides storage *anomaly* (rws), not absolute storage.
#    For reservoirs with sedimentation data, rws_wSediAdj is plotted instead.
#    The subtitle of each plot states which variable and method were used.

# 7. ESRI Shapefile 10-character field name truncations:
#    PrepSeason -> PrecSeason
#    SediLosErr -> SediLossEr

# 8. Citation for GLWS v1.1:
#    Yao, F., et al. (2023). Satellites reveal widespread decline in global
#    lake water storage. Science. https://doi.org/10.1126/science.abo2812