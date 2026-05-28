# =============================================================================
# SGI (Standardised Groundwater level Index) – Nechako River Basin
# Method: Bloomfield & Marchant (2013), doi:10.5194/hess-17-4769-2013
# Data:   GROW 1.0 (Bäthge et al. 2026)
# Output: One median SGI value per month across all qualifying wells
# =============================================================================

# ---- 0. Packages ------------------------------------------------------------
# install.packages(c("arrow","sf","tidyverse","lubridate","qpdf"))
library(arrow)       # read parquet files efficiently
library(sf)          # spatial operations
library(tidyverse)   # dplyr, tidyr, ggplot2, purrr, readr
library(lubridate)   # date helpers


# ---- 1. Configuration -------------------------------------------------------
cfg <- list(

  # Paths to the two GROW 1.0 files (edit to your local copies)
  path_attributes = "grow_attributes.parquet",
  path_timeseries = "grow_timeseries.parquet",

  # Basin boundary (KMZ = zipped KML; sf can read KML directly after unzip)
  basin_kmz = "Spatial/nechakoBound_dissolve.kmz",

  # Quality filters (adjust as needed)
  min_years     = 10,    # minimum record length in years
  max_gap_frac  = 0.20,  # maximum fraction of gaps allowed
  use_filled    = TRUE,  # TRUE  → use linearly gap-filled columns
                         # FALSE → use raw columns (NA where gaps exist)

  # Time series intervals to keep.
  # GROW stores: "d" = daily, "MS" = monthly, "YS" = yearly.
  # SGI requires monthly values; daily series will be aggregated to monthly.
  keep_intervals = c("d", "MS"),   # exclude yearly ("YS") – too coarse for SGI

  # Target output: a character vector of months you want SGI for.
  # Format "YYYY-MM". Leave NULL to compute all available months.
  target_months = NULL   # e.g. c("2010-07", "2015-08") or NULL for all
)


# ---- 2. Helper: read basin polygon from KMZ ---------------------------------
read_kmz <- function(kmz_path) {
  # KMZ is a ZIP archive containing a .kml file; sf cannot open it directly
  # on all platforms, so we unzip first.
  tmp_dir <- tempdir()
  utils::unzip(kmz_path, exdir = tmp_dir)
  kml_file <- list.files(tmp_dir, pattern = "\\.kml$", full.names = TRUE)[1]
  if (is.na(kml_file)) stop("No .kml found inside ", kmz_path)
  # sf::st_read reads all layers; we want the polygon layer
  layers <- sf::st_layers(kml_file)$name
  # prefer the largest polygon layer
  poly_layer <- layers[1]
  sf::st_read(kml_file, layer = poly_layer, quiet = TRUE) |>
    sf::st_transform(crs = 4326) |>       # ensure WGS-84
    sf::st_union() |>                      # dissolve sub-polygons if any
    sf::st_make_valid()
}


# ---- 3. Helper: normal-scores transform (Bloomfield & Marchant 2013) --------
# For n observations, pi_k = (2*rank_k - 1) / (2*n), k = 1 … n
# SGI_k = Φ⁻¹(pi_k)
# Ties get the same rank (average method).
# NAs are ignored in ranking and returned as NA in output.
normal_scores <- function(x) {
  n     <- sum(!is.na(x))
  if (n < 4) return(rep(NA_real_, length(x)))   # too few obs → skip
  r     <- rank(x, na.last = "keep", ties.method = "average")
  pi_k  <- (2 * r - 1) / (2 * n)
  qnorm(pi_k)   # inverse standard-normal CDF
}


# ---- 4. Load and spatially filter attributes --------------------------------
message("Loading GROW attributes …")
attrs_raw <- arrow::read_parquet(cfg$path_attributes)

# Spatial subset: keep wells inside the Nechako basin -------------------------
message("Reading basin boundary from KMZ …")
basin <- read_kmz(cfg$basin_kmz)

attrs_sf <- attrs_raw |>
  filter(!is.na(latitude), !is.na(longitude)) |>
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326,
               remove = FALSE)

in_basin <- sf::st_within(attrs_sf, basin, sparse = FALSE)[, 1]
attrs_basin <- attrs_sf[in_basin, ] |> sf::st_drop_geometry()

message(sprintf("  Wells with coordinates inside basin: %d", nrow(attrs_basin)))

# Quality filters --------------------------------------------------------------
attrs_filtered <- attrs_basin |>
  filter(
    interval  %in% cfg$keep_intervals,
    length_years >= cfg$min_years,
    gap_fraction  <= cfg$max_gap_frac
  )

message(sprintf("  Wells after quality filters: %d", nrow(attrs_filtered)))

if (nrow(attrs_filtered) == 0) {
  stop("No wells pass the quality filters. ",
       "Try relaxing min_years or max_gap_frac in cfg.")
}

keep_ids <- attrs_filtered$GROW_ID


# ---- 5. Load time series for qualifying wells only --------------------------
message("Loading time series (filtered by GROW_ID) …")

# arrow's lazy evaluation lets us push the filter to disk before reading
ts_raw <- arrow::open_dataset(cfg$path_timeseries) |>
  filter(GROW_ID %in% keep_ids) |>
  select(GROW_ID, interval, date, month,
         # raw columns
         groundwater_depth_from_ground_elevation_m,
         groundwater_water_level_elevation_m_asl,
         # gap-filled columns
         groundwater_filled_depth_from_ground_elevation_m,
         groundwater_filled_water_level_elevation_m_asl) |>
  collect()   # pull into memory

message(sprintf("  Records loaded: %d", nrow(ts_raw)))


# ---- 6. Harmonise groundwater column ----------------------------------------
# Each time series carries groundwater as EITHER depth-from-ground OR
# water-level-elevation (the other two columns are NA).
# For SGI we need a consistent direction: "higher value = more water".
# Convention used here:
#   • water_level_elevation → keep as-is (higher = wetter)
#   • depth_from_ground     → negate (smaller depth = wetter, so ×-1)
# This matches the sign convention in Bloomfield & Marchant (2013) and the
# slope inversion done in func_processing_gw_time_series.py.

ts <- ts_raw |>
  mutate(
    gw_raw = case_when(
      !is.na(groundwater_water_level_elevation_m_asl)
        ~  groundwater_water_level_elevation_m_asl,
      !is.na(groundwater_depth_from_ground_elevation_m)
        ~ -groundwater_depth_from_ground_elevation_m,
      TRUE ~ NA_real_
    ),
    gw_filled = case_when(
      !is.na(groundwater_filled_water_level_elevation_m_asl)
        ~  groundwater_filled_water_level_elevation_m_asl,
      !is.na(groundwater_filled_depth_from_ground_elevation_m)
        ~ -groundwater_filled_depth_from_ground_elevation_m,
      TRUE ~ NA_real_
    ),
    gw = if (cfg$use_filled) gw_filled else gw_raw
  )


# ---- 7. Aggregate daily series to monthly means -----------------------------
ts_monthly <- ts |>
  mutate(
    month_label = if_else(
      interval == "d",
      format(as.Date(date), "%Y-%m"),  # derive month from daily date
      as.character(month)              # already monthly in GROW
    )
  ) |>
  group_by(GROW_ID, month_label) |>
  summarise(
    gw        = mean(gw, na.rm = TRUE),
    n_records = sum(!is.na(gw)),
    .groups   = "drop"
  ) |>
  filter(!is.na(gw))   # drop months that ended up all-NA


# ---- 8. Compute SGI per well using the normal-scores transform --------------
# For each well × calendar month (1–12) separately, as in Bloomfield & Marchant.

message("Computing SGI …")

sgi_long <- ts_monthly |>
  mutate(
    cal_month = month(ym(month_label))   # integer 1–12
  ) |>
  group_by(GROW_ID, cal_month) |>
  mutate(
    sgi = normal_scores(gw)
  ) |>
  ungroup()


# ---- 9. Aggregate to one basin-wide value per calendar month ----------------
# Median across all qualifying wells (robust to outliers and unequal coverage)

sgi_basin <- sgi_long |>
  group_by(month_label) |>
  summarise(
    sgi_median = median(sgi, na.rm = TRUE),
    sgi_mean   = mean(sgi,   na.rm = TRUE),
    sgi_sd     = sd(sgi,     na.rm = TRUE),
    n_wells    = sum(!is.na(sgi)),
    .groups    = "drop"
  ) |>
  arrange(month_label) |>
  mutate(
    date        = ym(month_label),         # proper Date for plotting
    # Drought classification after McKee et al. (1993)
    drought_cat = case_when(
      sgi_median >  0.00               ~ "No drought",
      sgi_median >= -1.00              ~ "Minor drought",
      sgi_median >= -1.50              ~ "Moderate drought",
      sgi_median >= -2.00              ~ "Severe drought",
      TRUE                             ~ "Extreme drought"
    )
  )

# Optional: subset to requested months
if (!is.null(cfg$target_months)) {
  sgi_basin <- sgi_basin |> filter(month_label %in% cfg$target_months)
}


# ---- 10. Print results -------------------------------------------------------
cat("\n========================================================\n")
cat("  SGI – Nechako River Basin (median across wells)\n")
cat(sprintf("  Wells used: %d | min_years≥%d | gap_frac≤%.2f\n",
            n_distinct(sgi_long$GROW_ID),
            cfg$min_years, cfg$max_gap_frac))
cat("========================================================\n\n")

print(sgi_basin |>
        select(month_label, sgi_median, n_wells, drought_cat),
      n = 50)


# ---- 11. Save CSV output ----------------------------------------------------
output_path <- "sgi_nechako_monthly.csv"
write_csv(sgi_basin, output_path)
message(sprintf("\nResults saved to: %s", output_path))


# ---- 12. Quick diagnostic plot ----------------------------------------------
# Requires ggplot2 (part of tidyverse)

p <- ggplot(sgi_basin, aes(x = date, y = sgi_median)) +

  # Shade drought episodes (SGI < 0) in red, wet in blue
  geom_ribbon(
    aes(ymin = pmin(sgi_median, 0), ymax = 0),
    fill = "#d73027", alpha = 0.4
  ) +
  geom_ribbon(
    aes(ymin = 0, ymax = pmax(sgi_median, 0)),
    fill = "#4575b4", alpha = 0.4
  ) +

  geom_line(linewidth = 0.6, colour = "grey30") +
  geom_hline(yintercept = 0,   linetype = "solid",  colour = "black", linewidth = 0.4) +
  geom_hline(yintercept = -1,  linetype = "dashed", colour = "#fc8d59", linewidth = 0.4) +
  geom_hline(yintercept = -1.5,linetype = "dashed", colour = "#d73027", linewidth = 0.4) +
  geom_hline(yintercept = -2,  linetype = "dashed", colour = "#a50026", linewidth = 0.4) +

  # Annotate n_wells as secondary information
  geom_text(
    data = sgi_basin |> slice_min(date),
    aes(label = paste0("n wells = ", n_wells)),
    hjust = 0, vjust = -0.5, size = 3, colour = "grey40"
  ) +

  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_y_continuous(
    breaks = c(-2, -1.5, -1, 0, 1, 2),
    labels = c("-2\n(Extreme)", "-1.5\n(Severe)", "-1\n(Moderate)",
               "0", "1", "2")
  ) +
  labs(
    title    = "Standardised Groundwater level Index (SGI)\nNechako River Basin",
    subtitle = sprintf(
      "Median across %d wells  |  min_years ≥ %d  |  gap_fraction ≤ %.2f  |  %s groundwater",
      n_distinct(sgi_long$GROW_ID),
      cfg$min_years, cfg$max_gap_frac,
      if (cfg$use_filled) "gap-filled" else "raw"
    ),
    x        = NULL,
    y        = "SGI",
    caption  = "Data: GROW 1.0 (Bäthge et al. 2026) | Method: Bloomfield & Marchant (2013)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave("sgi_nechako_plot.png", p, width = 10, height = 4, dpi = 150)
message("Plot saved to: sgi_nechako_plot.png")

# ---- END -------------------------------------------------------------------
