#=============================================================================
# SGI (Standardised Groundwater level Index) – Nechako River Basin
# Method: Bloomfield & Marchant (2013), doi:10.5194/hess-17-4769-2013
# Multi-source groundwater data integration:
#   1. BC PGOWN (primary) via Environmental Reporting BC
#   2. Federal: NRCan GIN/NGWD + ECCC WSC (tidyhydat)
#   3. BC Real-Time Water Data (ArcGIS FeatureServer)
#   4. BC GWELLS API (static levels + occasional time series)
# Output: One median SGI value per month across all qualifying wells
#=============================================================================

#---- 0. Packages ------------------------------------------------------------
pkgs <- c("httr", "sf", "tidyverse", "lubridate", "jsonlite", "ckanr", "tidyhydat", "geosphere")
new <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if(length(new)) install.packages(new)
invisible(lapply(pkgs, library, character.only = TRUE))

#---- 1. Configuration -------------------------------------------------------
cfg <- list(
  ## === Source toggles ===
  src_pgown      = TRUE,   # BC PGOWN monthly medians (default primary)
  src_gin_ngwd   = FALSE,  # NRCan Groundwater Information System
  src_wsc        = FALSE,  # Water Survey of Canada (tidyhydat)
  src_bc_realtime= FALSE,  # BC Real-Time Water Data (ArcGIS REST)
  src_gwells     = FALSE,  # BC GWELLS API (static + sparse TS)
  
  ## === PGOWN (primary) ===
  url_pgown_ts   = "https://catalogue.data.gov.bc.ca/dataset/a74f1b97-17f7-499b-84e7-6455e169e425/resource/35162de7-803e-42f2-93d4-563330bb7dc9/download/monthlywell_ts.csv",
  url_pgown_attr = "https://catalogue.data.gov.bc.ca/dataset/a74f1b97-17f7-499b-84e7-6455e169e425/resource/a8933793-eadb-4a9c-992c-da4f6ac8ca51/download/gw_well_table.csv",
  
  ## === GIN/NGWD (Federal) ===
  gin_query_bbox = c(-126.5, -53.5, -124.0, -52.0),
  gin_cache      = "gin_nechako_subset.csv",
  
  ## === WSC (tidyhydat) ===
  wsc_radius_km  = 50,
  wsc_daily_agg  = "median",
  
  ## === BC Real-Time (ArcGIS) ===
  bc_rt_layer    = "https://maps.gov.bc.ca/arcserver/rest/services/pub/Water/BC_Water_Levels/MapServer/1",
  bc_rt_cache    = "bc_rt_groundwater.csv",
  
  ## === GWELLS API ===
  gwells_api_base= "https://apps.nrs.gov.bc.ca/gwells/api/v2",
  gwells_min_obs = 5,
  
  ## === Basin boundary ===
  basin_kmz = "Spatial/nechakoBound_dissolve.kmz",
  
  ## === Caching ===
  cache_dir = "cache/",
  
  ## === Quality filters ===
  min_years    = 10,
  max_gap_frac = 0.20,
  min_freq     = "monthly",
  
  ## === Output ===
  target_months = NULL,
  output_prefix = "sgi_nechako_multisource"
)

if(!dir.exists(cfg$cache_dir)) dir.create(cfg$cache_dir, recursive=TRUE)

#---- 2. Helpers: Basin, Fetch, Transform ------------------------------------
read_kmz <- function(kmz_path) {
  tmp_dir <- tempdir()
  utils::unzip(kmz_path, exdir=tmp_dir)
  kml_file <- list.files(tmp_dir, "\\.kml$", full.names=TRUE)[1]
  if(is.na(kml_file)) stop("No .kml in ", kmz_path)
  lyr <- sf::st_layers(kml_file)$name[1]
  sf::st_read(kml_file, layer=lyr, quiet=TRUE) |> 
    sf::st_transform(4326) |> 
    sf::st_union() |> 
    sf::st_make_valid()
}

fetch_csv <- function(url, cache_path, desc="") {
  if(!is.null(cache_path) && file.exists(cache_path)) {
    message(sprintf("  [cache] %s", cache_path))
    return(read_csv(cache_path, show_col_types=FALSE))
  }
  message(sprintf("  [download] %s …", desc))
  dest <- if(!is.null(cache_path)) cache_path else tempfile(fileext=".csv")
  resp <- GET(url, write_disk(dest, overwrite=TRUE), progress(), timeout(120))
  stop_for_status(resp, task=paste("download", desc))
  read_csv(dest, show_col_types=FALSE)
}

normal_scores <- function(x) {
  n <- sum(!is.na(x))
  if(n<4) return(rep(NA_real_, length(x)))
  r <- rank(x, na.last="keep", ties.method="average")
  qnorm((2*r - 1)/(2*n))
}

# Standardize any source dataframe to: well_num, date, gw (m, neg=deeper), source_id
standardize_ts <- function(df, src, well_col, date_col, val_col, negate=TRUE) {
  df |> 
    rename(well_num = all_of(well_col), date = all_of(date_col), gw_raw = all_of(val_col)) |>
    mutate(
      well_num = as.character(well_num),
      date = as.Date(date),
      gw = if(negate) -as.numeric(gw_raw) else as.numeric(gw_raw),
      source_id = src,
      month_label = format(date, "%Y-%m")
    ) |>
    filter(!is.na(date), !is.na(gw)) |>
    dplyr::select(well_num, date, gw, month_label, source_id)
}

#---- 3. Spatial Filter Setup ------------------------------------------------
message("Loading basin boundary …")
basin <- read_kmz(cfg$basin_kmz)
basin_centroid <- sf::st_centroid(basin) |> sf::st_coordinates()

#---- 4. Data Source Handlers ------------------------------------------------

## 4.1 PGOWN (Primary) ------------------------------------------------------
load_pgown <- function() {
  if(!cfg$src_pgown) return(NULL)
  message("Loading PGOWN (Environmental Reporting BC) …")
  
  detect_col <- function(df, patterns, label) {
    m <- grep(paste(patterns, collapse="|"), names(df), ignore.case=TRUE, value=TRUE)
    if(!length(m)) stop(sprintf("Cannot find '%s' column. Tried: %s", label, paste(patterns, collapse=", ")))
    m[1]
  }
  
  attrs_raw <- fetch_csv(cfg$url_pgown_attr, file.path(cfg$cache_dir, "pgown_attrs.csv"), "PGOWN attributes") |>
    rename_with(tolower)
  
  idx_cols <- grep("^\\.\\.\\.\\d+$", names(attrs_raw), value=TRUE)
  if(length(idx_cols)) attrs_raw <- dplyr::select(attrs_raw, -dplyr::all_of(idx_cols))
  
  col_well <- detect_col(attrs_raw, c("^well_num$", "well_id", "obs.*well"), "well ID")
  col_lat  <- detect_col(attrs_raw, c("^lat$", "^latitude$"), "latitude")
  col_lon  <- detect_col(attrs_raw, c("^lon$", "^long$", "^longitude$"), "longitude")
  
  attrs <- attrs_raw |>
    rename(well_num = all_of(col_well), lat = all_of(col_lat), lon = all_of(col_lon)) |>
    mutate(well_num = as.character(well_num)) |>
    filter(!is.na(lat), !is.na(lon))
  
  wells_sf <- st_as_sf(attrs, coords=c("lon","lat"), crs=4326)
  keep <- st_within(wells_sf, basin, sparse=FALSE)[,1]
  keep_wells <- attrs[keep, "well_num"]
  
  ts_raw <- fetch_csv(cfg$url_pgown_ts, file.path(cfg$cache_dir, "pgown_ts.csv"), "PGOWN monthly") |>
    rename_with(tolower)
  
  col_wn  <- detect_col(ts_raw, c("^well_num$", "well_id"), "well ID")
  col_dt  <- detect_col(ts_raw, c("^date$", "datetime"), "date")
  col_val <- detect_col(ts_raw, c("med_gwl", "gwl", "value", "level"), "GWL")
  
  ts <- ts_raw |>
    rename(well_num = all_of(col_wn), date = all_of(col_dt), gw_raw = all_of(col_val)) |>
    mutate(
      well_num = as.character(well_num),
      date = as.Date(date),
      gw = -as.numeric(gw_raw),
      month_label = format(date, "%Y-%m"),
      source_id = "PGOWN"
    ) |>
    filter(!is.na(date), !is.na(gw), well_num %in% keep_wells) |>
    dplyr::select(dplyr::all_of(c("well_num", "date", "gw", "month_label", "source_id")))  # FIXED
  
  message(sprintf("  PGOWN: %d wells, %d records", n_distinct(ts$well_num), nrow(ts)))
  ts
}

## 4.2 GIN/NGWD (Federal) ---------------------------------------------------
load_gin <- function() {
  if(!cfg$src_gin_ngwd) return(NULL)
  message("Loading GIN/NGWD (NRCan) …")
  if(file.exists(cfg$gin_cache)) {
    df <- read_csv(cfg$gin_cache, show_col_types=FALSE)
  } else {
    message("  [manual] Download GIN subset for Nechako bbox, save to ", cfg$gin_cache)
    return(NULL)
  }
  df_std <- df |> 
    rename_with(tolower) |>
    standardize_ts("GIN", 
                   well_col = matches("station_id|well_num|obs_id"),
                   date_col = matches("date|datetime"),
                   val_col = matches("level_m|water_level|gwl"),
                   negate = FALSE)
  
  # Add lat/lon for spatial filter (flexible column detection)
  lat_col <- grep("lat|latitude", names(df), ignore.case=TRUE, value=TRUE)[1]
  lon_col <- grep("lon|longitude", names(df), ignore.case=TRUE, value=TRUE)[1]
  
  df_std |>
    mutate(lat = as.numeric(df[[lat_col]]),
           lon = as.numeric(df[[lon_col]])) |>
    st_as_sf(coords=c("lon","lat"), crs=4326, remove=FALSE) |>
    filter(st_within(., basin, sparse=FALSE)[,1]) |>
    as_tibble() |>
    dplyr::select(well_num, date, gw, month_label, source_id)
}

## 4.3 WSC via tidyhydat ----------------------------------------------------
load_wsc <- function() {
  if(!cfg$src_wsc) return(NULL)
  message("Loading WSC groundwater (tidyhydat) …")
  stations <- tidyhydat::hydrometric_stations() |>
    filter(STATION_TYPE == "Groundwater") |>
    mutate(dist_km = geosphere::distHaversine(
      c(LONGITUDE, LATITUDE), c(basin_centroid[1], basin_centroid[2]))/1000) |>
    filter(dist_km <= cfg$wsc_radius_km)
  
  if(nrow(stations)==0) { message("  No WSC groundwater stations nearby"); return(NULL) }
  
  ts_list <- lapply(stations$STATION_NUMBER, function(id) {
    tryCatch({
      tidyhydat::hy_daily_levels(id) |>
        mutate(well_num = id, source_id = "WSC") |>
        filter(!is.na(Value)) |>
        group_by(well_num, month_label = format(Date, "%Y-%m")) |>
        summarise(gw = if(cfg$wsc_daily_agg=="median") median(Value) else mean(Value),
                  date = as.Date(paste0(month_label, "-01")), .groups="drop")
    }, error=function(e) NULL)
  }) |> compact()
  
  if(length(ts_list)==0) return(NULL)
  bind_rows(ts_list) |>
    mutate(gw = -gw) |>
    dplyr::select(well_num, date, gw, month_label, source_id)
}

## 4.4 BC Real-Time (ArcGIS FeatureServer) ----------------------------------
load_bc_realtime <- function() {
  if(!cfg$src_bc_realtime) return(NULL)
  message("Loading BC Real-Time groundwater (ArcGIS) …")
  basin_wkt <- sf::st_as_text(sf::st_transform(basin, 3857))
  query_url <- paste0(cfg$bc_rt_layer, "/query?where=1%3D1&outFields=*&f=json&geometry=",
                      URLencode(basin_wkt), "&geometryType=esriGeometryPolygon&inSR=3857")
  
  resp <- GET(query_url, timeout(60))
  if(status_code(resp)!=200) { message("  ArcGIS query failed"); return(NULL) }
  
  dat <- fromJSON(content(resp, "text", encoding="UTF-8"))$features
  if(is.null(dat) || length(dat)==0) return(NULL)
  
  parse_feat <- function(f) {
    a <- f$attributes
    tibble(
      well_num = as.character(a[[grep("well|id|station", names(a), ignore.case=TRUE, value=TRUE)[1]]]),
      date = as.Date(a[[grep("date|time|obs", names(a), ignore.case=TRUE, value=TRUE)[1]]]),
      gw_raw = as.numeric(a[[grep("level|depth|value|gwl", names(a), ignore.case=TRUE, value=TRUE)[1]]]),
      lat = as.numeric(a[[grep("y|lat", names(a), ignore.case=TRUE, value=TRUE)[1]]]),
      lon = as.numeric(a[[grep("x|lon", names(a), ignore.case=TRUE, value=TRUE)[1]]])
    )
  }
  df <- map_dfr(dat, parse_feat) |>
    standardize_ts("BC_RT", "well_num", "date", "gw_raw", negate=TRUE) |>
    st_as_sf(coords=c("lon","lat"), crs=4326, remove=FALSE) |>
    filter(st_within(., basin, sparse=FALSE)[,1]) |>
    as_tibble() |>
    dplyr::select(well_num, date, gw, month_label, source_id)
  
  message(sprintf("  BC Real-Time: %d records", nrow(df)))
  df
}

## 4.5 GWELLS API (static + sparse TS) --------------------------------------
load_gwells <- function() {
  if(!cfg$src_gwells) return(NULL)
  message("Loading GWELLS API (BC) …")
  basin_bbox <- st_bbox(st_transform(basin, 3857))
  api_url <- paste0(cfg$gwells_api_base, "/wells/?min_longitude=", basin_bbox$xmin,
                    "&max_longitude=", basin_bbox$xmax,
                    "&min_latitude=", basin_bbox$ymin,
                    "&max_latitude=", basin_bbox$ymax,
                    "&has_water_level=true&limit=1000")
  
  resp <- GET(api_url, add_headers("Accept"="application/json"), timeout(60))
  if(status_code(resp)!=200) return(NULL)
  
  dat <- fromJSON(content(resp, "text"))$results
  if(length(dat)==0) return(NULL)
  
  ts_list <- lapply(dat, function(w) {
    wl <- w$water_levels
    if(is.null(wl) || length(wl) < cfg$gwells_min_obs) return(NULL)
    tibble(
      well_num = as.character(w$well_tag_number),
      date = as.Date(sapply(wl, `[[`, "observation_date")),
      gw_raw = as.numeric(sapply(wl, `[[`, "water_level")),
      lat = w$latitude, lon = w$longitude
    )
  }) |> compact()
  
  if(length(ts_list)==0) return(NULL)
  bind_rows(ts_list) |>
    standardize_ts("GWELLS", "well_num", "date", "gw_raw", negate=TRUE) |>
    st_as_sf(coords=c("lon","lat"), crs=4326, remove=FALSE) |>
    filter(st_within(., basin, sparse=FALSE)[,1]) |>
    as_tibble() |>
    dplyr::select(well_num, date, gw, month_label, source_id)
}

#---- 5. Aggregate All Sources -----------------------------------------------
message("\n=== Loading groundwater time series ===")
sources <- list(
  PGOWN = load_pgown(),
  GIN   = load_gin(),
  WSC   = load_wsc(),
  BC_RT = load_bc_realtime(),
  GWELLS= load_gwells()
)

# Combine non-NULL sources
ts_all <- bind_rows(compact(sources))
if(nrow(ts_all)==0) stop("No data loaded from any source.")

# Deduplicate: same well_num + month_label → keep source priority order
source_priority <- c("PGOWN", "GIN", "WSC", "BC_RT", "GWELLS")
ts_dedup <- ts_all |>
  mutate(src_rank = match(source_id, source_priority)) |>
  group_by(well_num, month_label) |>
  slice_min(src_rank, with_ties=FALSE) |>
  ungroup() |> select(-src_rank)

message(sprintf("Combined: %d unique wells, %d monthly records (after dedup)",
                n_distinct(ts_dedup$well_num), nrow(ts_dedup)))

#---- 6. Quality Filters (Source-Adaptive) -----------------------------------
well_stats <- ts_dedup |>
  group_by(well_num) |>
  summarise(
    n_months = n(),
    date_min = min(date), date_max = max(date),
    record_years = as.numeric(difftime(max(date), min(date), units="days"))/365.25,
    n_expected = as.numeric(difftime(max(date), min(date), units="days"))/30.44 + 1,
    gap_frac = 1 - n_months / pmax(1, n_expected),
    sources_used = paste(unique(source_id), collapse=","),
    .groups="drop"
  )

wells_ok <- well_stats |>
  filter(
    record_years >= if(grepl("PGOWN", sources_used)) cfg$min_years else pmax(5, cfg$min_years*0.5),
    gap_frac <= cfg$max_gap_frac * if(grepl("PGOWN", sources_used)) 1 else 1.5
  )

message(sprintf("Wells after quality filters: %d", nrow(wells_ok)))
if(nrow(wells_ok)==0) stop("No wells pass filters. Relax cfg$min_years or cfg$max_gap_frac.")

ts_filtered <- ts_dedup |> filter(well_num %in% wells_ok$well_num)

#---- 7. Compute SGI ---------------------------------------------------------
message("Computing SGI (Bloomfield & Marchant 2013) …")
sgi_long <- ts_filtered |>
  mutate(cal_month = month(as.Date(paste0(month_label, "-01")))) |>
  group_by(well_num, cal_month) |>
  mutate(sgi = normal_scores(gw)) |>
  ungroup()

# Basin aggregate: median across wells
sgi_basin <- sgi_long |>
  group_by(month_label) |>
  summarise(
    sgi_median = median(sgi, na.rm=TRUE),
    sgi_mean = mean(sgi, na.rm=TRUE),
    sgi_sd = sd(sgi, na.rm=TRUE),
    n_wells = sum(!is.na(sgi)),
    sources_contrib = paste(sort(unique(source_id)), collapse=","),
    .groups="drop"
  ) |>
  arrange(month_label) |>
  mutate(
    date = as.Date(paste0(month_label, "-01")),
    drought_cat = case_when(
      sgi_median >  0.00 ~ "No drought",
      sgi_median >= -1.00 ~ "Minor drought",
      sgi_median >= -1.50 ~ "Moderate drought", 
      sgi_median >= -2.00 ~ "Severe drought",
      TRUE ~ "Extreme drought"
    )
  )

if(!is.null(cfg$target_months)) sgi_basin <- sgi_basin |> filter(month_label %in% cfg$target_months)

#---- 8. Output --------------------------------------------------------------
cat("\n========================================================\n")
cat("  SGI – Nechako River Basin (multi-source)\n")
cat(sprintf("  Sources active: %s\n", paste(names(compact(sources)), collapse=", ")))
cat(sprintf("  Wells: %d | Records: %d | Date range: %s to %s\n",
            n_distinct(sgi_long$well_num), nrow(sgi_long),
            min(sgi_long$date), max(sgi_long$date)))
cat("========================================================\n\n")
print(sgi_basin |> select(month_label, sgi_median, n_wells, sources_contrib, drought_cat), n=30)

# Save
out_csv <- paste0(cfg$output_prefix, ".csv")
write_csv(sgi_basin, out_csv)
message(sprintf("\nSaved: %s", out_csv))

#---- 9. Plot ----------------------------------------------------------------
p <- ggplot(sgi_basin, aes(x=date, y=sgi_median)) +
  geom_ribbon(aes(ymin=pmin(sgi_median,0), ymax=0), fill="#d73027", alpha=0.3) +
  geom_ribbon(aes(ymin=0, ymax=pmax(sgi_median,0)), fill="#4575b4", alpha=0.3) +
  geom_line(linewidth=0.5, colour="grey20") +
  geom_hline(yintercept=c(0,-1,-1.5,-2), linetype=c("solid","dashed","dashed","dashed"),
             colour=c("black","#fc8d59","#d73027","#a50026"), linewidth=0.3) +
  scale_x_date(date_breaks="2 years", date_labels="%Y") +
  scale_y_continuous(breaks=c(-2,-1.5,-1,0,1,2),
                     labels=c("-2\n(Extreme)","-1.5\n(Severe)","-1\n(Moderate)","0","1","2")) +
  labs(
    title="Standardised Groundwater Index (SGI) – Nechako Basin",
    subtitle=paste0("Multi-source | ", n_distinct(sgi_long$well_num), " wells | ",
                    paste(names(compact(sources)), collapse="+")),
    x=NULL, y="SGI",
    caption="Method: Bloomfield & Marchant (2013) | Data: PGOWN+GIN+WSC+BC_RT+GWELLS"
  ) +
  theme_bw(base_size=10) +
  theme(plot.title=element_text(face="bold"), panel.grid.minor=element_blank())

ggsave(paste0(cfg$output_prefix, ".png"), p, width=10, height=4, dpi=150)
message("Plot saved: ", cfg$output_prefix, ".png")

#---- END --------------------------------------------------------------------