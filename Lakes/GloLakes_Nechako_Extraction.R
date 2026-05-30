# ============================================================================
# GloLakes + GloRiver — NECHAKO BASIN EXTRACTION
# Based on confirmed NetCDF structure (files inspected February 2026)
#
# Dataset: GloLakes v2 / Global Water Monitor (ANU, Dr. Jiawei Hou)
# Reference: Hou et al. (2024), Earth Syst. Sci. Data, 16, 201–218
#            https://doi.org/10.5194/essd-16-201-2024  
# Data: https://doi.org/10.25914/K8ZF-6G46  
#
# ── CONFIRMED FILE STRUCTURE ─────────────────────────────────────────────────
#
# THREE LAKE FILES (storage in MCM = million cubic metres = Mm³):
#
#   (1) LandsatPlusSentinel2  — 4,054 lakes × 1,572 time steps
#       time encoding : "days since 1984-01-01"  (actual day offsets, likely
#                        irregular Landsat/S2 observation dates → resampled
#                        to monthly means in this script)
#       lake_name     : MOSTLY EMPTY — must use lat/lon bounding box
#       data_quality  : dims=[time] — global quality flag per time step,
#                        NOT per lake (Q1/Q2/Q3/Q4 describes the epoch/method)
#
#   (2) LandsatPlusGREALM   — 135 lakes × 506 months (1984-01 to 2026-02)
#       time encoding : "monthly since 1984-01-01"
#       lake_name     : ALL NAMED (Atitla, Junin, Naivasha ...)
#       Coverage      : large radar-altimetry lakes only; check if BC lakes included
#
#   (3) LandsatPlusICESat2  — 22,232 lakes × 506 months (1984-01 to 2026-02)
#       time encoding : "monthly since 1984-01-01"
#       lake_name     : MANY EMPTY — must use lat/lon bounding box
#       Coverage      : largest set; best chance of Nechako lakes
#
# UNITS CONFIRMED: MCM = million cubic metres = Mm³  → NO CONVERSION NEEDED
#   (SWSI_Garen1993_v2.R uses Mm³ throughout — direct plug-in)
#
# QUALITY FLAG (data_quality, dims=[time]):
#   Q1 — absolute volume via geostatistical model + satellite lake extents
#   Q2 — absolute volume via V-H (volume–height) relationship
#   Q3 — RELATIVE volume from satellite heights + extents  ← NOT absolute
#   Q4 — RELATIVE volume from heights + A-H derived extents ← NOT absolute
#   IMPORTANT: Q3/Q4 months give only anomalies, not total volumes.
#   They cannot be used directly in SWSI formula which needs absolute storage.
#   → This script flags Q3/Q4 months and excludes them from the SWSI output.
#
# TWO RIVER DISCHARGE FILES (see bottom of script — Section 9 — for details
#   on whether these are useful for your research):
#   MSGR: 58,264 stations, 8-day composites, 2000-present (satellite-only)
#   GSGR:  9,366 stations, 8-day composites, 2000-present (gauge-satellite fused)
#
# ── NECHAKO BASIN LAKE TARGETS ────────────────────────────────────────────────
#   Stuart Lake         ~825 km², 54.9°N, 124.1°W  — natural, unregulated
#   Nechako Reservoir   ~910 km², 53.5°N, 126.8°W  — regulated (Rio Tinto)
#   Fraser Lake         ~176 km², 54.0°N, 124.8°W  — natural
#   Babine Lake         ~479 km², 54.8°N, 126.0°W  — natural (check attribution)
#   Ootsa/Eutsuk/Whitesail — part of Nechako Reservoir complex
# ============================================================================

rm(list = ls())

# ── Packages ──────────────────────────────────────────────────────────────────
for (pkg in c("ncdf4", "tidyverse", "lubridate", "sf")) {
  if (!require(pkg, character.only=TRUE, quietly=TRUE)) {
    install.packages(pkg, repos="https://cran.rstudio.com/")
    library(pkg, character.only=TRUE)
  }
}


# ============================================================================
# SECTION 1: FILE PATHS AND SETTINGS
# ============================================================================

# Folder containing all five .nc files
# MODIFIED: Changed to user-specified path
INPUT_DIR  <- "D:/Nechako_Drought/Nechako/Lakes/Glolakes"

# ── Lake storage files ───────────────────────────────────────────────────────
FILE_S2     <- file.path(INPUT_DIR,
                         "Global_Lake_Absolute_Storage_LandsatPlusSentinel2 (1984-present).nc")
FILE_GREALM <- file.path(INPUT_DIR,
                         "Global_Lake_Absolute_Storage_LandsatPlusGREALM (1984-present).nc")
FILE_ICE    <- file.path(INPUT_DIR,
                         "Global_Lake_Absolute_Storage_LandsatPlusICESat2 (1984-present).nc")

# ── River discharge files ────────────────────────────────────────────────────
FILE_MSGR   <- file.path(INPUT_DIR, "Global_River_Discharge_MSGR (2000-present).nc")
FILE_GSGR   <- file.path(INPUT_DIR, "Global_River_Discharge_GSGR (2000-present).nc")

# ── Output directory ─────────────────────────────────────────────────────────
# Now automatically creates subfolder under the new path
OUTPUT_DIR  <- file.path(INPUT_DIR, "nechako_extracted")
dir.create(OUTPUT_DIR, showWarnings=FALSE, recursive=TRUE)

# ── MODIFIED: Nechako basin boundary shapefile (replaces simple bounding box) ─
NECHAKO_BOUNDARY <- "D:/Nechako_Drought/Nechako/Spatial/nechakoBound_dissolve.shp"

# ── Fallback bounding box (used for quick pre-filter before shapefile test) ─
# Derived from shapefile extent with small buffer for efficiency
BASIN_LAT_MIN  <- 52.0
BASIN_LAT_MAX  <- 56.5
BASIN_LON_MIN  <- -128.5
BASIN_LON_MAX  <- -119.5

# ── Name keywords (cross-check for lakes with centroids near bounding box edge)
NAME_KEYWORDS  <- c(
  "stuart", "nechako", "fraser", "ootsa", "eutsuk", "whitesail",
  "tahtsa", "natalkuz", "knewstubb", "babine", "tachie",
  "nadsilnich", "tchesinkut", "cheslatta", "endako", "manson"
)

# ── SWSI water year start ────────────────────────────────────────────────────
WY_START_MONTH <- 10L   # October

# ── Analysis date for reporting ─────────────────────────────────────────────
ANALYSIS_DATE <- format(Sys.Date(), "%Y-%m-%d")


# ============================================================================
# SECTION 1B: LOAD NECHAKO BASIN BOUNDARY SHAPEFILE
# ============================================================================
cat("\n", strrep("─",60), "\n")
cat("SECTION 1B: LOADING NECHAKO BASIN BOUNDARY SHAPEFILE\n")
cat(strrep("─",60), "\n")

if (file.exists(NECHAKO_BOUNDARY)) {
  # Load the shapefile
  basin_shp <- sf::st_read(NECHAKO_BOUNDARY, quiet=TRUE)
  
  # Ensure it's in WGS84 (EPSG:4326) for lat/lon matching
  if (!sf::st_is_longlat(basin_shp)) {
    basin_shp <- sf::st_transform(basin_shp, 4326)
    cat("  Transformed shapefile to WGS84 (EPSG:4326)\n")
  }
  
  # Extract bounding box from shapefile for quick pre-filter
  bbox_shp <- sf::st_bbox(basin_shp)
  BASIN_LON_MIN <- bbox_shp["xmin"]
  BASIN_LON_MAX <- bbox_shp["xmax"]
  BASIN_LAT_MIN <- bbox_shp["ymin"]
  BASIN_LAT_MAX <- bbox_shp["ymax"]
  
  cat(sprintf("  Shapefile loaded successfully: %s\n", NECHAKO_BOUNDARY))
  cat(sprintf("  Shapefile extent: %.2f to %.2f°N, %.2f to %.2f°W\n",
              BASIN_LAT_MIN, BASIN_LAT_MAX, BASIN_LON_MIN, BASIN_LON_MAX))
  cat(sprintf("  Number of polygons: %d\n", nrow(basin_shp)))
  
  # Store basin_shp in global environment for use in extraction function
  assign("BASIN_SHAPEFILE", basin_shp, envir=.GlobalEnv)
  
} else {
  cat(sprintf("  WARNING: Shapefile not found: %s\n", NECHAKO_BOUNDARY))
  cat("  Falling back to simple bounding box method only.\n")
  assign("BASIN_SHAPEFILE", NULL, envir=.GlobalEnv)
}


# ============================================================================
# SECTION 2: TIME DECODING UTILITIES
#
# The three lake files use TWO different time encodings:
#   S2 file    : "days since 1984-01-01"  — time values are actual day
#                offsets (e.g. 0, 16, 31, 47 ...) and may be IRREGULAR.
#                These represent individual Landsat/Sentinel-2 observation
#                dates.  This script resamples them to monthly medians.
#   GREALM/ICE : "monthly since 1984-01-01" — time values are sequential
#                integers 0, 1, 2 ... = month indices.
# ============================================================================

# NULL-coalescing helper
`%||%` <- function(a, b) if (!is.null(a)) a else b

decode_time <- function(nc, time_var="time") {
  raw       <- as.numeric(ncdf4::ncvar_get(nc, time_var))
  atts      <- ncdf4::ncatt_get(nc, time_var)
  units_str <- tolower(if (!is.null(atts$units)) atts$units else "")
  
  epoch_1984 <- as.Date("1984-01-01")
  
  if (grepl("days since 1984", units_str) ||
      grepl("days since 1984", tolower(atts$time_unit %||% ""))) {
    dates <- epoch_1984 + as.integer(round(raw))
    cat(sprintf("    Time: day offsets since 1984-01-01 | range %s to %s\n",
                min(dates), max(dates)))
    is_irregular <- (length(unique(diff(as.integer(sort(dates))))) > 3)
    if (is_irregular)
      cat("    NOTE: irregular observation dates — will resample to monthly\n")
    return(list(dates=dates, is_monthly=FALSE, is_irregular=is_irregular))
    
  } else if (grepl("monthly since 1984", units_str) ||
             grepl("monthly", units_str)) {
    dates <- epoch_1984 %m+% months(as.integer(raw))
    cat(sprintf("    Time: month index since 1984-01-01 | range %s to %s\n",
                min(dates), max(dates)))
    return(list(dates=dates, is_monthly=TRUE, is_irregular=FALSE))
    
  } else {
    cat(sprintf("    Time units not recognised: '%s'\n  Trying month-index fallback\n",
                units_str))
    dates <- epoch_1984 %m+% months(as.integer(raw))
    cat(sprintf("    Fallback dates: %s to %s\n", min(dates), max(dates)))
    return(list(dates=dates, is_monthly=TRUE, is_irregular=FALSE))
  }
}

# Robust extraction of string variables from NetCDF
extract_string_var <- function(nc, var_name, n_expected) {
  if (!(var_name %in% names(nc$var))) {
    return(rep(NA_character_, n_expected))
  }
  v <- ncdf4::ncvar_get(nc, var_name)
  
  # Case 1: 1D array or vector of strings (NetCDF4 NC_STRING type)
  if ((is.null(dim(v)) || length(dim(v)) == 1) && length(v) == n_expected) {
    return(trimws(as.character(v)))
  }
  
  # Case 2: 2D array where second dimension matches n_expected (C-order strings)
  if (length(dim(v)) == 2 && dim(v)[2] == n_expected) {
    if (is.character(v)) {
      # Character matrix: each cell is a single character (or empty)
      return(apply(v, 2, function(col) trimws(paste(col, collapse = ""))))
    } else if (is.raw(v) || is.integer(v)) {
      # Raw or integer matrix: convert each column from bytes to string
      return(apply(v, 2, function(col) {
        bytes <- as.raw(col[col > 0 & col < 128])   # keep only printable ASCII
        if (length(bytes) == 0) "" else rawToChar(bytes)
      }))
    }
  }
  
  # Case 3: 2D array where first dimension matches n_expected (Fortran-order strings)
  if (length(dim(v)) == 2 && dim(v)[1] == n_expected) {
    if (is.character(v)) {
      return(apply(v, 1, function(row) trimws(paste(row, collapse = ""))))
    }
  }
  
  # Updated warning to print the actual dimensions for easier debugging
  warning(sprintf("Unexpected format for variable '%s' (dims: %s) – returning NAs", 
                  var_name, paste(dim(v), collapse="x")))
  rep(NA_character_, n_expected)
}

# ── NEW: Extract optional numeric metadata (safe, non-breaking) ───────────
extract_optional_numeric_var <- function(nc, var_name, n_expected) {
  if (!(var_name %in% names(nc$var))) return(rep(NA_real_, n_expected))
  tryCatch({
    v <- as.numeric(ncdf4::ncvar_get(nc, var_name))
    if (length(v) == n_expected) return(v)
    else return(rep(NA_real_, n_expected))
  }, error = function(e) rep(NA_real_, n_expected))
}

# ── NEW: Point-in-polygon test using shapefile boundary ────────────────────
check_point_in_basin <- function(lat, lon, basin_shp) {
  if (is.null(basin_shp)) {
    # Fallback to bounding box if no shapefile
    return(lon >= BASIN_LON_MIN & lon <= BASIN_LON_MAX &
             lat >= BASIN_LAT_MIN & lat <= BASIN_LAT_MAX)
  }
  
  # Create sf points from coordinates
  points_sf <- sf::st_as_sf(data.frame(lon=lon, lat=lat),
                            coords=c("lon", "lat"),
                            crs=4326)
  
  # Test which points fall within the basin polygon
  in_basin <- sf::st_intersects(points_sf, basin_shp, sparse=FALSE)
  
  # Return logical vector (TRUE = inside basin)
  return(as.logical(in_basin[,1]))
}


# ============================================================================
# SECTION 3: GENERIC LAKE EXTRACTION FUNCTION
#
# Searches one NetCDF file for Nechako lakes using:
#   (a) Shapefile boundary (PRIMARY) — point-in-polygon test
#   (b) Bounding box (pre-filter for efficiency)
#   (c) lake_name keyword matching (handles empty-name lakes gracefully)
#
# Extracts the full lake_storage time series for matched lakes.
# Resamples to monthly if the S2 file's irregular observation dates are found.
# Returns a tidy data frame and a basin-total monthly series.
# ============================================================================

extract_nechako_lakes <- function(nc_path, product_label) {
  
  cat(sprintf("\n%s\nEXTRACTING: %s\n%s\n",
              strrep("─",60), product_label, strrep("─",60)))
  
  if (!file.exists(nc_path)) {
    cat(sprintf("  File not found: %s\n  Skipping.\n", nc_path))
    return(NULL)
  }
  
  nc <- ncdf4::nc_open(nc_path)
  on.exit(ncdf4::nc_close(nc))
  
  n_lakes <- nc$dim$ID$len
  cat(sprintf("  Lakes in file: %d\n", n_lakes))
  
  # Coordinates
  lat_all  <- as.numeric(ncdf4::ncvar_get(nc, "latitude"))
  lon_all  <- as.numeric(ncdf4::ncvar_get(nc, "longitude"))
  
  # Lake names (robust extraction)
  lake_names <- extract_string_var(nc, "lake_name", n_lakes)
  
  # Country, state, basin names
  country_names <- extract_string_var(nc, "country_name", n_lakes)
  state_names   <- extract_string_var(nc, "state_name", n_lakes)
  basin_names   <- extract_string_var(nc, "basin_name", n_lakes)
  
  # ── Optional lake metadata (extract if available) ───────────────────────
  lake_area_km2    <- extract_optional_numeric_var(nc, "lake_area_km2", n_lakes)
  lake_elevation_m <- extract_optional_numeric_var(nc, "elevation", n_lakes)
  lake_id          <- extract_string_var(nc, "lake_id", n_lakes)
  
  # ── Time ──────────────────────────────────────────────────────────────────
  time_info <- decode_time(nc)
  dates     <- time_info$dates
  n_time    <- length(dates)
  
  # ── Quality flag ─────────────────────────────────────────────────────────
  # data_quality has dims=[time] — one value per time step for ALL lakes
  qual_raw   <- as.character(ncdf4::ncvar_get(nc, "data_quality"))
  # Map q1/q2/q3/q4 strings to integers
  qual_int   <- dplyr::case_when(
    tolower(qual_raw) == "q1" ~ 1L,
    tolower(qual_raw) == "q2" ~ 2L,
    tolower(qual_raw) == "q3" ~ 3L,
    tolower(qual_raw) == "q4" ~ 4L,
    TRUE                      ~ NA_integer_
  )
  # Months where quality is absolute (Q1 or Q2 — usable in SWSI)
  is_absolute_month <- qual_int %in% c(1L, 2L)
  cat(sprintf("  Quality: Q1=%d, Q2=%d, Q3=%d, Q4=%d time steps\n",
              sum(qual_int==1L,na.rm=T), sum(qual_int==2L,na.rm=T),
              sum(qual_int==3L,na.rm=T), sum(qual_int==4L,na.rm=T)))
  cat(sprintf("  Absolute (Q1+Q2) months: %d / %d\n",
              sum(is_absolute_month), n_time))
  
  # ── SEARCH 0: Quick bounding box pre-filter (for efficiency) ─────────────
  bbox_idx <- which(
    lat_all >= BASIN_LAT_MIN & lat_all <= BASIN_LAT_MAX &
      lon_all >= BASIN_LON_MIN & lon_all <= BASIN_LON_MAX &
      !is.na(lat_all) & !is.na(lon_all)
  )
  cat(sprintf("  Bounding box pre-filter hits: %d\n", length(bbox_idx)))
  
  # ── SEARCH 1: Shapefile boundary (PRIMARY METHOD) ────────────────────────
  if (!is.null(get0("BASIN_SHAPEFILE", envir=.GlobalEnv))) {
    basin_shp <- get("BASIN_SHAPEFILE", envir=.GlobalEnv)
    cat("  Testing points against shapefile boundary...\n")
    
    # Only test lakes that passed the bounding box pre-filter
    in_basin <- rep(FALSE, n_lakes)
    if (length(bbox_idx) > 0) {
      in_basin[bbox_idx] <- check_point_in_basin(lat_all[bbox_idx], 
                                                 lon_all[bbox_idx], 
                                                 basin_shp)
    }
    shapefile_idx <- which(in_basin)
    cat(sprintf("  Shapefile boundary hits: %d\n", length(shapefile_idx)))
  } else {
    shapefile_idx <- bbox_idx
    cat("  Using bounding box only (no shapefile available)\n")
  }
  
  # ── SEARCH 2: Name keywords ────────────────────────────────────────────────
  # Extra safety: also search country=Canada to avoid false positives
  pattern    <- paste(NAME_KEYWORDS, collapse="|")
  name_idx   <- grep(pattern, lake_names, ignore.case=TRUE)
  canada_idx <- grep("Canada", country_names, ignore.case=TRUE)
  # Name hits that are in Canada (avoids e.g. "Fraser" in Australia)
  name_canada_idx <- intersect(name_idx, canada_idx)
  cat(sprintf("  Name+Canada hits: %d\n", length(name_canada_idx)))
  
  # All candidate indices (union of shapefile + name matching)
  all_idx <- sort(unique(c(shapefile_idx, name_canada_idx)))
  cat(sprintf("  Combined candidates: %d\n", length(all_idx)))
  
  if (length(all_idx) == 0) {
    cat("  *** NONE FOUND — no Nechako lakes in this product ***\n")
    return(NULL)
  }
  
  # ── Build metadata table ────────────────────────────────────────────────
  meta <- data.frame(
    file_idx         = all_idx,               # 1-based R index in the file
    product          = product_label,
    lake_name        = lake_names[all_idx],
    lake_id          = lake_id[all_idx],              # NEW
    lat              = lat_all[all_idx],
    lon              = lon_all[all_idx],
    lake_area_km2    = lake_area_km2[all_idx],        # NEW
    elevation_m      = lake_elevation_m[all_idx],     # NEW
    country          = country_names[all_idx],
    state            = state_names[all_idx],
    basin            = basin_names[all_idx],
    found_by      = dplyr::case_when(
      all_idx %in% shapefile_idx & all_idx %in% name_canada_idx ~ "shapefile+name",
      all_idx %in% shapefile_idx                                ~ "shapefile",
      all_idx %in% bbox_idx                                     ~ "bbox",
      TRUE                                                       ~ "name+Canada"
    ),
    stringsAsFactors = FALSE
  )
  
  cat("\n  CANDIDATE LAKES (review before trusting):\n")
  # Replace NA with empty string for cleaner display
  meta_display <- meta[, c("file_idx","lake_name","lat","lon","lake_area_km2","state","basin","found_by")]
  meta_display[] <- lapply(meta_display, function(x) ifelse(is.na(x), "", x))
  print(meta_display, row.names=FALSE)
  
  # ── Extract lake_storage [ID × time] for matched indices ──────────────────
  n_match <- nrow(meta)
  storage_mat <- matrix(NA_real_, nrow=n_match, ncol=n_time)
  
  for (k in seq_len(n_match)) {
    idx_k <- meta$file_idx[k]
    # Read one lake at a time (start=[lake_idx,1], count=[1, n_time])
    raw <- as.numeric(ncdf4::ncvar_get(nc, "lake_storage",
                                       start=c(idx_k, 1),
                                       count=c(1, n_time)))
    raw[is.nan(raw) | is.na(raw)] <- NA_real_
    storage_mat[k, ] <- raw
    cat(sprintf("    [%d/%d] %s (idx %d): %d non-NA values, mean=%.1f Mm³\n",
                k, n_match,
                if (!is.na(meta$lake_name[k]) && nzchar(meta$lake_name[k])) meta$lake_name[k] else "(unnamed)",
                idx_k,
                sum(!is.na(raw)),
                mean(raw, na.rm=TRUE)))
  }
  
  # ── For S2 file: resample irregular dates to monthly medians ─────────────
  if (!time_info$is_monthly) {
    cat("\n  Resampling irregular observation dates to calendar months...\n")
    yr_mo     <- format(dates, "%Y-%m")
    months_u  <- sort(unique(yr_mo))
    monthly_dates   <- as.Date(paste0(months_u, "-01"))
    monthly_storage <- matrix(NA_real_, nrow=n_match, ncol=length(months_u))
    monthly_qual    <- integer(length(months_u))
    
    for (j in seq_along(months_u)) {
      t_idx <- which(yr_mo == months_u[j])
      # Storage: median of all observations within the month
      for (k in seq_len(n_match)) {
        vals <- storage_mat[k, t_idx]
        monthly_storage[k, j] <- if (all(is.na(vals))) NA_real_
        else median(vals, na.rm=TRUE)
      }
      # Quality: worst (highest) quality flag in the month
      monthly_qual[j] <- if (all(is.na(qual_int[t_idx]))) NA_integer_
      else max(qual_int[t_idx], na.rm=TRUE)
    }
    
    # Replace originals with monthly resampled versions
    dates           <- monthly_dates
    storage_mat     <- monthly_storage
    qual_int        <- monthly_qual
    is_absolute_month <- qual_int %in% c(1L, 2L)
    n_time          <- length(dates)
    cat(sprintf("  Resampled to %d monthly values (%s to %s)\n",
                n_time, min(dates), max(dates)))
  }
  
  # ── Build long-format data frame ─────────────────────────────────────────
  long_rows <- list()
  for (k in seq_len(n_match)) {
    df_k <- data.frame(
      date             = dates,
      year             = year(dates),
      month            = month(dates),
      water_year       = ifelse(month(dates) >= WY_START_MONTH,
                                year(dates), year(dates) - 1L),
      storage_Mm3      = storage_mat[k, ],
      quality          = qual_int,
      is_absolute      = is_absolute_month,  # TRUE = Q1/Q2 = usable in SWSI
      file_idx         = meta$file_idx[k],
      lake_name        = meta$lake_name[k],
      lake_id          = meta$lake_id[k],              # NEW
      lat              = meta$lat[k],
      lon              = meta$lon[k],
      lake_area_km2    = meta$lake_area_km2[k],        # NEW
      elevation_m      = meta$elevation_m[k],          # NEW
      country          = meta$country[k],
      state            = meta$state[k],
      basin            = meta$basin[k],
      product          = product_label,
      stringsAsFactors = FALSE
    )
    long_rows[[k]] <- df_k
  }
  long_df <- bind_rows(long_rows)
  
  # ── Basin-total monthly storage (Q1+Q2 only) ─────────────────────────────
  # Sum all matched lakes; require all lakes to have a value in that month.
  # Months where ANY lake is NA or Q3/Q4 are flagged but included as partial.
  basin_total <- long_df %>%
    group_by(date, year, month, water_year, quality, is_absolute) %>%
    summarise(
      total_storage_Mm3    = sum(storage_Mm3, na.rm=TRUE),
      n_lakes_contributing = sum(!is.na(storage_Mm3)),
      n_lakes_total        = n_match,
      all_lakes_present    = sum(!is.na(storage_Mm3)) == n_match,
      .groups = "drop"
    ) %>%
    mutate(
      # Set total to NA if no lakes contributed or quality is relative-only
      total_storage_Mm3 = case_when(
        n_lakes_contributing == 0 ~ NA_real_,
        !is_absolute              ~ NA_real_,  # Q3/Q4: exclude from SWSI
        TRUE                      ~ total_storage_Mm3
      ),
      product = product_label
    )
  
  n_usable <- sum(!is.na(basin_total$total_storage_Mm3))
  cat(sprintf("\n  Basin total: %d usable months (Q1+Q2), %d Q3/Q4 excluded\n",
              n_usable,
              sum(!basin_total$is_absolute)))
  if (n_usable > 0) {
    cat(sprintf("  Storage range: %.0f to %.0f Mm³  (mean %.0f Mm³)\n",
                min(basin_total$total_storage_Mm3, na.rm=TRUE),
                max(basin_total$total_storage_Mm3, na.rm=TRUE),
                mean(basin_total$total_storage_Mm3, na.rm=TRUE)))
  }
  
  return(list(
    meta        = meta,
    long_df     = long_df,
    basin_total = basin_total,
    dates       = dates,
    storage_mat = storage_mat,
    qual_int    = qual_int,
    n_lakes     = n_match,
    time_range  = c(min(dates), max(dates)),
    product     = product_label
  ))
}


# ============================================================================
# SECTION 4: RUN EXTRACTION ON ALL THREE LAKE FILES
# ============================================================================
cat("\n", strrep("═",65), "\n")
cat("GLOLAKES NECHAKO EXTRACTION — All three lake products\n")
cat(strrep("═",65), "\n")

result_s2     <- extract_nechako_lakes(FILE_S2,     "LandsatPlusSentinel2")
result_grealm <- extract_nechako_lakes(FILE_GREALM, "LandsatPlusGREALM")
result_ice    <- extract_nechako_lakes(FILE_ICE,    "LandsatPlusICESat2")

# Collect non-NULL results
results_all <- Filter(Negate(is.null), list(result_s2, result_grealm, result_ice))
cat(sprintf("\nProducts with Nechako lakes found: %d / 3\n", length(results_all)))


# ============================================================================
# SECTION 5: SAVE PER-PRODUCT OUTPUTS
# ============================================================================
for (res in results_all) {
  lbl   <- gsub("[^A-Za-z0-9]", "_", res$product)
  
  # Long CSV (with blanks instead of "NA")
  write.csv(res$long_df,
            file.path(OUTPUT_DIR, paste0(lbl, "_long.csv")),
            row.names=FALSE, na = "")
  
  # Wide CSV (one column per lake, labelled by lat/lon if unnamed)
  wide_df <- data.frame(date = res$dates,
                        year = year(res$dates),
                        month = month(res$dates),
                        water_year = ifelse(month(res$dates) >= WY_START_MONTH,
                                            year(res$dates), year(res$dates) - 1L),
                        quality = res$qual_int,
                        is_absolute = res$qual_int %in% c(1L, 2L))
  
  for (k in seq_len(res$n_lakes)) {
    # Build a clean lake label
    if (!is.na(res$meta$lake_name[k]) && nzchar(res$meta$lake_name[k])) {
      nm <- gsub(" ", "_", res$meta$lake_name[k])
    } else {
      # Use coordinates (ensure longitude sign is consistent with metadata)
      nm <- sprintf("lake_%.4fN_%.4fW", res$meta$lat[k], abs(res$meta$lon[k]))
    }
    # Final safety: if the label is somehow still NA, fall back to an index
    if (is.na(nm)) nm <- paste0("lake_", k)
    
    wide_df[[nm]] <- res$storage_mat[k, ]
  }
  write.csv(wide_df,
            file.path(OUTPUT_DIR, paste0(lbl, "_wide.csv")),
            row.names=FALSE, na = "")
  
  # Basin total CSV
  write.csv(res$basin_total,
            file.path(OUTPUT_DIR, paste0(lbl, "_basin_total_Mm3.csv")),
            row.names=FALSE, na = "")
  
  # RDS
  saveRDS(res, file.path(OUTPUT_DIR, paste0(lbl, ".rds")))
  cat(sprintf("  Saved: %s  (%d lakes, %d months)\n",
              lbl, res$n_lakes, nrow(res$basin_total)))
}


# ============================================================================
# SECTION 5B: CONSOLIDATED LAKE SUMMARY REPORT  [NEW - Minimal addition]
# ============================================================================
cat("\n", strrep("─",60), "\n")
cat("SECTION 5B: CREATING CONSOLIDATED LAKE SUMMARY\n")
cat(strrep("─",60), "\n")

# Combine metadata from all products
all_lake_meta <- bind_rows(lapply(results_all, function(res) res$meta))

# Add derived reporting fields
all_lake_meta <- all_lake_meta %>%
  mutate(
    lake_uid = sprintf("%.4fN_%.4fW", lat, abs(lon)),
    is_target_lake = grepl(paste(NAME_KEYWORDS, collapse="|"), 
                           lake_name, ignore.case=TRUE)
  ) %>%
  select(file_idx, product, lake_uid, lake_name, lake_id,
         lat, lon, lake_area_km2, elevation_m,
         country, state, basin, found_by, is_target_lake)

# Save comprehensive summary
write.csv(all_lake_meta,
          file.path(OUTPUT_DIR, "ALL_Nechako_Lakes_Summary.csv"),
          row.names=FALSE, na = "")

# Save target lakes only
target_lakes <- all_lake_meta %>% 
  filter(is_target_lake | grepl("Nechako|Stuart|Fraser|Babine", lake_name, ignore.case=TRUE))
write.csv(target_lakes,
          file.path(OUTPUT_DIR, "TARGET_Lakes_Summary.csv"),
          row.names=FALSE, na = "")

# Console summary
cat(sprintf("\n  TOTAL UNIQUE LAKES FOUND: %d\n", nrow(all_lake_meta)))
cat(sprintf("  Target-named lakes: %d\n", sum(all_lake_meta$is_target_lake)))
if (nrow(target_lakes) > 0) {
  cat("\n  TARGET LAKES DETAILS:\n")
  print(target_lakes[, c("lake_name","lat","lon","lake_area_km2","product")], row.names=FALSE)
}


# ============================================================================
# SECTION 5C: KEY METRICS SUMMARY REPORT  [NEW - Added per user request]
# ============================================================================
cat("\n", strrep("─",60), "\n")
cat("SECTION 5C: GENERATING KEY METRICS SUMMARY\n")
cat(strrep("─",60), "\n")

# Calculate key metrics from the best result (ICESat-2)
if (!is.null(result_ice)) {
  
  # Total lakes in basin (from best product)
  total_lakes_in_basin <- result_ice$n_lakes
  
  # Calculate average volume per lake (in km³)
  # Storage is in Mm³ (million cubic metres), convert to km³ by dividing by 1e6
  # Average across all time steps and all lakes
  all_storage_values <- as.numeric(result_ice$storage_mat)
  valid_storage <- all_storage_values[!is.na(all_storage_values)]
  
  if (length(valid_storage) > 0) {
    # Average storage per lake per time step (in Mm³)
    avg_storage_per_lake_Mm3 <- mean(valid_storage) / total_lakes_in_basin
    # Convert to km³ (1 km³ = 1,000,000 Mm³)
    avg_volume_per_lake_km3 <- avg_storage_per_lake_Mm3 / 1e6
  } else {
    avg_volume_per_lake_km3 <- NA_real_
  }
  
  # Also calculate total basin average volume (all lakes combined)
  total_basin_avg_Mm3 <- mean(result_ice$basin_total$total_storage_Mm3, na.rm=TRUE)
  total_basin_avg_km3 <- total_basin_avg_Mm3 / 1e6
  
  # ── TEMPORAL SPAN CALCULATIONS ─────────────────────────────────────────────
  # Time range
  time_start <- result_ice$time_range[1]
  time_end <- result_ice$time_range[2]
  
  # Total months in the series
  total_months <- length(result_ice$dates)
  
  # Total years (including fractional)
  total_years <- as.numeric(difftime(time_end, time_start, units="days")) / 365.25
  
  # Count usable months (Q1+Q2 only)
  usable_months <- sum(result_ice$basin_total$is_absolute, na.rm=TRUE)
  
  # Count excluded months (Q3/Q4)
  excluded_months <- total_months - usable_months
  
  # Data completeness (percentage of usable months)
  data_completeness_pct <- round(100 * usable_months / total_months, 1)
  
  # Identify any gaps in the time series (months with NA storage)
  na_months <- sum(is.na(result_ice$basin_total$total_storage_Mm3))
  
  # First and last usable dates
  first_usable_date <- result_ice$basin_total$date[which(!is.na(result_ice$basin_total$total_storage_Mm3) & result_ice$basin_total$is_absolute)[1]]
  last_usable_date <- result_ice$basin_total$date[which(!is.na(result_ice$basin_total$total_storage_Mm3) & result_ice$basin_total$is_absolute)[length(which(!is.na(result_ice$basin_total$total_storage_Mm3) & result_ice$basin_total$is_absolute))]]
  
  # Create summary metrics data frame
  key_metrics <- data.frame(
    Metric = c("Analysis_Date",
               "Total_Lakes_In_Basin",
               "Total_Average_Volume_LAKES_km3",
               "Basin_Average_Volume_km3",
               "Temporal_Span_Start",
               "Temporal_Span_End",
               "Temporal_Span_Total_Months",
               "Temporal_Span_Total_Years",
               "Usable_Months_Q1_Q2",
               "Excluded_Months_Q3_Q4",
               "Data_Completeness_Pct",
               "NA_Months_In_Series",
               "First_Usable_Date",
               "Last_Usable_Date",
               "Data_Product",
               "Boundary_Method"),
    Value = c(ANALYSIS_DATE,
              as.character(total_lakes_in_basin),
              sprintf("%.6f", avg_volume_per_lake_km3),
              sprintf("%.6f", total_basin_avg_km3),
              as.character(time_start),
              as.character(time_end),
              as.character(total_months),
              sprintf("%.2f", total_years),
              as.character(usable_months),
              as.character(excluded_months),
              paste0(data_completeness_pct, "%"),
              as.character(na_months),
              as.character(first_usable_date),
              as.character(last_usable_date),
              result_ice$product,
              ifelse(!is.null(get0("BASIN_SHAPEFILE", envir=.GlobalEnv)),
                     "Shapefile (nechakoBound_dissolve.shp)",
                     "Bounding Box Only")),
    stringsAsFactors = FALSE
  )
  
  # Save to CSV
  write.csv(key_metrics,
            file.path(OUTPUT_DIR, "KEY_METRICS_Summary.csv"),
            row.names=FALSE, na = "")
  
  # Save to text file for easy reading
  sink(file.path(OUTPUT_DIR, "KEY_METRICS_Summary.txt"))
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("           NECHAKO BASIN LAKE ANALYSIS — KEY METRICS          \n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  cat(sprintf("  Analysis_Date                      : %s\n", ANALYSIS_DATE))
  cat(sprintf("  Total_Lakes_In_Basin               : %d\n", total_lakes_in_basin))
  cat(sprintf("  Total_Average_Volume_LAKES_km3     : %.6f\n", avg_volume_per_lake_km3))
  cat(sprintf("  Basin_Average_Volume_km3           : %.6f\n", total_basin_avg_km3))
  cat("\n")
  cat("  ── TEMPORAL SPAN OF AVAILABLE DATA ──────────────────────────\n")
  cat(sprintf("  Temporal_Span_Start                : %s\n", time_start))
  cat(sprintf("  Temporal_Span_End                  : %s\n", time_end))
  cat(sprintf("  Temporal_Span_Total_Months         : %d\n", total_months))
  cat(sprintf("  Temporal_Span_Total_Years          : %.2f\n", total_years))
  cat(sprintf("  Usable_Months_Q1_Q2                : %d\n", usable_months))
  cat(sprintf("  Excluded_Months_Q3_Q4              : %d\n", excluded_months))
  cat(sprintf("  Data_Completeness_Pct              : %.1f%%\n", data_completeness_pct))
  cat(sprintf("  NA_Months_In_Series                : %d\n", na_months))
  cat(sprintf("  First_Usable_Date                  : %s\n", first_usable_date))
  cat(sprintf("  Last_Usable_Date                   : %s\n", last_usable_date))
  cat("\n")
  cat("  ── DATA SOURCE ──────────────────────────────────────────────\n")
  cat(sprintf("  Data_Product                       : %s\n", result_ice$product))
  cat(sprintf("  Boundary_Method                    : %s\n",
              ifelse(!is.null(get0("BASIN_SHAPEFILE", envir=.GlobalEnv)),
                     "Shapefile (nechakoBound_dissolve.shp)",
                     "Bounding Box Only")))
  cat(sprintf("  Shapefile_Path                     : %s\n", NECHAKO_BOUNDARY))
  cat("\n═══════════════════════════════════════════════════════════════\n")
  cat("  Output Directory: %s\n", OUTPUT_DIR)
  cat("═══════════════════════════════════════════════════════════════\n")
  sink()
  
  # Display in console
  cat("\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("           NECHAKO BASIN LAKE ANALYSIS — KEY METRICS          \n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  cat(sprintf("  Analysis_Date                      : %s\n", ANALYSIS_DATE))
  cat(sprintf("  Total_Lakes_In_Basin               : %d\n", total_lakes_in_basin))
  cat(sprintf("  Total_Average_Volume_LAKES_km3     : %.6f\n", avg_volume_per_lake_km3))
  cat(sprintf("  Basin_Average_Volume_km3           : %.6f\n", total_basin_avg_km3))
  cat("\n")
  cat("  ── TEMPORAL SPAN OF AVAILABLE DATA ──────────────────────────\n")
  cat(sprintf("  Temporal_Span_Start                : %s\n", time_start))
  cat(sprintf("  Temporal_Span_End                  : %s\n", time_end))
  cat(sprintf("  Temporal_Span_Total_Months         : %d\n", total_months))
  cat(sprintf("  Temporal_Span_Total_Years          : %.2f\n", total_years))
  cat(sprintf("  Usable_Months_Q1_Q2                : %d\n", usable_months))
  cat(sprintf("  Excluded_Months_Q3_Q4              : %d\n", excluded_months))
  cat(sprintf("  Data_Completeness_Pct              : %.1f%%\n", data_completeness_pct))
  cat(sprintf("  NA_Months_In_Series                : %d\n", na_months))
  cat(sprintf("  First_Usable_Date                  : %s\n", first_usable_date))
  cat(sprintf("  Last_Usable_Date                   : %s\n", last_usable_date))
  cat("\n")
  cat("  ── DATA SOURCE ──────────────────────────────────────────────\n")
  cat(sprintf("  Data_Product                       : %s\n", result_ice$product))
  cat(sprintf("  Boundary_Method                    : %s\n",
              ifelse(!is.null(get0("BASIN_SHAPEFILE", envir=.GlobalEnv)),
                     "Shapefile (nechakoBound_dissolve.shp)",
                     "Bounding Box Only")))
  cat("\n═══════════════════════════════════════════════════════════════\n")
  cat("  Saved: KEY_METRICS_Summary.csv\n")
  cat("  Saved: KEY_METRICS_Summary.txt\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  
} else {
  cat("  WARNING: No lake data available for metrics calculation.\n")
  key_metrics <- NULL
}


# ============================================================================
# SECTION 6: USE THE BEST DATASET (ICESat‑2) FOR THE BASIN‑TOTAL SERIES
# (No merging – only the best product is retained)
# ============================================================================
cat("\n", strrep("─",60), "\n")
cat("SECTION 6: USING BEST LAKE PRODUCT (ICESat‑2) FOR BASIN TOTAL\n")
cat(strrep("─",60), "\n")

if (length(results_all) == 0) {
  cat("  No lake data found – check bounding box and file paths.\n")
  best_basin <- NULL
} else {
  # Use the ICESat‑2 product (the only one that found lakes)
  best_result <- result_ice
  best_basin <- best_result$basin_total
  
  cat(sprintf("  Using: %s\n", best_result$product))
  cat(sprintf("  %d months (%s to %s)\n",
              nrow(best_basin), min(best_basin$date), max(best_basin$date)))
  cat(sprintf("  Storage range: %.0f – %.0f Mm³  (mean %.0f Mm³)\n",
              min(best_basin$total_storage_Mm3, na.rm=TRUE),
              max(best_basin$total_storage_Mm3, na.rm=TRUE),
              mean(best_basin$total_storage_Mm3, na.rm=TRUE)))
  
  # Save a copy as the best estimate (for compatibility with SWSI script)
  write.csv(best_basin,
            file.path(OUTPUT_DIR, "BEST_ESTIMATE_basin_total_Mm3.csv"),
            row.names=FALSE, na = "")
  saveRDS(best_result, file.path(OUTPUT_DIR, "BEST_ESTIMATE.rds"))
  cat("  Saved: BEST_ESTIMATE_basin_total_Mm3.csv\n")
}


# ============================================================================
# SECTION 7: DIAGNOSTIC PLOTS
# ============================================================================
if (!is.null(best_basin) && nrow(best_basin) > 0) {
  
  # Plot 1: Basin total through time
  p1 <- ggplot(best_basin, aes(x=date, y=total_storage_Mm3)) +
    geom_line(linewidth=0.5, alpha=0.8, colour = "#1f78b4") +
    geom_point(size=0.8, alpha=0.6, colour = "#1f78b4") +
    labs(
      title    = "Nechako Basin — GloLakes Basin Storage (ICESat‑2)",
      subtitle = sprintf("Sum of %d lakes | Q1+Q2 (absolute) only | Shapefile boundary",
                         best_result$n_lakes),
      x=NULL, y="Total basin lake storage (Mm³)",
      caption="Units: Mm³ = MCM (million m³) | GloLakes v2, Hou et al. 2024"
    ) +
    theme_bw() + theme(legend.position="bottom")
  
  ggsave(file.path(OUTPUT_DIR, "plot1_basin_total_storage.png"),
         p1, width=14, height=5, dpi=150)
  
  # Plot 2: Individual lakes
  if (!is.null(result_ice)) {
    all_long_plot <- result_ice$long_df %>%
      filter(!is.na(storage_Mm3)) %>%
      mutate(lake_label = if_else(!is.na(lake_name) & nzchar(lake_name),
                                  lake_name,
                                  sprintf("(%.2fN, %.2fW)", lat, abs(lon))))
    
    if (nrow(all_long_plot) > 0) {
      p2 <- ggplot(all_long_plot,
                   aes(x=date, y=storage_Mm3, colour = as.factor(file_idx))) +
        geom_line(linewidth=0.4, show.legend = FALSE) +
        facet_wrap(~lake_label, scales="free_y") +
        labs(title = "Per‑Lake Storage Time Series — Nechako Basin (ICESat‑2)",
             subtitle = "Lakes selected using shapefile boundary",
             x = NULL, y = "Storage (Mm³)") +
        theme_bw(base_size=9) + theme(legend.position = "none")
      
      ggsave(file.path(OUTPUT_DIR, "plot2_per_lake_timeseries.png"),
             p2, width=14, height=10, dpi=150)
    }
  }
  
  # Plot 3: Seasonal cycle
  seasonal <- best_basin %>%
    group_by(month) %>%
    summarise(mean = mean(total_storage_Mm3, na.rm = TRUE),
              q10  = quantile(total_storage_Mm3, 0.10, na.rm = TRUE),
              q90  = quantile(total_storage_Mm3, 0.90, na.rm = TRUE),
              .groups = "drop")
  
  p3 <- ggplot(seasonal, aes(x=month)) +
    geom_ribbon(aes(ymin=q10, ymax=q90), fill="#1f78b4", alpha=0.2) +
    geom_line(aes(y=mean), colour="#1f78b4", linewidth=1.2) +
    geom_point(aes(y=mean), colour="#1f78b4", size=2.5) +
    scale_x_continuous(breaks=1:12, labels=month.abb) +
    labs(title = "Nechako Basin — GloLakes Seasonal Storage Cycle (ICESat‑2)",
         subtitle = "Monthly mean ± 10th/90th percentile across all years",
         x = "Month", y = "Basin storage (Mm³)") +
    theme_bw()
  
  ggsave(file.path(OUTPUT_DIR, "plot3_seasonal_cycle.png"),
         p3, width=10, height=5, dpi=150)
  
  cat("\n  Plots saved to:", OUTPUT_DIR, "\n")
}


# ============================================================================
# SECTION 8: PLUG-IN CODE FOR SWSI_Garen1993_v2.R
#
# Paste this into Section 8 of SWSI_Garen1993_v2.R, replacing the stub.
# Then set  RESERVOIR_DATA_AVAILABLE <- TRUE  at the top of that script.
# ============================================================================
cat("\n", strrep("═",65), "\n")
cat("SECTION 8: COPY THIS INTO SWSI_Garen1993_v2.R — SECTION 8\n")
cat(strrep("═",65), "\n")
cat('
# ── Replace load_reservoir_storage() stub in SWSI_Garen1993_v2.R ──────────

load_reservoir_storage <- function() {
  # GloLakes v2 — best-estimate Nechako basin lake storage (ICESat‑2)
  # Source : Hou et al. (2024) ESSD 16, 201-218, doi:10.5194/essd-16-201-2024
  # Units  : Mm3 (MCM) — already in the units used throughout this script
  # Quality: Q1+Q2 only (absolute volume); Q3/Q4 (relative) excluded
  # Boundary: Nechako basin shapefile (nechakoBound_dissolve.shp)

  f <- "D:/Nechako_Drought/Nechako/Lakes/Glolakes/nechako_extracted/BEST_ESTIMATE_basin_total_Mm3.csv"

  if (!file.exists(f))
    stop(paste("GloLakes best-estimate file not found:", f,
               "\\nRun GloLakes_Nechako_Extraction.R first."))

  raw <- read.csv(f, stringsAsFactors=FALSE)

  reservoir <- raw %>%
    mutate(
      date       = as.Date(date),
      storage    = as.numeric(total_storage_Mm3),
      year       = year(date),
      month      = month(date),
      water_year = if_else(month >= WATER_YEAR_START_MONTH, year, year - 1L)
    ) %>%
    filter(!is.na(storage), storage >= 0, n_lakes_contributing >= 1) %>%
    dplyr::select(date, storage, year, month, water_year,
                   n_lakes_contributing)

  cat(sprintf("  GloLakes storage loaded: %d months (%s to %s)\\n",
              nrow(reservoir), min(reservoir$date), max(reservoir$date)))
  cat(sprintf("  Storage range: %.0f – %.0f Mm3 | mean %.0f Mm3\\n",
              min(reservoir$storage), max(reservoir$storage),
              mean(reservoir$storage)))
  return(reservoir)
}

# Also change the flag at the top of the script:
RESERVOIR_DATA_AVAILABLE <- TRUE
')


# ============================================================================
# SECTION 9: RIVER DISCHARGE EXTRACTION (with blank instead of "NA")
# ============================================================================
cat("\n", strrep("═",65), "\n")
cat("SECTION 9: EXTRACTING RIVER DISCHARGE FOR NECHAKO\n")
cat(strrep("═",65), "\n")

extract_river_discharge <- function(nc_path, product_label) {
  cat(sprintf("\n%s\nEXTRACTING: %s\n%s\n",
              strrep("─",60), product_label, strrep("─",60)))
  if (!file.exists(nc_path)) {
    cat(sprintf("  File not found: %s\n", nc_path)); return(NULL)
  }
  
  nc <- ncdf4::nc_open(nc_path)
  on.exit(ncdf4::nc_close(nc))
  
  n_stations <- nc$dim$ID$len
  n_time     <- nc$dim$time$len
  cat(sprintf("  Stations: %d | Time steps: %d (8-day composites, 2000-present)\n",
              n_stations, n_time))
  
  lat_all <- as.numeric(ncdf4::ncvar_get(nc, "latitude"))
  lon_all <- as.numeric(ncdf4::ncvar_get(nc, "longitude"))
  
  # ── Locate Nechako stations ───────────────────────────────────────────────
  # Use shapefile boundary if available, otherwise bounding box
  if (!is.null(get0("BASIN_SHAPEFILE", envir=.GlobalEnv))) {
    basin_shp <- get("BASIN_SHAPEFILE", envir=.GlobalEnv)
    in_basin <- check_point_in_basin(lat_all, lon_all, basin_shp)
    bbox_idx <- which(in_basin)
    cat("  Using shapefile boundary for station selection\n")
  } else {
    bbox_idx <- which(
      lat_all >= BASIN_LAT_MIN & lat_all <= BASIN_LAT_MAX &
        lon_all >= BASIN_LON_MIN & lon_all <= BASIN_LON_MAX &
        !is.na(lat_all) & !is.na(lon_all)
    )
    cat("  Using bounding box for station selection\n")
  }
  
  cat(sprintf("  Stations in basin: %d\n", length(bbox_idx)))
  if (length(bbox_idx) == 0) return(NULL)
  
  # Time: 8-day composites since 2000-01-01
  time_raw <- as.numeric(ncdf4::ncvar_get(nc, "time"))
  epoch_2000 <- as.Date("2000-01-01")
  dates_8day <- epoch_2000 + as.integer(round(time_raw))
  cat(sprintf("  Date range: %s to %s\n", min(dates_8day), max(dates_8day)))
  
  # ── Read station metadata using robust string extraction ─────────────────
  river_names   <- extract_string_var(nc, "river_name", n_stations)
  station_names <- extract_string_var(nc, "station_name", n_stations)
  country_names <- extract_string_var(nc, "country_name", n_stations)
  basin_names   <- extract_string_var(nc, "basin_name", n_stations)
  
  qmean <- as.numeric(ncdf4::ncvar_get(nc, "riv_qmean"))
  qmax  <- as.numeric(ncdf4::ncvar_get(nc, "riv_qmax"))
  
  meta_riv <- data.frame(
    file_idx     = bbox_idx,
    lat          = lat_all[bbox_idx],
    lon          = lon_all[bbox_idx],
    river_name   = river_names[bbox_idx],
    station_name = station_names[bbox_idx],
    country      = country_names[bbox_idx],
    basin        = basin_names[bbox_idx],
    qmean_cms    = round(qmean[bbox_idx], 2),
    qmax_cms     = round(qmax[bbox_idx],  2),
    product      = product_label,
    stringsAsFactors = FALSE
  )
  
  cat("\n  Stations found:\n")
  # Replace NA with empty string for cleaner display
  meta_display <- meta_riv[, c("file_idx","river_name","station_name","lat","lon",
                               "qmean_cms","qmax_cms","basin")]
  meta_display[] <- lapply(meta_display, function(x) ifelse(is.na(x), "", x))
  print(meta_display, row.names=FALSE)
  
  # ── Extract 8-day discharge, aggregate to monthly mean ────────────────────
  n_match <- nrow(meta_riv)
  discharge_mat <- matrix(NA_real_, nrow=n_match, ncol=n_time)
  
  # Dimension order in river files is [time × ID]
  for (k in seq_len(n_match)) {
    idx_k <- meta_riv$file_idx[k]
    raw   <- as.numeric(ncdf4::ncvar_get(nc, "river_discharge",
                                         start=c(1, idx_k),
                                         count=c(n_time, 1)))
    raw[is.nan(raw)] <- NA_real_
    discharge_mat[k, ] <- raw
  }
  
  # Resample 8-day to monthly mean
  yr_mo         <- format(dates_8day, "%Y-%m")
  months_u      <- sort(unique(yr_mo))
  monthly_dates <- as.Date(paste0(months_u, "-01"))
  monthly_Q     <- matrix(NA_real_, nrow=n_match, ncol=length(months_u))
  
  for (j in seq_along(months_u)) {
    t_idx <- which(yr_mo == months_u[j])
    for (k in seq_len(n_match)) {
      vals <- discharge_mat[k, t_idx]
      monthly_Q[k, j] <- if (all(is.na(vals))) NA_real_ else mean(vals, na.rm=TRUE)
    }
  }
  
  # Build long-format data frame
  rows <- list()
  for (k in seq_len(n_match)) {
    rows[[k]] <- data.frame(
      date         = monthly_dates,
      year         = year(monthly_dates),
      month        = month(monthly_dates),
      water_year   = ifelse(month(monthly_dates)>=WY_START_MONTH,
                            year(monthly_dates), year(monthly_dates)-1L),
      discharge_cms = monthly_Q[k, ],
      file_idx     = meta_riv$file_idx[k],
      river_name   = meta_riv$river_name[k],
      station_name = meta_riv$station_name[k],
      lat          = meta_riv$lat[k],
      lon          = meta_riv$lon[k],
      product      = product_label,
      stringsAsFactors = FALSE
    )
  }
  long_riv <- bind_rows(rows)
  
  # Convert to monthly volume (Mm³) for comparison with WSC data
  long_riv <- long_riv %>%
    mutate(
      days_in_mo   = days_in_month(date),
      volume_Mm3   = discharge_cms * days_in_mo * 86400 / 1e6
    )
  
  # Save
  lbl <- gsub("[^A-Za-z0-9]", "_", product_label)
  write.csv(meta_riv,
            file.path(OUTPUT_DIR, paste0(lbl, "_stations.csv")),
            row.names = FALSE, na = "")
  write.csv(long_riv,
            file.path(OUTPUT_DIR, paste0(lbl, "_monthly.csv")),
            row.names = FALSE, na = "")
  saveRDS(list(meta = meta_riv, long_df = long_riv),
          file.path(OUTPUT_DIR, paste0(lbl, ".rds")))
  
  cat(sprintf("\n  Saved: %s_monthly.csv  (%d stations, %d months)\n",
              lbl, n_match, length(months_u)))
  return(list(meta = meta_riv, long_df = long_riv))
}

# Extract both river discharge products
result_msgr <- extract_river_discharge(FILE_MSGR, "MSGR_satellite")
result_gsgr <- extract_river_discharge(FILE_GSGR, "GSGR_gauge_fused")


# ============================================================================
# SECTION 10: FINAL SUMMARY
# ============================================================================
cat("\n", strrep("═",65), "\n")
cat("EXTRACTION COMPLETE\n")
cat(strrep("═",65), "\n")
cat(sprintf("Output directory: %s\n\n", OUTPUT_DIR))

cat("LAKE STORAGE OUTPUTS:\n")
# Modified pattern to include Summary files
for (f in list.files(OUTPUT_DIR, pattern="(BEST|Landsat|Summary|METRICS).*\\.csv$",
                     full.names=FALSE)) {
  sz <- file.info(file.path(OUTPUT_DIR,f))$size / 1024
  cat(sprintf("  %-55s  %.0f KB\n", f, sz))
}

cat("\nRIVER DISCHARGE OUTPUTS:\n")
for (f in list.files(OUTPUT_DIR, pattern="(MSGR|GSGR).*\\.csv$",
                     full.names=FALSE)) {
  sz <- file.info(file.path(OUTPUT_DIR,f))$size / 1024
  cat(sprintf("  %-55s  %.0f KB\n", f, sz))
}

cat("\nKEY FILES FOR ANALYSIS:\n")
cat("  KEY_METRICS_Summary.csv        ← Analysis date, lake count, avg volume, TEMPORAL SPAN\n")
cat("  KEY_METRICS_Summary.txt        ← Human-readable metrics report\n")
cat("  ALL_Nechako_Lakes_Summary.csv  ← Complete lake metadata\n")
cat("  TARGET_Lakes_Summary.csv       ← Focus lakes (Stuart, Nechako, etc.)\n")
cat("  BEST_ESTIMATE_basin_total_Mm3.csv\n")
cat("  → Plug into load_reservoir_storage() in SWSI_Garen1993_v2.R\n")
cat("  → Then set RESERVOIR_DATA_AVAILABLE <- TRUE\n")
cat("\nNEXT STEP:\n")
cat("  1. Inspect the wide CSVs in Excel to verify the lake names/locations.\n")
cat("  2. Cross-check plot1_basin_total_storage.png for data gaps.\n")
cat("  3. Check plot2_per_lake_timeseries.png for any suspiciously flat\n")
cat("     or noisy individual lake series.\n")
cat("  4. Review KEY_METRICS_Summary.txt for quick reference.\n")
cat("  5. Verify shapefile boundary was used (check found_by column in CSVs).\n")
cat("  6. Check temporal span for data coverage and completeness.\n")
# cat("  4. If Q3/Q4 periods create large gaps in the SWSI storage series,\n")
# cat("     consider WaterGAP 2.2d model output to fill those gaps.\n")