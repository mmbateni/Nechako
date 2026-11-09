# Nechako 1977 vs 1979: standalone, reproducible station/basin evidence.
# Run from Nechako project root. Requires data.table, nlme, terra.
# SOURCE DATA: Public historical BC manual snow surveys, BC manual snow survey
#   locations and ECCC AHCCD monthly station observations. External retrieval
#   is OFF by default; explicitly opt in to download missing public inputs
#   from official government services to late1970s_reference_data/. Cached files are reused;
#   previously supplied observation files are not overwritten.
#   Keep AUTO_DOWNLOAD_OBSERVATIONS <- FALSE to prohibit observation downloads.
# UPSTREAM FRAMEWORK PRODUCERS:
#   Monthly_ERALand_download.R -> monthly_data_direct/*.nc (ERA5-Land)
#   2b_PET_ERALand.R          -> PET and meteorological control NetCDFs
#   11pr_pet_trends.r          -> basin annual, monthly and native mask files
#   (11 is run in a separate R process ONLY when its outputs are absent.)
# THIS SCRIPT BUILDS ITS OWN: ERA5-AHCCD station pairs, snow-course cohort,
#   QC/complete winters, site models, 1968-1988 fixed network, and summaries.
# No preexisting *_station_pairs.csv or *_station_tests.csv is needed.
# ADDITIONAL TESTS: pooled station fixed effects and optional station slopes,
# year-level HAC and synchronized block-wild score intervals, directional
# artifact-size sensitivity, independent AHCCD observation timing, fixed
# neighbor checks, and exploratory alternative-date profiles.
# A monthly ERA5-vs-survey SWE diagnostic is NOT used for attribution because
# monthly grid means are not collocated with actual snow-course visit dates.
# Does not modify existing observations or Scripts 0-13; no PDO backcast.
# Dependency jsonlite is required ONLY when an observation must be downloaded.

suppressPackageStartupMessages({
  library(data.table)
  library(nlme)
  library(terra)
})

ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
STATIONS <- file.path(ROOT, "trend_analysis_pr_pet", "station_difference_test")
BASIN <- file.path(ROOT, "trend_analysis_pr_pet", "swe_homogeneity")
OUT <- file.path(STATIONS, "late1970s_single_script_evidence")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
BASIN_ID <- "Nechako_Basin"
SPLITS <- c(1977L, 1979L)
FIXED_YEARS <- 1968:1988
ALPHA <- 0.05
MIN_SIDE_COHORT <- 10L
AUTO_RUN_SCRIPT11_IF_MISSING <- FALSE # do not execute unreviewed upstream scripts implicitly
AUTO_DOWNLOAD_OBSERVATIONS <- FALSE # opt in explicitly to external retrieval
OBS_DIR <- file.path(ROOT, "late1970s_reference_data")
OBS_SOURCES <- c(
  allmss_archive.csv = "https://www.env.gov.bc.ca/wsd/data_searches/snow/asws/data/allmss_archive.csv",
  snow_course_locations.csv = paste0(
    "https://delivery.maps.gov.bc.ca/arcgis/rest/services/",
    "mpcm/bcgwpub/MapServer/496"),
  ahccd_monthly_bbox.csv = "https://api.weather.gc.ca/collections/ahccd-monthly/items"
)
OBS_DOWNLOADED <- character()

# Observation downloads use only public BC/ECCC services, send no station data
# or user files and do not install packages. Geographic bbox in the AHCCD query
# is derived from the project's ERA5-Land grid (terra::ext(pr_r)).
read_public_json <- function(url) {
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("Install jsonlite to retrieve missing public observations: install.packages('jsonlite')")
  tmp <- tempfile(fileext = ".json")
  on.exit(unlink(tmp), add = TRUE)
  old_timeout <- getOption("timeout")
  options(timeout = max(120L, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)
  stat <- tryCatch(suppressWarnings(utils::download.file(
    url, tmp, method = "libcurl", mode = "wb", quiet = TRUE)),
    error = function(e) stop("Public API request failed: ", conditionMessage(e), call. = FALSE))
  if (!identical(as.integer(stat), 0L) || !file.exists(tmp) || file.size(tmp) == 0)
    stop("Failed to download API JSON: ", sub("\\?.*$", "", url), call. = FALSE)
  x <- tryCatch(jsonlite::fromJSON(tmp, flatten = TRUE),
                error = function(e) stop("Invalid API JSON: ", conditionMessage(e), call. = FALSE))
  if (!is.null(x$error)) stop("Government API error: ", paste(x$error$message, collapse = " "))
  x
}

save_observation <- function(x, destination, required_cols) {
  if (!data.table::is.data.table(x)) x <- data.table::as.data.table(x)
  if (!all(required_cols %in% names(x)))
    stop("Downloaded ", basename(destination), " missing columns: ",
         paste(setdiff(required_cols, names(x)), collapse = ", "), call. = FALSE)
  if (!nrow(x)) stop("Downloaded ", basename(destination), " contains zero records")
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  temp_out <- tempfile(pattern = ".gov_snapshot_", tmpdir = dirname(destination),
                       fileext = ".csv")
  on.exit(unlink(temp_out), add = TRUE)
  data.table::fwrite(x, temp_out, na = "")
  if (file.exists(destination)) stop("Refusing to replace an existing observation: ", destination)
  if (!file.rename(temp_out, destination)) stop("Cannot move downloaded data to ", destination)
  OBS_DOWNLOADED <<- c(OBS_DOWNLOADED, basename(destination))
  message("[OBS] Saved ", nrow(x), " records: ", destination)
  destination
}

fetch_manual_snow <- function(destination) {
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = ".bc_snow_", tmpdir = dirname(destination), fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  old_timeout <- getOption("timeout")
  options(timeout = max(120L, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)
  stat <- tryCatch(suppressWarnings(utils::download.file(
    OBS_SOURCES[["allmss_archive.csv"]], tmp,
    mode = "wb", method = "libcurl", quiet = TRUE)),
    error = function(e) stop("BC snow archive download failed: ", conditionMessage(e)))
  if (!identical(as.integer(stat), 0L) || !file.exists(tmp))
    stop("Unable to download BC manual snow survey archive")
  x <- data.table::fread(tmp, showProgress = FALSE)
  save_observation(x, destination,
                   c("Number", "Date of Survey", "Water Equiv. mm", "Survey Period"))
}

fetch_snow_locations <- function(destination) {
  # Official BC Geographic Warehouse ArcGIS layer (manual surveys, active and inactive).
  # Field names checked against the published layer schema.
  endpoint <- paste0(OBS_SOURCES[["snow_course_locations.csv"]], "/query")
  pages <- list(); offset <- 0L; page_size <- 1000L
  for (i in seq_len(1000L)) {
    url <- sprintf(paste0(
      "%s?where=1%%3D1&outFields=LOCATION_ID%%2CLONGITUDE%%2CLATITUDE",
      "&returnGeometry=false&f=json&orderByFields=OBJECTID",
      "&resultOffset=%d&resultRecordCount=%d"),
      endpoint, offset, page_size)
    js <- read_public_json(url)
    f <- js$features
    if (is.null(f) || !NROW(f)) break
    d <- data.table::as.data.table(f)
    wanted <- paste0("attributes.", c("LOCATION_ID", "LONGITUDE", "LATITUDE"))
    if (!all(wanted %in% names(d)))
      stop("BC snow-location API changed its response schema")
    z <- d[, ..wanted]
    data.table::setnames(z, c("id", "lon", "lat"))
    pages[[length(pages) + 1L]] <- z
    offset <- offset + nrow(z)
    message("  BC snow locations: ", offset, " rows retrieved")
    if (nrow(z) < page_size && !isTRUE(js$exceededTransferLimit)) break
  }
  if (!length(pages)) stop("BC snow-location API returned no records")
  loc <- data.table::rbindlist(pages, fill = TRUE)
  loc[, `:=`(id = toupper(trimws(as.character(id))),
             lon = suppressWarnings(as.numeric(lon)),
             lat = suppressWarnings(as.numeric(lat)))]
  loc <- unique(loc[nzchar(id) & is.finite(lon) & is.finite(lat) &
                      lon >= -180 & lon <= 180 & lat >= -90 & lat <= 90], by = "id")
  save_observation(loc, destination, c("id", "lon", "lat"))
}

fetch_ahccd <- function(destination) {
  if (!exists("pr_r", inherits = TRUE))
    stop("ERA5 precipitation raster must be loaded before retrieving AHCCD")
  if (!terra::is.lonlat(pr_r)) stop("ERA5 raster is not in geographic coordinates")
  e <- as.vector(terra::ext(pr_r))  # xmin, xmax, ymin, ymax
  bbox <- paste(sprintf("%.5f", c(e[1L], e[3L], e[2L], e[4L])), collapse = ",")
  api <- OBS_SOURCES[["ahccd_monthly_bbox.csv"]]
  keep <- c("station_id__id_station", "date", "lat__lat", "lon__long",
            "temp_mean__temp_moyenne", "total_precip__precip_totale")
  parts <- list(); offset <- 0L; page_size <- 10000L
  for (i in seq_len(1000L)) {
    url <- sprintf("%s?f=json&bbox=%s&limit=%d&offset=%d",
                   api, bbox, page_size, offset)
    js <- read_public_json(url)
    if (is.null(js$features) || !NROW(js$features)) break
    d <- data.table::as.data.table(js$features)
    data.table::setnames(d, sub("^properties\\.", "", names(d)))
    if (!all(c("station_id__id_station", "date") %in% names(d)))
      stop("ECCC AHCCD API changed its response schema")
    for (nm in setdiff(keep, names(d))) d[, (nm) := NA_real_]
    parts[[length(parts) + 1L]] <- d[, ..keep]
    offset <- offset + nrow(d)
    message("  ECCC AHCCD: ", offset, " monthly rows retrieved")
    matched <- suppressWarnings(as.integer(js$numberMatched))
    if (length(matched) == 1L && is.finite(matched) && offset >= matched) break
    if (nrow(d) < page_size) break
  }
  if (!length(parts)) stop("No AHCCD station records overlap the ERA5 grid bbox")
  a <- data.table::rbindlist(parts, fill = TRUE)
  data.table::setnames(a, keep, c("station", "date", "lat", "lon", "tmean", "precip"))
  a[, `:=`(station = trimws(as.character(station)),
           date = as.character(date))]
  if (anyDuplicated(a, by = c("station", "date")))
    stop("AHCCD API returned duplicate station-month records; inspect pagination")
  save_observation(a, destination, c("station", "date", "lat", "lon", "tmean", "precip"))
}

# All observation inputs use the original expected file names. Cached files in
# late1970s_reference_data/, station_difference_test/, or ROOT are reused.
find_observation <- function(filename) {
  paths <- file.path(c(OBS_DIR, STATIONS, ROOT), filename)
  found <- paths[file.exists(paths)]
  if (length(found)) return(found[1L])
  if (!AUTO_DOWNLOAD_OBSERVATIONS)
    stop("Missing independent observation file: ", filename,
         "; automatic external retrieval is disabled", call. = FALSE)
  destination <- file.path(OBS_DIR, filename)
  switch(filename,
         "allmss_archive.csv" = fetch_manual_snow(destination),
         "snow_course_locations.csv" = fetch_snow_locations(destination),
         "ahccd_monthly_bbox.csv" = fetch_ahccd(destination),
         stop("Unknown observational file: ", filename))
}

# Script 11 resets its own environment; run in a SEPARATE process, never source
# it into the current analysis. Only launch when one of its outputs is missing.
require_upstream11 <- function() {
  needed <- c(
    file.path(BASIN, "swe_homogeneity_annual_series.csv"),
    file.path(ROOT, "trend_analysis_pr_pet", "regional", "regional_avg_monthly.csv"),
    file.path(ROOT, "trend_analysis_pr_pet", "native_basin_mask.rds"))
  missing <- needed[!file.exists(needed)]
  if (!length(missing)) return(invisible(FALSE))
  if (!AUTO_RUN_SCRIPT11_IF_MISSING) stop(
    "Script 11 products missing: ", paste(missing, collapse = "; "))
  producer <- file.path(ROOT, "11pr_pet_trends.r")
  if (!file.exists(producer)) stop(
    "Missing Script 11 producer: ", producer,
    "\nCopy 11pr_pet_trends.r and DROUGHT_ANALYSIS_utils.R into the project root.",
    call. = FALSE)
  # Do not invoke Monthly_ERALand_download.R: it contacts CDS and requires
  # credentials. Its required local ERA5 data must already have been acquired.
  required_nc <- file.path(ROOT, "monthly_data_direct", c(
    "total_precipitation_monthly.nc", "snow_depth_water_equivalent_monthly.nc",
    "potential_evapotranspiration_monthly.nc",
    "potential_evaporation_era5land_monthly.nc"))
  temp_candidates <- file.path(ROOT, "monthly_data_direct", c(
    "2m_temperature_monthly.nc", "t2m_monthly.nc",
    "air_temperature_2m_monthly.nc", "ERA5Land_T2m_monthly.nc"))
  absent_nc <- required_nc[!file.exists(required_nc)]
  if (!any(file.exists(temp_candidates))) absent_nc <- c(absent_nc,
                                                         paste0("one temperature file: ", paste(temp_candidates, collapse = " OR ")))
  if (length(absent_nc)) stop(
    "Cannot generate Script 11 basin products; upstream ERA5/PET files are missing:\n",
    paste("  ", absent_nc, collapse = "\n"),
    "\nRun Monthly_ERALand_download.R (ERA5 acquisition) and ",
    "2b_PET_ERALand.R (PET/control products) separately as needed.",
    call. = FALSE)
  exe <- file.path(R.home("bin"), if (.Platform$OS.type == "windows")
    "Rscript.exe" else "Rscript")
  if (!file.exists(exe)) stop("Rscript executable unavailable: ", exe)
  log_path <- file.path(OUT, "script11_auto_generation.log")
  message("Script 11 output(s) missing. Running 11pr_pet_trends.r in a separate R process. ",
          "Full original Script 11 analysis may be computationally expensive.\nLog: ", log_path)
  status <- suppressWarnings(system2(exe, shQuote(producer),
                                     stdout = log_path, stderr = log_path,
                                     wait = TRUE))
  if (!is.numeric(status) || status != 0L) stop(
    "11pr_pet_trends.r failed (exit status ", status, "). Inspect: ", log_path,
    call. = FALSE)
  if (!all(file.exists(needed))) stop(
    "Script 11 finished without all required outputs: ",
    paste(needed[!file.exists(needed)], collapse = "; "), call. = FALSE)
  invisible(TRUE)
}

read_csv <- function(path, cols, ...) {
  if (!file.exists(path)) stop("Missing input: ", path, call. = FALSE)
  x <- fread(path, ...)
  cleaned_names <- trimws(names(x))
  if (anyDuplicated(cleaned_names)) stop("Duplicate column names after trimming: ", basename(path), call. = FALSE)
  setnames(x, cleaned_names)
  missing <- setdiff(cols, names(x))
  if (length(missing)) stop("Missing fields in ", basename(path), ": ", paste(missing, collapse = ", "), call. = FALSE)
  x
}
write_table <- function(x, filename) fwrite(x, file.path(OUT, filename), na = "")

# ML AR(1) GLS; same nonmissing years in all three competing models.
# Using actual integer years preserves gaps in the AR(1) correlation distance.
fit_one <- function(x, term = NA_integer_) {
  if (!is.na(term)) x[, step := as.integer(year >= term)]
  formula <- if (is.na(term)) y ~ time else y ~ time + step
  fit <- tryCatch(
    nlme::gls(formula, data = as.data.frame(x), method = "ML",
              correlation = nlme::corAR1(form = ~ year), na.action = na.fail,
              control = nlme::glsControl(maxIter = 150L, msMaxIter = 200L)),
    error = function(e) NULL
  )
  if (is.null(fit)) return(data.table(
    converged = FALSE, aicc = NA_real_, step = NA_real_,
    ci95_lo = NA_real_, ci95_hi = NA_real_, p_nominal = NA_real_))
  k <- as.numeric(attr(logLik(fit), "df"))
  n <- nrow(x)
  aicc <- if (n > k + 1) as.numeric(AIC(fit)) + 2 * k * (k + 1) / (n - k - 1) else NA_real_
  if (is.na(term)) return(data.table(
    converged = TRUE, aicc = aicc, step = NA_real_,
    ci95_lo = NA_real_, ci95_hi = NA_real_, p_nominal = NA_real_))
  z <- summary(fit)$tTable["step", ]
  b <- unname(z["Value"])
  se <- unname(z["Std.Error"])
  crit <- qt(0.975, df = n - length(coef(fit)))
  data.table(converged = TRUE, aicc = aicc, step = b,
             ci95_lo = b - crit * se, ci95_hi = b + crit * se,
             p_nominal = unname(z["p-value"]))
}
fit_three <- function(year, values, dataset, variable, response, station, scope,
                      min_per_side = 10L) {
  x <- data.table(year = as.integer(year), y = as.numeric(values))
  x <- x[is.finite(year) & is.finite(y)]
  setorder(x, year)
  if (anyDuplicated(x$year)) stop("Duplicate years: ", dataset, "/", variable, "/", station)
  if (nrow(x) < 2L * min_per_side + 2L ||
      any(vapply(SPLITS, function(s) min(sum(x$year < s), sum(x$year >= s)) < min_per_side,
                 logical(1L)))) return(data.table())
  x[, time := year - min(year)]
  specs <- list(trend_only = NA_integer_, trend_step_1977 = 1977L,
                trend_step_1979 = 1979L)
  z <- rbindlist(lapply(names(specs), function(nm) {
    r <- fit_one(copy(x), specs[[nm]])
    r[, `:=`(model = nm, n = nrow(x),
             dataset = dataset, variable = variable, response = response,
             station = station, scope = scope)]
    r
  }), use.names = TRUE)
  setcolorder(z, c("dataset", "variable", "response", "station", "scope",
                   "model", "n", "converged", "aicc", "step", "ci95_lo",
                   "ci95_hi", "p_nominal"))
  z
}
model_comparison <- function(z) {
  keys <- c("dataset", "variable", "response", "station", "scope")
  if (!nrow(z)) return(data.table())
  w <- dcast(z, dataset + variable + response + station + scope ~ model,
             value.var = "aicc")
  needed <- c("trend_only", "trend_step_1977", "trend_step_1979")
  for (nm in setdiff(needed, names(w))) w[, (nm) := NA_real_]
  w[, `:=`(aicc_1979_minus_1977 = trend_step_1979 - trend_step_1977,
           aicc_1977_minus_trend = trend_step_1977 - trend_only,
           aicc_1979_minus_trend = trend_step_1979 - trend_only)]
  w[]
}

# 1. UPSTREAM BASIN PRODUCTS AND NATIVE GRID.
# The local Script 11 outputs can be regenerated without any dependency on the
# prior difference-test scripts. These files are the authoritative basin metrics.
require_upstream11()
monthly_dir <- file.path(ROOT, "monthly_data_direct")
pr_path <- file.path(monthly_dir, "total_precipitation_monthly.nc")
t_path <- file.path(monthly_dir, "2m_temperature_monthly.nc")
if (!file.exists(t_path)) {
  candidates <- file.path(monthly_dir, c("t2m_monthly.nc",
                                         "air_temperature_2m_monthly.nc", "ERA5Land_T2m_monthly.nc"))
  candidate <- candidates[file.exists(candidates)]
  if (!length(candidate)) stop("No ERA5-Land 2 m temperature NetCDF found")
  t_path <- candidate[1L]
}
swe_path <- file.path(monthly_dir, "snow_depth_water_equivalent_monthly.nc")
if (!all(file.exists(c(pr_path, t_path, swe_path)))) stop(
  "Missing monthly ERA5-Land precipitation, T2m or SWE NetCDF. ",
  "Run the existing monthly producer first.")
pr_r <- terra::rast(pr_path)
t_r <- terra::rast(t_path)
swe_r <- terra::rast(swe_path)
if (!terra::is.lonlat(pr_r)) stop("ERA5 precipitation grid must use geographic coordinates")
if (!terra::compareGeom(pr_r[[1]], t_r[[1]], swe_r[[1]], stopOnError = FALSE))
  stop("ERA5 precipitation, temperature, and SWE grids have different geometry")
if (length(unique(c(terra::nlyr(pr_r), terra::nlyr(t_r), terra::nlyr(swe_r)))) != 1L)
  stop("Monthly ERA5 variables have different numbers of layers")

monthly_dates <- function(r) {
  d <- tryCatch(as.Date(terra::time(r)), error = function(e) NULL)
  if (is.null(d) || length(d) != terra::nlyr(r) || anyNA(d)) {
    if (terra::nlyr(r) %% 12L != 0L)
      stop("Undated monthly NetCDF with non-annual layer count")
    d <- seq(as.Date("1950-01-01"), by = "month", length.out = terra::nlyr(r))
    warning("ERA5 NetCDF has no usable time axis; applying Script 11's 1950-01 calendar convention")
  }
  as.Date(format(d, "%Y-%m-01"))
}
dates <- monthly_dates(pr_r)
if (!identical(monthly_dates(t_r), dates) || !identical(monthly_dates(swe_r), dates))
  stop("ERA5-Land monthly input calendars do not agree")
if (anyDuplicated(dates) || any(diff(as.integer(format(dates, "%Y")) * 12L +
                                     as.integer(format(dates, "%m"))) != 1L))
  stop("ERA5-Land monthly calendar is not strictly consecutive")

# Use exactly the native basin membership written by Script 11, not a guessed
# bounding box or a coarse polygon substitute.
nmask <- readRDS(file.path(ROOT, "trend_analysis_pr_pet", "native_basin_mask.rds"))
if (!is.list(nmask) || is.null(nmask$template) || is.null(nmask$basin_mask))
  stop("Script 11 native_basin_mask.rds has an unknown structure")
mask_r <- if (inherits(nmask$template, "PackedSpatRaster"))
  terra::unwrap(nmask$template) else nmask$template
if (!inherits(mask_r, "SpatRaster")) stop("Cannot restore Script 11 basin raster")
if (length(nmask$basin_mask) != terra::ncell(mask_r))
  stop("Script 11 native basin mask has wrong cell count")
terra::values(mask_r) <- as.numeric(nmask$basin_mask)
if (!terra::same.crs(mask_r, pr_r)) stop("Script 11 native mask and ERA5 grid have different CRS")

site_in_basin <- function(lon, lat) {
  z <- terra::extract(mask_r, cbind(as.numeric(lon), as.numeric(lat)))
  v <- as.numeric(z[[ncol(z)]])
  is.finite(v) & v == 1
}
site_cells <- function(lon, lat) {
  terra::cellFromXY(pr_r[[1]], cbind(as.numeric(lon), as.numeric(lat)))
}
extract_cells <- function(r, cells) {
  z <- as.matrix(terra::extract(r, as.integer(cells)))
  if (ncol(z) == terra::nlyr(r) + 1L) z <- z[, -1L, drop = FALSE]
  if (ncol(z) != terra::nlyr(r) || nrow(z) != length(cells))
    stop("Unexpected terra::extract dimensions")
  storage.mode(z) <- "double"
  z
}
winter_dt <- function(vec, fun) {
  mon <- as.integer(format(dates, "%m"))
  yr <- as.integer(format(dates, "%Y"))
  choose <- mon %in% c(11L, 12L, 1L, 2L, 3L)
  z <- data.table(year = yr[choose] + as.integer(mon[choose] >= 11L),
                  mon = mon[choose], v = as.numeric(vec[choose]))
  z[, .(value = if (.N == 5L && uniqueN(mon) == 5L && all(is.finite(v)))
    fun(v) else NA_real_), by = year]
}

# Read independent AHCCD observations from supplied snapshot, validate and
# reconstruct complete five-month station winters before computing differences.
raw <- read_csv(find_observation("ahccd_monthly_bbox.csv"),
                c("station", "date", "lat", "lon", "tmean", "precip"),
                colClasses = c(station = "character"))
raw[, `:=`(station = trimws(station),
           date = as.Date(paste0(substr(as.character(date), 1L, 7L), "-01")),
           lat = as.numeric(lat), lon = as.numeric(lon),
           tmean = suppressWarnings(as.numeric(tmean)),
           precip = suppressWarnings(as.numeric(precip)))]
if (anyNA(raw$date) || anyNA(raw$station) || any(!nzchar(raw$station)) ||
    anyDuplicated(raw, by = c("station", "date")))
  stop("AHCCD archive contains invalid IDs/dates or repeated station-months")
raw[, `:=`(t_missing = is.finite(tmean) & tmean <= -999,
           p_missing = is.finite(precip) & precip < 0)]
write_table(data.table(variable = c("temperature", "precipitation"),
                       n_sentinel_months = c(sum(raw$t_missing), sum(raw$p_missing))),
            "AHCCD_missing_value_audit.csv")
raw[t_missing == TRUE, tmean := NA_real_]
raw[p_missing == TRUE, precip := NA_real_]
if (raw[is.finite(tmean), any(tmean < -65 | tmean > 55)] ||
    raw[is.finite(precip), any(precip > 2000)])
  stop("Other implausible AHCCD values found; inspect input")
raw[, `:=`(cyear = as.integer(format(date, "%Y")),
           mon = as.integer(format(date, "%m")))]

# Coordinates are taken from the cached station record (median of valid sites).
# AHCCD monthly records generally repeat station metadata: coordinate stability
# cannot establish that a station was never relocated. This flags only observable
# variation in coordinates within this snapshot, not relocation history.
ahccd_locs <- raw[is.finite(lon) & is.finite(lat),
                  .(lon = median(lon), lat = median(lat),
                    lon_range = diff(range(lon)), lat_range = diff(range(lat))), by = station]
ahccd_locs[, `:=`(cell = site_cells(lon, lat),
                  in_basin = site_in_basin(lon, lat),
                  coordinate_variation_flag = lon_range > 0.02 | lat_range > 0.02)]
write_table(ahccd_locs, "AHCCD_station_grid_and_coordinate_variation_QA.csv")
ahccd_locs <- ahccd_locs[!is.na(cell)]
if (!nrow(ahccd_locs)) stop("No AHCCD station lies inside the native ERA5 domain")
# Unlike snow courses, AHCCD sites are selected within the ERA5 domain, not
# exclusively inside the native Nechako basin mask; retain this distinction.

# Convert ERA5 as in Script 11 / original station comparison:
# precipitation mean daily rate [mm/day] -> monthly total [mm];
# temperature [K] -> [degC]. Avoid unit guessing from station differences.
pr_sample <- as.numeric(terra::values(pr_r[[which(as.integer(format(dates, "%m")) == 7L)[1L]]], mat = FALSE))
pr_sample <- pr_sample[is.finite(pr_sample)]
# Both ends matter: an input in m/day can otherwise pass a maximum-only check.
# The July median threshold is a heuristic guard, not a metadata-based unit test.
if (!length(pr_sample) || stats::median(pr_sample) > 50 ||
    stats::median(pr_sample) <= 0.1) stop("ERA5 July precipitation not plausibly in mm/day")
# Month-length calculation uses date differences (including leap years).
next_month <- seq(as.Date(dates[1L]), by = "month", length.out = length(dates) + 1L)[-1L]
days_in_month <- as.numeric(next_month - dates)
if (any(!is.finite(days_in_month) | days_in_month < 28 | days_in_month > 31))
  stop("Failed to infer calendar month lengths")
make_era_winter <- function(locs) {
  cells <- locs$cell
  p <- sweep(extract_cells(pr_r, cells), 2L, days_in_month, `*`)
  t <- extract_cells(t_r, cells)
  # Explicit physical units: the project ERA5-Land 2m_temperature is Kelvin.
  t <- t - 273.15
  if (any(t[is.finite(t)] < -100 | t[is.finite(t)] > 70))
    stop("ERA5 temperature conversion outside plausible degC range")
  list(
    temperature = rbindlist(lapply(seq_len(nrow(locs)), function(i) {
      d <- winter_dt(t[i, ], mean); d[, station := locs$station[i]]; d
    })),
    precipitation = rbindlist(lapply(seq_len(nrow(locs)), function(i) {
      d <- winter_dt(p[i, ], sum); d[, station := locs$station[i]]; d
    })))
}
era_winter <- make_era_winter(ahccd_locs)

make_station_pairs <- function(variable, obs_col, fun) {
  old <- copy(era_winter[[variable]])
  setnames(old, "value", "value_era")
  m <- raw[mon %in% c(11L, 12L, 1L, 2L, 3L),
           .(station, mon, year = cyear + as.integer(mon >= 11L),
             v = get(obs_col))]
  obs <- m[, .(n_valid_months = sum(is.finite(v)),
               observed = if (.N == 5L && uniqueN(mon) == 5L &&
                              all(is.finite(v))) fun(v) else NA_real_),
           by = .(station, year)]
  z <- merge(old, obs, by = c("station", "year"), all.x = TRUE)
  z[, `:=`(difference = NA_real_,
           eligible = is.finite(value_era) & is.finite(observed),
           variable = variable)]
  if (variable == "precipitation") {
    z[eligible == TRUE & (value_era <= 0 | observed <= 0), eligible := FALSE]
    z[eligible == TRUE, difference := log(value_era / observed)]
  } else {
    z[eligible == TRUE, difference := value_era - observed]
  }
  if (anyDuplicated(z, by = c("station", "year")))
    stop("Derived duplicate station-year pairs: ", variable)
  write_table(z, paste0(variable, "_generated_NovMar_station_pairs_QC.csv"))
  message(variable, ": generated ", sum(z$eligible), " complete station-winter pairs")
  z[eligible == TRUE]
}
tpair <- make_station_pairs("temperature", "tmean", mean)
ppair <- make_station_pairs("precipitation", "precip", sum)
pairs <- rbindlist(list(tpair, ppair), use.names = TRUE, fill = TRUE)
if (!nrow(pairs)) stop("No complete AHCCD-ERA5 station-winter pairs")

# 2. BC MANUAL SNOW COURSES: derive cohort from the native basin mask and
# sufficient independently observed pre/post coverage, not old test outputs.
sc <- read_csv(find_observation("allmss_archive.csv"),
               c("Number", "Date of Survey", "Water Equiv. mm", "Survey Period"))
sc[, `:=`(station = toupper(gsub("\\s+", "", as.character(Number))),
          date = as.Date(`Date of Survey`, format = "%Y/%m/%d"),
          swe_mm = suppressWarnings(as.numeric(`Water Equiv. mm`)),
          period = trimws(as.character(`Survey Period`)))]
sc[, `:=`(year = as.integer(format(date, "%Y")),
          target_month = as.integer(ifelse(period == "01-Mar", 3L,
                                           ifelse(period == "01-Apr", 4L, NA_integer_))))]
sc <- sc[year %in% 1950:2025 & !is.na(target_month) &
           is.finite(swe_mm) & swe_mm >= 0 & nzchar(station)]
sc[, target_date := as.Date(sprintf("%04d-%02d-01", year, target_month))]
sc[, offset := abs(as.integer(date - target_date))]
sc <- sc[offset <= 10L]
setorder(sc, station, year, target_month, offset, date)
sc <- sc[, .SD[1L], by = .(station, year, target_month)]
snow <- sc[, .(n_periods = uniqueN(target_month),
               swe = if (.N == 2L && uniqueN(target_month) == 2L)
                 mean(swe_mm) else NA_real_), by = .(station, year)]
snow <- snow[n_periods == 2L & is.finite(swe)]

loc <- read_csv(find_observation("snow_course_locations.csv"),
                c("id", "lon", "lat"))
loc[, `:=`(station = toupper(gsub("\\s+", "", as.character(id))),
           lon = as.numeric(lon), lat = as.numeric(lat))]
loc <- unique(loc[is.finite(lon) & is.finite(lat),
                  .(station, lon, lat)], by = "station")
loc[, `:=`(cell = site_cells(lon, lat), in_basin = site_in_basin(lon, lat))]
write_table(loc, "snow_course_spatial_membership_QA.csv")
loc <- loc[!is.na(cell) & in_basin == TRUE]
snow <- snow[station %in% loc$station]
cohort <- snow[, .(before = sum(year < 1979L),
                   after = sum(year >= 1979L),
                   first = min(year), last = max(year)), by = station]
cohort[, eligible := before >= MIN_SIDE_COHORT & after >= MIN_SIDE_COHORT]
write_table(cohort, "snowcourse_generated_cohort_coverage.csv")
snow <- snow[station %in% cohort[eligible == TRUE, station]]
if (!nrow(snow)) stop("No in-basin snow courses have complete Mar/Apr visits on both sides of 1979")
write_table(snow, "snowcourse_complete_marapr.csv")
snow_fits <- rbindlist(lapply(split(snow, by = "station", keep.by = TRUE), function(x)
  fit_three(x$year, x$swe, "Snow_course", "SWE_MarApr", "observed",
            x$station[1L], "site", 10L)), fill = TRUE)
if (!nrow(snow_fits)) stop("Snow-course station coverage too short for pre/post tests")
snow_fits[, p_bh := p.adjust(p_nominal, method = "BH"), by = model]
write_table(snow_fits, "snowcourse_AR1_models.csv")
snow_cmp <- model_comparison(snow_fits)
write_table(snow_cmp, "snowcourse_AICc_comparison.csv")

# Optional snow-course fixed network: full 1968-1988 and pre-1977 standardization.
sc_fixed <- snow[year %in% FIXED_YEARS]
sc_keep <- sc_fixed[, .(n_years = uniqueN(year)), by = station][n_years == length(FIXED_YEARS), station]
snow_fixed_fits <- data.table()
if (length(sc_keep) >= 2L) {
  z <- copy(sc_fixed[station %in% sc_keep])
  z[, `:=`(baseline_mean = mean(swe[year <= 1976]),
           baseline_sd = sd(swe[year <= 1976])), by = station]
  z <- z[is.finite(baseline_sd) & baseline_sd > 0]
  if (uniqueN(z$station) >= 2L) {
    z[, standardized := (swe - baseline_mean) / baseline_sd]
    C <- z[, .(value = mean(standardized), n_sites = .N), by = year]
    C <- C[n_sites == uniqueN(z$station)]
    snow_fixed_fits <- fit_three(C$year, C$value, "Snow_course", "SWE_MarApr",
                                 "observed_z", "FIXED_1968_1988", "fixed", 7L)
  }
}
write_table(snow_fixed_fits, "snowcourse_fixed_network_models.csv")
write_table(data.table(n_complete_stations = length(sc_keep),
                       estimable = nrow(snow_fixed_fits) == 3L),
            "snowcourse_fixed_network_status.csv")

# 3. ERA5-AHCCD DIFFERENCE TESTS (freshly generated station pairs above).
# Site-level differences retain the heterogeneous-site diagnostic.
site_diff <- rbindlist(lapply(split(pairs, by = c("variable", "station"), keep.by = TRUE), function(x)
  fit_three(x$year, x$difference, "ERA5_minus_AHCCD", x$variable[1L],
            "difference", x$station[1L], "site", 10L)), fill = TRUE)
if (nrow(site_diff)) {
  site_diff[, p_bh := p.adjust(p_nominal, method = "BH"), by = .(variable, model)]
  write_table(site_diff, "AHCCD_site_difference_AR1_models.csv")
}

# Fixed 21-winter network excludes membership changes, but has limited power.
fixed_models <- list()
fixed_member <- list()
for (v in c("temperature", "precipitation")) {
  d <- pairs[variable == v & year %in% FIXED_YEARS]
  coverage <- d[, .(n_years = uniqueN(year)), by = station]
  keep <- coverage[n_years == length(FIXED_YEARS), station]
  fixed_member[[v]] <- data.table(variable = v, station = coverage$station,
                                  n_years = coverage$n_years,
                                  included = coverage$station %in% keep)
  if (length(keep) < 2L) {
    warning(v, ": fewer than 2 stations complete in 1968-1988")
    next
  }
  d <- copy(d[station %in% keep])
  for (col in c("value_era", "observed", "difference")) {
    nc <- paste0(col, "_centred")
    d[, (nc) := get(col) - mean(get(col)[year <= 1976L]), by = station]
  }
  C <- d[, .(ERA5 = mean(value_era_centred),
             AHCCD = mean(observed_centred),
             difference = mean(difference_centred),
             n_sites = .N), by = year]
  if (nrow(C) != length(FIXED_YEARS) || any(C$n_sites != length(keep)))
    stop("Fixed-network coverage assertion failed: ", v)
  for (name in c("ERA5", "AHCCD", "difference")) {
    fixed_models[[paste(v, name)]] <- fit_three(
      C$year, C[[name]], "ERA5_AHCCD", v, name,
      "FIXED_1968_1988", "fixed", 7L)
  }
}
write_table(rbindlist(fixed_member, fill = TRUE), "AHCCD_fixed_network_membership.csv")
fixed_fits <- rbindlist(fixed_models, fill = TRUE)
write_table(fixed_fits, "AHCCD_fixed_network_AR1_models.csv")
fixed_cmp <- model_comparison(fixed_fits)
write_table(fixed_cmp, "AHCCD_fixed_network_AICc_comparison.csv")

# 4. ERA5-LAND BASIN annual SWE and related series (internal evidence only).
# PET_ERA5_annual is excluded pending confirmation of its physical variable
# provenance; no conclusion is drawn about that series here.
annual <- read_csv(file.path(BASIN, "swe_homogeneity_annual_series.csv"),
                   c("region_id", "series", "year", "value"))
annual <- annual[region_id == BASIN_ID & series %in%
                   c("SWE_MarApr", "SWE_Peak", "Precip_NovMar", "Temp_NovMar") &
                   year %in% 1950:2025, .(series, year, value)]
extra_path <- file.path(ROOT, "trend_analysis_pr_pet", "homogeneity",
                        "homogeneity_annual_series_all.csv")
if (file.exists(extra_path)) {
  extra <- read_csv(extra_path, c("region_id", "series", "year", "value"))
  extra <- extra[region_id == BASIN_ID & series %in%
                   c("AET_annual", "AET_resid_P_PET",
                     "Snowfall_annual", "SnowDepth_MarApr") &
                   year %in% 1950:2025, .(series, year, value)]
  annual <- rbindlist(list(annual, extra), use.names = TRUE)
}
if (!nrow(annual)) stop("No matching basin annual records")
annual_fits <- rbindlist(lapply(split(annual, by = "series", keep.by = TRUE), function(x)
  fit_three(x$year, x$value, "ERA5_Land", x$series[1L], "basin",
            BASIN_ID, "basin", 12L)), fill = TRUE)
write_table(annual_fits, "ERA5_basin_annual_AR1_models.csv")
annual_cmp <- model_comparison(annual_fits)
write_table(annual_cmp, "ERA5_basin_annual_AICc_comparison.csv")

# 5. Calendar-month fingerprint, avoiding annual/melt-season conflation.
# Monthly SWE in mostly snow-free months is a presence/absence proxy and
# Gaussian GLS inference is exploratory. FDR is within variable and split.
monthly_path <- file.path(ROOT, "trend_analysis_pr_pet", "regional", "regional_avg_monthly.csv")
monthly_fits <- data.table()
monthly_QA <- data.table()
if (file.exists(monthly_path)) {
  mm <- read_csv(monthly_path,
                 c("region_id", "date", "precip_mm_month", "tair_degC_month", "swe_mm_month"))
  mm <- mm[region_id == BASIN_ID]
  mm[, date := as.Date(date)]
  if (anyNA(mm$date) || anyDuplicated(mm, by = "date"))
    stop("Regional monthly dates invalid or duplicated")
  mm[, `:=`(year = as.integer(format(date, "%Y")),
            month = as.integer(format(date, "%m")))]
  mm <- mm[year %in% 1950:2025]
  month_cols <- c("swe_mm_month", "precip_mm_month", "tair_degC_month")
  L <- melt(mm[, c("year", "month", month_cols), with = FALSE],
            id.vars = c("year", "month"), measure.vars = month_cols,
            variable.name = "variable", value.name = "value")
  L[, value := as.numeric(value)]
  monthly_QA <- L[variable == "swe_mm_month" & is.finite(value),
                  .(n_years = .N, frac_exact_zero = mean(value == 0),
                    frac_near_zero_0p1mm = mean(value <= 0.1),
                    median_swe_mm = median(value)), by = month]
  write_table(monthly_QA, "SWE_monthly_near_zero_QA.csv")
  monthly_fits <- rbindlist(lapply(split(L, by = c("variable", "month"), keep.by = TRUE), function(x)
    fit_three(x$year, x$value,
              "ERA5_Land", x$variable[1L], "monthly",
              paste0("month_", sprintf("%02d", x$month[1L])), "calendar_month", 12L)),
    fill = TRUE)
  if (nrow(monthly_fits)) {
    monthly_fits[, p_bh := p.adjust(p_nominal, method = "BH"),
                 by = .(variable, model)]
    write_table(monthly_fits, "ERA5_calendar_month_AR1_models.csv")
  }
} else {
  warning("Monthly fingerprint skipped: ", monthly_path, " not found")
}

# 6. POOLED STATION EVIDENCE WITH YEAR-DEPENDENT UNCERTAINTY.
# A shared regional signal makes station-years non-independent. A station-FE
# panel is fitted using the whole eligible record, but uncertainty comes from
# year-aggregated score vectors with a Bartlett HAC (Newey-West) correction.
# All sites in a given year contribute to the SAME score vector. This handles
# arbitrary station entry/exit without treating a station-year as a new year.
# A synchronized 5-year block wild-score bootstrap is a sensitivity check;
# it keeps the original missingness pattern and multiplies all site scores in
# each calendar-year block by one common Rademacher random variable.
# Neither adjustment remedies differential site-specific nonlinear trends.

PANEL_MIN_EACH_SIDE <- 5L
PANEL_MIN_TOTAL <- 16L
PANEL_HAC_LAG <- 5L
PANEL_BLOCK_LENGTH <- 5L
PANEL_BOOT_REPS <- 999L
PANEL_PLACEBO_YEARS <- 1960:1995
PANEL_SEED <- 20261105L
set.seed(PANEL_SEED)

panel_cohort <- function(d, min_side = PANEL_MIN_EACH_SIDE) {
  if (!all(c("station", "year", "y") %in% names(d)))
    stop("Pooled input must have station, year, y")
  x <- copy(d[, .(station = as.character(station),
                  year = as.integer(year), y = as.numeric(y))])
  x <- x[!is.na(station) & nzchar(station) & is.finite(year) &
           is.finite(y) & year >= 1950L & year <= 2025L]
  if (anyDuplicated(x, by = c("station", "year")))
    stop("Duplicated station-year in pooled analysis")
  cov <- x[, .(n = .N, first_year = min(year), last_year = max(year),
               pre77 = sum(year < 1977L), post77 = sum(year >= 1977L),
               pre79 = sum(year < 1979L), post79 = sum(year >= 1979L)),
           by = station]
  keep <- cov[n >= PANEL_MIN_TOTAL & pre77 >= min_side & post77 >= min_side &
                pre79 >= min_side & post79 >= min_side, station]
  setorder(x, year, station)
  x <- x[station %in% keep]
  x[, station_f := factor(station)]
  x[, time := year - 1950L]
  attr(x, "coverage") <- cov[, eligible := station %in% keep]
  x
}

# This uses a constant model matrix across candidate step years, with actual
# years rather than observation-row indices as the time axis.
panel_fit <- function(x, split = NA_integer_, site_trends = FALSE,
                      hac_lag = PANEL_HAC_LAG, boot_reps = 0L) {
  empty <- data.table(n_station = uniqueN(x$station), n_years = uniqueN(x$year),
                      n_station_years = nrow(x), converged = FALSE,
                      step = NA_real_, se_year_hac = NA_real_,
                      ci95_lo = NA_real_, ci95_hi = NA_real_,
                      p_year_hac = NA_real_, ci90_lo = NA_real_, ci90_hi = NA_real_,
                      year_hac_lag = hac_lag, iid_aicc_descriptive = NA_real_,
                      ci95_lo_block_wild = NA_real_, ci95_hi_block_wild = NA_real_,
                      ci90_lo_block_wild = NA_real_, ci90_hi_block_wild = NA_real_,
                      n_blocks = NA_integer_, df_years = NA_integer_)
  if (nrow(x) < 30L || uniqueN(x$station) < 2L || uniqueN(x$year) < 25L)
    return(empty)
  d <- copy(x)
  if (!is.na(split)) {
    if (min(sum(d$year < split), sum(d$year >= split)) < 8L ||
        min(uniqueN(d$year[d$year < split]),
            uniqueN(d$year[d$year >= split])) < 8L) return(empty)
    d[, step_indicator := as.integer(year >= split)]
  }
  formula <- if (site_trends) {
    if (is.na(split)) y ~ 0 + station_f + station_f:time else
      y ~ 0 + station_f + station_f:time + step_indicator
  } else {
    if (is.na(split)) y ~ 0 + station_f + time else
      y ~ 0 + station_f + time + step_indicator
  }
  mod <- tryCatch(lm(formula, data = d, na.action = na.fail),
                  error = function(e) NULL)
  if (is.null(mod)) return(empty)
  X <- model.matrix(mod)
  if (qr(X)$rank != ncol(X)) return(empty)
  k <- ncol(X) + 1L # count residual variance for Gaussian AICc
  n <- nrow(X)
  empty[, iid_aicc_descriptive := if (n > k + 1L)
    AIC(mod) + 2 * k * (k + 1) / (n - k - 1) else NA_real_]
  empty[, converged := TRUE]
  if (is.na(split)) return(empty)
  j <- match("step_indicator", colnames(X))
  if (is.na(j) || !is.finite(coef(mod)[j])) return(empty)
  
  # Aggregate all station residual score contributions into a *single* year.
  # Calendar gaps are represented by zero scores; the lag is in years.
  Urow <- sweep(X, 1L, as.numeric(residuals(mod)), `*`)
  U_by_year <- rowsum(Urow, group = d$year, reorder = TRUE)
  yrs <- seq.int(min(d$year), max(d$year))
  U <- matrix(0, nrow = length(yrs), ncol = ncol(X),
              dimnames = list(as.character(yrs), colnames(X)))
  U[rownames(U_by_year), ] <- U_by_year
  Tobs <- uniqueN(d$year)
  TT <- nrow(U)
  lag <- min(as.integer(hac_lag), TT - 1L)
  S <- crossprod(U)
  if (lag > 0L) for (h in seq_len(lag)) {
    A <- crossprod(U[(h + 1L):TT, , drop = FALSE],
                   U[seq_len(TT - h), , drop = FALSE])
    S <- S + (1 - h / (lag + 1L)) * (A + t(A))
  }
  bread <- tryCatch(chol2inv(chol(crossprod(X))), error = function(e) NULL)
  if (is.null(bread)) return(empty)
  # Finite-year correction; the time unit for degrees of freedom is years.
  vc <- bread %*% S %*% bread * (Tobs / (Tobs - 1L)) *
    ((n - 1L) / (n - ncol(X)))
  b <- unname(coef(mod)[j])
  se <- sqrt(max(vc[j, j], 0))
  if (!is.finite(se) || se <= 0) return(empty)
  df <- Tobs - 1L
  crit95 <- qt(0.975, df)
  crit90 <- qt(0.95, df)
  empty[, `:=`(step = b, se_year_hac = se,
               ci95_lo = b - crit95 * se, ci95_hi = b + crit95 * se,
               ci90_lo = b - crit90 * se, ci90_hi = b + crit90 * se,
               p_year_hac = 2 * pt(-abs(b / se), df = df),
               df_years = df)]
  
  if (boot_reps > 0L) {
    # Year-synchronous block wild-score bootstrap. Resamples regression score
    # shocks, not rows or residuals from mismatched site/year observations.
    # The observation design and station coverage remain unchanged.
    beta_draws <- vapply(seq_len(boot_reps), function(bi) {
      shift <- sample.int(PANEL_BLOCK_LENGTH, 1L) - 1L
      groups <- as.integer((seq_len(TT) + shift - 1L) %/% PANEL_BLOCK_LENGTH)
      G <- rowsum(U, group = groups, reorder = TRUE)
      mult <- sample(c(-1, 1), nrow(G), replace = TRUE)
      b + as.numeric((bread %*% colSums(sweep(G, 1L, mult, `*`)))[j])
    }, numeric(1L))
    qs <- quantile(beta_draws, probs = c(.025, .05, .95, .975),
                   names = FALSE, na.rm = TRUE)
    empty[, `:=`(ci95_lo_block_wild = qs[1L], ci90_lo_block_wild = qs[2L],
                 ci90_hi_block_wild = qs[3L], ci95_hi_block_wild = qs[4L],
                 n_blocks = ceiling(TT / PANEL_BLOCK_LENGTH))]
  }
  empty
}

# Directional null tests exclude one *specified* hypothesized artifact size.
# This is NOT a two-sided equivalence claim. A separate symmetric TOST is
# calculated only as a sensitivity test for a user-defined +/-margin.
margin_test <- function(row, margin) {
  if (nrow(row) != 1L || !is.finite(margin) || margin == 0 ||
      !is.finite(row$step) || !is.finite(row$se_year_hac) ||
      row$se_year_hac <= 0 || !is.finite(row$df_years)) {
    return(data.table(artifact_margin = margin,
                      p_exclude_directional_hac = NA_real_,
                      exclude_directional_95_hac = NA,
                      exclude_directional_95_block_wild = NA,
                      p_tost_symmetric_hac = NA_real_,
                      equivalent_within_pm_margin_hac = NA))
  }
  b <- row$step; se <- row$se_year_hac; df <- row$df_years
  p_dir <- if (margin < 0)
    pt((b - margin) / se, df = df, lower.tail = FALSE) else
      pt((b - margin) / se, df = df, lower.tail = TRUE)
  b_dir <- if (margin < 0)
    is.finite(row$ci90_lo_block_wild) && row$ci90_lo_block_wild > margin else
      is.finite(row$ci90_hi_block_wild) && row$ci90_hi_block_wild < margin
  a <- abs(margin)
  p_tost <- max(pt((b + a) / se, df = df, lower.tail = FALSE),
                pt((b - a) / se, df = df, lower.tail = TRUE))
  data.table(artifact_margin = margin,
             p_exclude_directional_hac = p_dir,
             exclude_directional_95_hac = is.finite(p_dir) && p_dir < ALPHA,
             exclude_directional_95_block_wild = b_dir,
             p_tost_symmetric_hac = p_tost,
             equivalent_within_pm_margin_hac = is.finite(p_tost) && p_tost < ALPHA)
}

# Observation and difference series are kept in their native units.
# Pool physical snow-course SWE in mm (site-specific levels absorbed by FE);
# also report site-standardized SWE as a secondary sensitivity analysis.
snow_z <- copy(snow)
snow_z[, `:=`(pre_mean = mean(swe[year < 1977L]),
              pre_sd = sd(swe[year < 1977L])), by = station]
snow_z <- snow_z[is.finite(pre_sd) & pre_sd > 0]
snow_z[, y := (swe - pre_mean) / pre_sd]

pooled_sets <- list(
  observed_snowcourse_mm = snow[, .(station, year, y = swe)],
  observed_snowcourse_pre1977_sd = snow_z[, .(station, year, y)],
  observed_AHCCD_temperature_C = pairs[variable == "temperature",
                                       .(station, year, y = observed)],
  observed_AHCCD_precip_log_mm = pairs[variable == "precipitation",
                                       .(station, year, y = log(observed))],
  ERA5_station_precip_log_mm = pairs[variable == "precipitation",
                                     .(station, year, y = log(value_era))],
  ERA5_minus_AHCCD_temperature_C = pairs[variable == "temperature",
                                         .(station, year, y = difference)],
  ERA5_minus_AHCCD_precip_log_ratio = pairs[variable == "precipitation",
                                            .(station, year, y = difference)]
)
pooled_sets <- lapply(pooled_sets, panel_cohort)
pooled_cov <- rbindlist(lapply(names(pooled_sets), function(nm) {
  d <- pooled_sets[[nm]]
  z <- attr(d, "coverage")
  if (is.null(z)) return(data.table())
  z[, set := nm]
  z
}), fill = TRUE)
write_table(pooled_cov, "pooled_station_eligibility.csv")

pooled_rows <- list()
for (name in names(pooled_sets)) {
  d <- pooled_sets[[name]]
  for (slope_type in c("common_trend", "site_specific_trends")) {
    varying <- identical(slope_type, "site_specific_trends")
    for (sp in c(NA_integer_, SPLITS)) {
      # Also bootstrap the site-specific-trend 1979 fit needed by margin_test.
      # Preserve the original common-trend random-number stream: use a separate
      # reproducible seed for the newly added site-specific bootstrap only.
      extra_boot <- varying && !is.na(sp) && sp == 1979L
      b_reps <- if (!is.na(sp) && (!varying || extra_boot)) PANEL_BOOT_REPS else 0L
      if (extra_boot) {
        saved_seed <- .Random.seed
        set.seed(PANEL_SEED + match(name, names(pooled_sets)))
      }
      r <- panel_fit(d, sp, site_trends = varying,
                     hac_lag = PANEL_HAC_LAG, boot_reps = b_reps)
      if (extra_boot) assign(".Random.seed", saved_seed, envir = globalenv())
      r[, `:=`(set = name, trend_structure = slope_type,
               split = if (is.na(sp)) NA_integer_ else sp,
               model = if (is.na(sp)) "trend_only" else "trend_step")]
      pooled_rows[[length(pooled_rows) + 1L]] <- r
    }
  }
}
pooled <- rbindlist(pooled_rows, fill = TRUE)
setcolorder(pooled, c("set", "trend_structure", "model", "split",
                      setdiff(names(pooled), c("set", "trend_structure", "model", "split"))))
write_table(pooled, "pooled_year_HAC_and_block_wild_models.csv")

# Deterministic 1979 temperature leave-one-station-out influence analysis.
# Retain the original eligible panel's station-years and refit without the
# omitted station. No bootstrap draws are taken, so the RNG stream is unchanged.
temp_panel <- pooled_sets[["ERA5_minus_AHCCD_temperature_C"]]
loo_ids <- unique(as.character(temp_panel$station))
loo_T <- if (length(loo_ids)) rbindlist(lapply(loo_ids, function(s) {
  d <- copy(temp_panel[station != s])
  d[, station_f := droplevels(station_f)]
  rbindlist(lapply(c("common_trend", "site_specific_trends"), function(tr) {
    r <- panel_fit(d, 1979L, site_trends = tr == "site_specific_trends",
                   boot_reps = 0L)
    r[, `:=`(dropped_station = s, trend_structure = tr, split = 1979L)]
    r
  }), fill = TRUE)
}), fill = TRUE) else data.table(
  dropped_station = character(), trend_structure = character(),
  split = integer(), converged = logical(), step = numeric(),
  ci95_lo = numeric(), ci95_hi = numeric())
write_table(loo_T, "ERA5_minus_AHCCD_temperature_1979_leave_one_station_out.csv")

# The conventional iid Gaussian AICc from an LM panel is *descriptive* only:
# cross-site dependence violates its likelihood assumption. It must not be
# used as an inferential Bayes-factor or independent confirmation.
pooled_cmp <- pooled[, .(
  aicc_trend = iid_aicc_descriptive[model == "trend_only"][1L],
  aicc_1977 = iid_aicc_descriptive[model == "trend_step" & split == 1977L][1L],
  aicc_1979 = iid_aicc_descriptive[model == "trend_step" & split == 1979L][1L]),
  by = .(set, trend_structure)]
pooled_cmp[, `:=`(aicc_1979_minus_1977 = aicc_1979 - aicc_1977,
                  aicc_1977_minus_trend = aicc_1977 - aicc_trend,
                  aicc_1979_minus_trend = aicc_1979 - aicc_trend)]
write_table(pooled_cmp, "pooled_iid_AICc_DESCRIPTIVE_ONLY.csv")

# The 21- and 55-mm benchmarks came from DIFFERENT earlier backcast models;
# neither is a known error magnitude. Convert a hypothetical fixed mm bias
# using each paired ERA5 station cell's own PRE-1979 Nov-Mar mean, not basin P.
# Project that station-specific, pre-1979 log perturbation onto the actual
# pooled regression design. This supplies a model-matched coefficient margin,
# allowing for unequal station coverage and optional station-specific trends.
# It remains a sensitivity approximation: annual precipitation is replaced by
# each site's pre-1979 mean in the nonlinear mm-to-log transformation.
precip_panel <- pooled_sets[["ERA5_minus_AHCCD_precip_log_ratio"]]
precip_pre <- merge(
  precip_panel[year < 1979L, .(station, year)],
  pairs[variable == "precipitation", .(station, year, value_era)],
  by = c("station", "year"), all.x = TRUE)
if (any(!is.finite(precip_pre$value_era)))
  stop("Missing ERA5 station-cell precipitation in the paired pre-1979 reference")
precip_station_ref <- precip_pre[, .(
  n_pre1979_winters = .N,
  pre1979_station_ERA5_NovMar_mm = mean(value_era)), by = station]
for (bias_mm in c(21, 55)) {
  nm <- paste0("hypothetical_", bias_mm, "mm_step_log_ratio")
  precip_station_ref[, (nm) := NA_real_]
  precip_station_ref[pre1979_station_ERA5_NovMar_mm > bias_mm,
                     (nm) := -log(pre1979_station_ERA5_NovMar_mm /
                                    (pre1979_station_ERA5_NovMar_mm - bias_mm))]
}
write_table(precip_station_ref, "pre1979_station_cell_precip_margin_QA.csv")

project_precip_margin <- function(bias_mm, trend_structure) {
  d <- merge(copy(precip_panel),
             precip_station_ref[, .(station, pre1979_station_ERA5_NovMar_mm)],
             by = "station", all.x = TRUE)
  p <- d$pre1979_station_ERA5_NovMar_mm
  if (!nrow(d) || any(!is.finite(p)) || any(p <= bias_mm)) return(NA_real_)
  d[, `:=`(step_indicator = as.integer(year >= 1979L),
           y = fifelse(year < 1979L, log(p / (p - bias_mm)), 0))]
  f <- if (trend_structure == "site_specific_trends")
    y ~ 0 + station_f + station_f:time + step_indicator else
      y ~ 0 + station_f + time + step_indicator
  mod <- tryCatch(lm(f, data = d, na.action = na.fail), error = function(e) NULL)
  if (is.null(mod) || is.na(coef(mod)["step_indicator"])) return(NA_real_)
  unname(coef(mod)["step_indicator"])
}
Ppre_basin <- annual[series == "Precip_NovMar" & year < 1979L,
                     mean(value, na.rm = TRUE)] # context only; NOT used in margins
margin_specs <- rbindlist(lapply(c("common_trend", "site_specific_trends"), function(tr) {
  data.table(
    trend_structure = tr,
    set = c(rep("ERA5_minus_AHCCD_precip_log_ratio", 2L),
            "ERA5_minus_AHCCD_temperature_C"),
    benchmark = c("hypothetical_21mm_wet_pre1979",
                  "hypothetical_55mm_wet_pre1979",
                  "hypothetical_0p5degC_cold_pre1979"),
    bias_mm = c(21, 55, NA_real_),
    margin = c(project_precip_margin(21, tr),
               project_precip_margin(55, tr), 0.5))
}), fill = TRUE)
margin_rows <- lapply(seq_len(nrow(margin_specs)), function(i) {
  spec <- margin_specs[i]
  z <- pooled[set == spec$set & trend_structure == spec$trend_structure &
                split == 1979L & model == "trend_step"]
  if (nrow(z) != 1L) stop("Missing pooled 1979 model: ", spec$set,
                          " / ", spec$trend_structure)
  cbind(spec, margin_test(z, spec$margin))
})
margin_evidence <- rbindlist(margin_rows, fill = TRUE)
margin_evidence[, pre1979_ERA5_basin_precip_context_mm := Ppre_basin]
write_table(margin_evidence, "directional_artifact_size_and_equivalence_sensitivity.csv")

# Placebo scan: descriptive AICc ranks only, with identical selected sites and
# observed rows within a set. No multiple-testing-adjusted change-point p value.
placebo_rows <- list()
for (name in c("observed_snowcourse_mm", "observed_AHCCD_precip_log_mm",
               "ERA5_minus_AHCCD_temperature_C",
               "ERA5_minus_AHCCD_precip_log_ratio")) {
  d <- pooled_sets[[name]]
  for (sp in PANEL_PLACEBO_YEARS) {
    z <- panel_fit(d, sp, site_trends = FALSE, boot_reps = 0L)
    placebo_rows[[length(placebo_rows) + 1L]] <- data.table(
      set = name, split = sp, n_years = z$n_years,
      step = z$step, p_year_hac_exploratory = z$p_year_hac,
      iid_aicc_descriptive = z$iid_aicc_descriptive)
  }
}
placebo <- rbindlist(placebo_rows, fill = TRUE)
placebo[, `:=`(
  delta_aicc_from_best = iid_aicc_descriptive - min(iid_aicc_descriptive,
                                                    na.rm = TRUE),
  rank_aicc_descriptive = frank(iid_aicc_descriptive, ties.method = "min",
                                na.last = "keep")), by = set]
write_table(placebo, "pooled_placebo_split_scan_DESCRIPTIVE.csv")

# 7. MONTHLY ERA5 vs discrete snow-course surveys: exploratory QA only.
# Survey visits nominally on Mar 1 and Apr 1 are NOT coincident with monthly
# ERA5 SWE means. A step in this pairwise difference cannot be attributed to
# ERA5 until ERA5 SWE is sampled on the actual survey dates.
snow_loc <- loc[station %in% unique(snow$station)]
if (nrow(snow_loc)) {
  sv <- extract_cells(swe_r, snow_loc$cell)
  uu <- unique(trimws(as.character(terra::units(swe_r))))
  uu <- uu[!is.na(uu) & nzchar(uu)]
  detected <- if (!length(uu)) "m" else
    if (all(grepl("^mm( |$)", uu, ignore.case = TRUE))) "mm" else
      if (all(grepl("^m( |$)", uu, ignore.case = TRUE))) "m" else "unknown"
  if (identical(detected, "unknown")) stop("Unrecognized ERA5 SWE units")
  if (!length(uu)) warning("SWE units absent: assuming metres as in Script 11")
  if (detected == "m") sv <- sv * 1000
  s_yr <- as.integer(format(dates, "%Y"))
  s_mon <- as.integer(format(dates, "%m"))
  tmp <- rbindlist(lapply(seq_len(nrow(snow_loc)), function(i) {
    x <- data.table(station = snow_loc$station[i], year = s_yr,
                    mon = s_mon, era_swe_mm = sv[i, ])[mon %in% c(3L, 4L)]
    x[, .(n_marapr_months = uniqueN(mon),
          era_monthly_marapr_mean_mm = if (.N == 2L && all(is.finite(era_swe_mm)))
            mean(era_swe_mm) else NA_real_), by = .(station, year)]
  }), fill = TRUE)
  snow_monthly_pairs <- merge(snow[, .(station, year, snowcourse_survey_mean_mm = swe)],
                              tmp, by = c("station", "year"))
  snow_monthly_pairs[, `:=`(
    difference_mm = era_monthly_marapr_mean_mm - snowcourse_survey_mean_mm,
    temporally_matched = FALSE)]
  write_table(snow_monthly_pairs,
              "ERA5_minus_snowcourse_monthly_vs_survey_QA_ONLY.csv")
}

# 8. Spatial reference check for temperature stations whose ERA5-minus-AHCCD
# 1979 step is BH-significant. Neighbours are selected ONCE by proximity and
# complete paired overlap; their composition is never allowed to change yearly.
# Stations BH-flagged for EITHER 1977 or 1979 are excluded from each focal
# neighbour set, to avoid using potentially discontinuous reference stations.
# AHCCD is reference-homogenized using neighbouring stations, so this check is
# not an independent test of station versus reanalysis inhomogeneity.
haversine_km <- function(lon1, lat1, lon2, lat2) {
  rad <- pi / 180
  dphi <- (lat2 - lat1) * rad
  dlambda <- (lon2 - lon1) * rad
  a <- sin(dphi / 2)^2 + cos(lat1 * rad) * cos(lat2 * rad) * sin(dlambda / 2)^2
  6371.0088 * 2 * atan2(sqrt(pmax(a, 0)), sqrt(pmax(1 - a, 0)))
}
sig_t <- if (nrow(site_diff)) unique(site_diff[
  variable == "temperature" & model == "trend_step_1979" &
    is.finite(p_bh) & p_bh < ALPHA, station]) else character()
flagged_t_all <- if (nrow(site_diff)) unique(site_diff[
  variable == "temperature" &
    model %in% c("trend_step_1977", "trend_step_1979") &
    is.finite(p_bh) & p_bh < ALPHA, station]) else character()
obsT <- pairs[variable == "temperature", .(station, year, observed)]
site_checks <- list()
for (s in sig_t) {
  focal <- ahccd_locs[station == s]
  if (nrow(focal) != 1L) next
  cand <- copy(ahccd_locs[!(station %in% flagged_t_all) & !coordinate_variation_flag &
                            is.finite(lon) & is.finite(lat)])
  cand[, km := haversine_km(focal$lon[1L], focal$lat[1L], lon, lat)]
  setorder(cand, km)
  cand <- cand[km <= 250]
  central <- obsT[station == s, .(year, focal = observed)]
  # Each accepted neighbour must preserve the COMMON intersection with all
  # previously accepted neighbours, for BOTH candidate years and n >= 16.
  wide <- copy(central)
  good <- character()
  for (nb in cand$station) {
    trial <- merge(wide, obsT[station == nb, .(year, v = observed)],
                   by = "year")
    ok <- nrow(trial) >= 16L && all(vapply(SPLITS, function(sp)
      min(sum(trial$year < sp), sum(trial$year >= sp)) >= 7L,
      logical(1L)))
    if (ok) {
      setnames(trial, "v", paste0("nb", length(good) + 1L))
      wide <- trial
      good <- c(good, nb)
    }
    if (length(good) == 3L) break
  }
  if (length(good) < 2L) {
    site_checks[[length(site_checks) + 1L]] <- data.table(
      station = s, n_neighbours = length(good),
      neighbours = paste(good, collapse = ";"),
      n_common_winters = nrow(wide),
      note = if (!nrow(cand))
        "No eligible neighbour candidates after QA, proximity and flagged-station exclusions" else
          "Coverage too short: fewer than two fixed neighbours preserve common winters for both 1977 and 1979")
    next
  }
  neighbor_cols <- paste0("nb", seq_along(good))
  wide[, station_minus_fixed_neighbours := focal - rowMeans(.SD),
       .SDcols = neighbor_cols]
  f <- fit_three(wide$year, wide$station_minus_fixed_neighbours,
                 "AHCCD_station_minus_fixed_neighbours", "temperature",
                 "observed_neighbor_difference", s, "site", 7L)
  z <- if (nrow(f)) f[model == "trend_step_1979"] else data.table()
  ref <- site_diff[variable == "temperature" &
                     model == "trend_step_1979" & station == s]
  has_steps <- nrow(z) == 1L && nrow(ref) == 1L &&
    is.finite(z$step[1L]) && is.finite(ref$step[1L])
  era_nb <- if (has_steps) ref$step[1L] + z$step[1L] else NA_real_
  # This identity is only approximate for fitted step coefficients: the
  # two regressions use different annual coverage and regression weights.
  note <- if (!has_steps) {
    "Model fit failed or produced a nonfinite 1979 step despite adequate common-year coverage"
  } else if (sign(z$step[1L]) == -sign(ref$step[1L]) &&
             abs(z$step[1L]) >= 0.5 * abs(ref$step[1L])) {
    "Substantial opposing station-neighbour step: station-relative shift plausible; AHCCD neighbour-reference homogenization limits independence"
  } else {
    "ERA5-relative shift persists against neighbours: ERA5-side or regional signal possible; AHCCD references limit independent attribution"
  }
  site_checks[[length(site_checks) + 1L]] <- data.table(
    station = s, n_neighbours = length(good),
    neighbours = paste(good, collapse = ";"),
    n_common_winters = nrow(wide),
    station_minus_neighbours_step_C = if (nrow(z)) z$step else NA_real_,
    neighbor_ci95_lo = if (nrow(z)) z$ci95_lo else NA_real_,
    neighbor_ci95_hi = if (nrow(z)) z$ci95_hi else NA_real_,
    neighbor_p_ar1 = if (nrow(z)) z$p_nominal else NA_real_,
    ERA5_minus_station_step_C = ref$step[1L],
    ERA5_minus_station_p_bh = ref$p_bh[1L],
    ERA5_minus_neighbours_step_C_approx = era_nb,
    note = note)
}
site_check <- if (length(site_checks)) rbindlist(site_checks, fill = TRUE) else
  data.table(station = character(), n_neighbours = integer(),
             note = character())
write_table(site_check, "AHCCD_1979_temperature_fixed_neighbour_checks.csv")


# 9. Preliminary descriptive results; final synthesis uses pooled year-HAC evidence.
# Nothing is labelled an 'artifact'
# based solely on a nonsignificant ERA5-minus-AHCCD test.
snow77 <- snow_fits[model == "trend_step_1977"]
snow79 <- snow_fits[model == "trend_step_1979"]
site_evidence <- merge(
  snow77[, .(station, step, p_nominal, p_bh)],
  snow_cmp[, .(station, aicc_1979_minus_1977, aicc_1977_minus_trend)],
  by = "station", all.x = TRUE)
site_evidence[, supports_1977 := is.finite(p_bh) & p_bh < ALPHA &
                is.finite(step) & step < 0 &
                is.finite(aicc_1979_minus_1977) & aicc_1979_minus_1977 >= 2 &
                is.finite(aicc_1977_minus_trend) & aicc_1977_minus_trend < 0]
write_table(site_evidence, "snowcourse_1977_evidence_sites.csv")

if (nrow(fixed_fits) && nrow(fixed_cmp)) {
  fixed_diff <- fixed_fits[response == "difference" &
                             model %in% c("trend_step_1977", "trend_step_1979")]
  fixed_diff <- merge(fixed_diff,
                      fixed_cmp[response == "difference",
                                .(variable, aicc_1979_minus_1977,
                                  aicc_1977_minus_trend, aicc_1979_minus_trend)],
                      by = "variable", all.x = TRUE)
} else {
  fixed_diff <- data.table(variable = character(), model = character(),
                           n = integer(), step = numeric(), ci95_lo = numeric(),
                           ci95_hi = numeric(), p_nominal = numeric(),
                           aicc_1977_minus_trend = numeric(),
                           aicc_1979_minus_trend = numeric())
}
write_table(fixed_diff, "AHCCD_fixed_network_difference_decision_table.csv")

fmt <- function(x, digits = 2L) ifelse(is.finite(x),
                                       formatC(x, digits = digits, format = "f"), "NA")

n_snow <- nrow(snow_cmp)
n_pref77 <- snow_cmp[is.finite(aicc_1979_minus_1977),
                     sum(aicc_1979_minus_1977 > 0)]
n_delta2 <- snow_cmp[is.finite(aicc_1979_minus_1977),
                     sum(aicc_1979_minus_1977 >= 2)]
n_snow_sig <- snow77[is.finite(p_bh), sum(p_bh < ALPHA & step < 0)]
n_snow79_sig <- snow79[is.finite(p_bh), sum(p_bh < ALPHA)]
n_both <- sum(site_evidence$supports_1977, na.rm = TRUE)
n_annual <- annual_cmp[is.finite(aicc_1979_minus_1977), .N]
n_annual_pref77 <- annual_cmp[is.finite(aicc_1979_minus_1977),
                              sum(aicc_1979_minus_1977 > 0)]
n_annual_delta2 <- annual_cmp[is.finite(aicc_1979_minus_1977),
                              sum(aicc_1979_minus_1977 >= 2)]
n_annual_step_vs_trend <- annual_cmp[
  is.finite(aicc_1979_minus_1977),
  sum((is.finite(aicc_1977_minus_trend) & aicc_1977_minus_trend <= -2) |
        (is.finite(aicc_1979_minus_trend) & aicc_1979_minus_trend <= -2))]
annual_trend_only_best <- annual_cmp[
  is.finite(aicc_1977_minus_trend) & is.finite(aicc_1979_minus_trend) &
    aicc_1977_minus_trend >= 0 & aicc_1979_minus_trend >= 0,
  variable]
swe_basin <- annual_cmp[variable == "SWE_MarApr"]
swe_date_delta <- if (nrow(swe_basin)) swe_basin$aicc_1979_minus_1977[1L] else NA_real_
swe_baseline_delta <- if (nrow(swe_basin)) swe_basin$aicc_1977_minus_trend[1L] else NA_real_

cat("\nSNOW COURSES: independent observed March-Apr SWE\n")
cat("Comparable sites: ", n_snow, "; AICc favors 1977 at ", n_pref77,
    "; delta AICc >= 2 at ", n_delta2, "; negative BH-significant 1977 steps: ",
    n_snow_sig, "; BH-significant 1979 steps: ", n_snow79_sig,
    "; sites meeting all prespecified 1977 support criteria: ", n_both, "\n", sep = "")
print(site_evidence[supports_1977 == TRUE, .(station, step, p_bh,
                                             aicc_1979_minus_1977, aicc_1977_minus_trend)])
if (nrow(snow_fixed_fits)) {
  cat("\nSNOW COURSE fixed-network standardized SWE (sensitivity check)\n")
  print(snow_fixed_fits[, .(model, n, aicc, step, ci95_lo, ci95_hi, p_nominal)])
} else {
  cat("No estimable fixed snow-course network; do not infer a fixed-network result.\n")
}

cat("\nERA5-LAND annual basin: timing is not attribution\n")
cat(n_annual_pref77, "/", n_annual,
    " series prefer 1977 OVER 1979 by AR(1)-GLS AICc; ",
    n_annual_step_vs_trend, "/", n_annual,
    " also improve on trend-only by at least 2 AICc units; ",
    n_annual_delta2, " have delta(1979-1977) >= 2. SWE_MarApr delta(1979-1977) = ",
    fmt(swe_date_delta), "; delta(1977-trend) = ", fmt(swe_baseline_delta), "\n", sep = "")
print(annual_cmp[, .(variable, aicc_1979_minus_1977, aicc_1977_minus_trend)])

cat("\nFIXED-NETWORK ERA5-minus-AHCCD differences (1968-1988)\n")
print(fixed_diff[, .(variable, model, n, step, ci95_lo, ci95_hi,
                     p_nominal, aicc_1977_minus_trend, aicc_1979_minus_trend)])
if (nrow(site_diff)) {
  cat("\nAHCCD individual-site differences (BH within variable and model)\n")
  print(site_diff[model %in% c("trend_step_1977", "trend_step_1979"),
                  .(n_tested = sum(is.finite(p_bh)),
                    n_bh_significant = sum(p_bh < ALPHA, na.rm = TRUE)),
                  by = .(variable, model)])
}
if (nrow(monthly_fits)) {
  cat("\nERA5 monthly exploratory BH-significant months\n")
  print(monthly_fits[model %in% c("trend_step_1977", "trend_step_1979") &
                       is.finite(p_bh) & p_bh < ALPHA,
                     .(months = paste(sort(as.integer(sub("month_", "", station))), collapse = ","),
                       n_months = .N), by = .(variable, model)])
}

# 10. DATA-CONDITIONAL INTERPRETATION. Do not equate failure to detect a step
# with physical homogeneity, and never treat a selected break-date p-value as
# if that date had been prespecified.
ref_not_confirmed <- nrow(fixed_diff) == 4L &&
  all(is.finite(fixed_diff$aicc_1977_minus_trend)) &&
  all(is.finite(fixed_diff$aicc_1979_minus_trend)) &&
  all(fixed_diff$aicc_1977_minus_trend > 0) &&
  all(fixed_diff$aicc_1979_minus_trend > 0)

pool_row <- function(s, y, trend = "common_trend") {
  pooled[set == s & split == y & trend_structure == trend &
           model == "trend_step"]
}
snow77_pool <- pool_row("observed_snowcourse_mm", 1977L)
snow79_pool <- pool_row("observed_snowcourse_mm", 1979L)
snow77_slopes <- pool_row("observed_snowcourse_mm", 1977L,
                          "site_specific_trends")
obsP77 <- pool_row("observed_AHCCD_precip_log_mm", 1977L)
obsP79 <- pool_row("observed_AHCCD_precip_log_mm", 1979L)
refP77 <- pool_row("ERA5_minus_AHCCD_precip_log_ratio", 1977L)
refP79 <- pool_row("ERA5_minus_AHCCD_precip_log_ratio", 1979L)
eraP77 <- pool_row("ERA5_station_precip_log_mm", 1977L)
refT77 <- pool_row("ERA5_minus_AHCCD_temperature_C", 1977L)
refT79 <- pool_row("ERA5_minus_AHCCD_temperature_C", 1979L)
refT79_slopes <- pool_row("ERA5_minus_AHCCD_temperature_C", 1979L,
                          "site_specific_trends")
snow_support <- n_both >= 1L && is.finite(swe_date_delta) &&
  swe_date_delta > 0
pool_snow77_negative <- nrow(snow77_pool) == 1L &&
  is.finite(snow77_pool$ci95_hi) && snow77_pool$ci95_hi < 0
pool_snow77_sensitivity <- nrow(snow77_slopes) == 1L &&
  is.finite(snow77_slopes$ci95_hi) && snow77_slopes$ci95_hi < 0

print_pooled <- pooled[trend_structure == "common_trend" &
                         model == "trend_step" & split %in% SPLITS,
                       .(set, split, n_station, n_years, step,
                         ci95_lo, ci95_hi, p_year_hac,
                         ci95_lo_block_wild, ci95_hi_block_wild)]
cat("\nPOOLED STATIONS: YEAR-HAC INFERENCE (sites share climate)\n")
print(print_pooled)
cat("\nHYPOTHETICAL ARTIFACT-SIZE EXCLUSION: DIRECTIONAL, NOT GENERAL EQUIVALENCE\n")
print(margin_evidence[, .(trend_structure, benchmark, bias_mm, artifact_margin,
                          p_exclude_directional_hac,
                          exclude_directional_95_hac,
                          exclude_directional_95_block_wild,
                          p_tost_symmetric_hac,
                          equivalent_within_pm_margin_hac)])
cat("\nTEMPERATURE STATION vs FIXED NEIGHBOURS (diagnostic, not proof)\n")
print(site_check)

fmt_pooled <- function(z, unit = "") {
  if (!nrow(z) || !is.finite(z$step) ||
      !is.finite(z$ci95_lo) || !is.finite(z$ci95_hi)) return("not estimable")
  sprintf("%+.3f%s (95%% year-HAC CI [%+.3f, %+.3f]; p=%.4f; %d stations, %d years)",
          z$step, unit, z$ci95_lo, z$ci95_hi, z$p_year_hac,
          z$n_station, z$n_years)
}
fmt_block_wild <- function(z) {
  if (nrow(z) != 1L || !is.finite(z$ci95_lo_block_wild) ||
      !is.finite(z$ci95_hi_block_wild)) return("95% block-wild interval unavailable")
  sprintf("95%% block-wild CI [%+.3f, %+.3f]",
          z$ci95_lo_block_wild, z$ci95_hi_block_wild)
}
margin_line <- function(name, slope_model = "common_trend") {
  r <- margin_evidence[benchmark == name & trend_structure == slope_model]
  if (nrow(r) != 1L) return(paste(name, "not evaluated"))
  if (!is.finite(r$artifact_margin))
    return(sprintf("%s (%s): not estimable%s.", name, slope_model,
                   if (is.na(r$bias_mm)) "" else
                     " (one or more station pre-1979 means do not exceed the specified mm bias)"))
  label <- if (is.na(r$bias_mm)) "specified hypothetical step" else
    "model-projected step from a specified mm bias"
  sprintf("%s (%s): %s %+.4f; one-sided year-HAC p=%.4f; year-HAC exclusion=%s; block-wild sensitivity exclusion=%s; symmetric TOST p=%.4f.",
          name, slope_model, label, r$artifact_margin, r$p_exclude_directional_hac,
          as.character(r$exclude_directional_95_hac),
          as.character(r$exclude_directional_95_block_wild),
          r$p_tost_symmetric_hac)
}

snow_statement <- if (pool_snow77_negative && snow_support) {
  if (pool_snow77_sensitivity)
    "Independent survey observations show a negative 1977 pooled step under both common-trend and site-specific-trend models; these and the basin SWE preference are consistent with a regional late-1970s decline." else
      "Independent survey observations show a negative 1977 pooled step under the common-trend model; the site-specific-trend sensitivity is inconclusive, so evidence of an abrupt common step is model-dependent."
} else if (snow_support) {
  "Several independent snow-course sites and basin SWE models favor 1977, but the pooled year-HAC 1977 effect does not meet the strict negative-interval criterion; evidence for a common abrupt step is mixed."
} else {
  "The independent and basin results do not jointly support a precise 1977 snowpack step under the declared criteria."
}

ref_strong <- margin_evidence[
  trend_structure == "common_trend" &
    benchmark %in% c("hypothetical_21mm_wet_pre1979",
                     "hypothetical_0p5degC_cold_pre1979")]
ref_exclusion <- nrow(ref_strong) == 2L &&
  all(ref_strong$exclude_directional_95_hac %in% TRUE) &&
  all(ref_strong$exclude_directional_95_block_wild %in% TRUE)
temperature_step_detected <- function(z) {
  nrow(z) == 1L && is.finite(z$step) && z$step < 0 &&
    is.finite(z$p_year_hac) && z$p_year_hac < ALPHA &&
    is.finite(z$ci95_hi) && z$ci95_hi < 0 &&
    is.finite(z$ci95_hi_block_wild) && z$ci95_hi_block_wild < 0
}
temperature_step_statement <- if (temperature_step_detected(refT79)) {
  paste0(if (temperature_step_detected(refT77))
    "Negative 1977 and 1979 steps are BOTH DETECTED in pooled ERA5-minus-AHCCD temperature " else
      "A negative 1979 step is DETECTED in pooled ERA5-minus-AHCCD temperature ",
    "under the station-FE/common-trend model (both year-HAC and 95% ",
    "block-wild intervals exclude zero). The 1979 sign is opposite to a ",
    "hypothesized pre-1979 ERA5 cold bias: the fitted ERA5 relative ",
    "temperature drops after 1978, so its inferred warming step is ",
    "smaller than the station-observed step. This does not locate the ",
    "source of the discrepancy or establish ERA5-only inhomogeneity.")
} else {
  "The 1979 pooled ERA5-minus-AHCCD temperature step is not confirmed as negative by both year-HAC and 95% block-wild criteria under the common-trend model."
}
ref_margin_statement <- if (ref_exclusion) {
  "For the full-record station-FE/common-trend specification, the illustrative +21-mm pre-1979 wet-bias and -0.5-C pre-1979 cold-bias scenarios are directionally excluded by both year-HAC and block-wild intervals. The precipitation benchmark uses station-cell pre-1979 means and the matched panel design; site-specific trends are evaluated separately. This is not general equivalence or proof of homogeneity."
} else if (ref_not_confirmed) {
  "The balanced 1968-1988 station-difference composites do not favor fixed steps over their trend-only baseline, but the full-record directional artifact-size exclusions are incomplete or model-dependent. Absence of a detected step is NOT proof of equivalence."
} else {
  "The reference comparisons do not justify excluding the prespecified hypothetical artifact magnitudes under both uncertainty calculations."
}
# Year-HAC-only LOO sensitivity. Missing/failed refits are labelled separately;
# the full-record block-wild interval must not be imputed to these refits.
loo_interpretation <- function(tr) {
  z <- loo_T[trend_structure == tr]
  if (!nrow(z)) return(sprintf("%s: no eligible omission fits", tr))
  failed <- z[converged != TRUE | !is.finite(step) |
                !is.finite(ci95_lo) | !is.finite(ci95_hi), dropped_station]
  lost <- z[converged == TRUE & is.finite(step) & is.finite(ci95_hi) &
              (step >= 0 | ci95_hi >= 0), dropped_station]
  if (length(failed)) return(sprintf(
    "%s: nonestimable after omitting %s%s", tr,
    paste(failed, collapse = ", "),
    if (length(lost)) paste0("; loses the negative 95% year-HAC exclusion after omitting ",
                             paste(lost, collapse = ", ")) else ""))
  if (length(lost)) return(sprintf(
    "%s: the negative 1979 step is NOT leave-one-station-out robust by 95%% year-HAC CI; omission of %s removes its exclusion of zero",
    tr, paste(lost, collapse = ", ")))
  sprintf("%s: the negative 1979 step retains a 95%% year-HAC interval below zero under all %d station omissions",
          tr, nrow(z))
}
loo_common_statement <- loo_interpretation("common_trend")
loo_slopes_statement <- loo_interpretation("site_specific_trends")

# Distinguish the independent observed 1977 precipitation step from an
# ERA5-minus-AHCCD 1979 reference diagnostic at a DIFFERENT break date.
negative_both_intervals <- function(z) nrow(z) == 1L &&
  is.finite(z$step) && z$step < 0 &&
  is.finite(z$ci95_hi) && z$ci95_hi < 0 &&
  is.finite(z$ci95_hi_block_wild) && z$ci95_hi_block_wild < 0
includes_zero_hac <- function(z) nrow(z) == 1L &&
  is.finite(z$ci95_lo) && is.finite(z$ci95_hi) &&
  z$ci95_lo <= 0 && z$ci95_hi >= 0
precip_statement <- if (negative_both_intervals(obsP77)) {
  observed_drop <- 100 * (1 - exp(obsP77$step[1L]))
  if (negative_both_intervals(eraP77) && includes_zero_hac(refP77)) {
    sprintf(paste0(
      "Independent AHCCD station observations show an approximately %.0f%% ",
      "Nov-Mar precipitation decline at the 1977 split (year-HAC and block-wild ",
      "95%% intervals below zero). Matched ERA5 station cells also have a ",
      "negative 1977 step, and the 1977 ERA5-minus-AHCCD year-HAC CI includes ",
      "zero. This corroborates the station-cell decline without proving ",
      "equivalence or a basin-wide step."), observed_drop)
  } else {
    sprintf(paste0(
      "Independent AHCCD station observations show an approximately %.0f%% ",
      "Nov-Mar precipitation decline at the 1977 split (year-HAC and block-wild ",
      "95%% intervals below zero). %s A 1979 ERA5-minus-AHCCD comparison alone ",
      "cannot verify ERA5 reproduction of the 1977 step; the observations are ",
      "from station cells, not a basin mean."), observed_drop,
      if (includes_zero_hac(refP79))
        "The 1979 ERA5-minus-AHCCD year-HAC interval includes zero."
      else "The 1979 ERA5-minus-AHCCD interval does not meet the zero-inclusion criterion.")
  }
} else {
  "The independent AHCCD 1977 precipitation step does not meet the criterion of negative year-HAC and block-wild 95% intervals; evidence for a robust common abrupt decline remains inconclusive."
}

# A 21-winter balanced composite can share the full-record sign without
# selecting a step by AR(1)-GLS AICc; the sample sizes/power differ.
fixed_t79 <- fixed_diff[variable == "temperature" & model == "trend_step_1979"]
fixed_reconciliation <- if (nrow(fixed_t79) == 1L &&
                            is.finite(fixed_t79$step) &&
                            is.finite(fixed_t79$ci95_lo) &&
                            is.finite(fixed_t79$ci95_hi)) {
  sprintf(paste0(
    "The fixed 1968-1988 temperature-difference network estimates a 1979 ",
    "step of %+.2f C (95%% AR(1)-GLS CI [%+.2f, %+.2f], %d winters). ",
    "Its sign is %s the full-record estimate; its shorter record provides ",
    "less information and a step need not beat trend-only by AICc. ",
    "This is not a contradiction between the two analyses."),
    fixed_t79$step, fixed_t79$ci95_lo, fixed_t79$ci95_hi,
    fixed_t79$n,
    if (nrow(refT79) == 1L && is.finite(refT79$step) &&
        sign(fixed_t79$step) == sign(refT79$step)) "the same as" else
          "different from")
} else {
  "The fixed 1968-1988 temperature network is not estimable; it cannot be used to confirm or contradict the full-record pooled result."
}

ref_statement <- paste(temperature_step_statement, loo_common_statement,
                       ref_margin_statement)

summary_text <- c(
  "# Nechako: 1976/77 climate-shift and 1979 ERA5-Land artifact evidence",
  "",
  "## Interpretation",
  snow_statement,
  precip_statement,
  ref_statement,
  fixed_reconciliation,
  "No result establishes that the PDO caused the change or rules out an additional ERA5 forcing inhomogeneity.",
  "",
  "## Independent BC snow courses",
  sprintf("- Site-level AR(1): %d/%d stations prefer 1977 to 1979 by AICc; %d have a preference of at least 2 AICc units; %d sites pass the negative-1977/BH/AICc joint criterion.",
          n_pref77, n_snow, n_delta2, n_both),
  paste0("- Full-record pooled observed SWE, step after 1976: ",
         fmt_pooled(snow77_pool, " mm SWE"), "."),
  paste0("- Full-record pooled observed SWE, step after 1978: ",
         fmt_pooled(snow79_pool, " mm SWE"), "."),
  paste0("- Site-specific slope sensitivity, step after 1976: ",
         fmt_pooled(snow77_slopes, " mm SWE"), "."),
  "- Fixed 1968-1988 network sensitivity is recorded in snowcourse_fixed_network_models.csv.",
  "",
  "## Independent AHCCD precipitation and ERA5-relative differences",
  paste0("- Pooled observed AHCCD Nov-Mar log precipitation, 1977 step: ",
         fmt_pooled(obsP77, " log(mm)"), "; ", fmt_block_wild(obsP77), "."),
  paste0("- Pooled observed AHCCD Nov-Mar log precipitation, 1979 step: ",
         fmt_pooled(obsP79, " log(mm)"), "."),
  paste0("- Matched ERA5 station-cell Nov-Mar log precipitation, 1977 step: ",
         fmt_pooled(eraP77), "; ", fmt_block_wild(eraP77), "."),
  paste0("- ERA5 minus AHCCD Nov-Mar log precipitation ratio, 1977 step: ",
         fmt_pooled(refP77), "; ", fmt_block_wild(refP77), "."),
  paste0("- ERA5 minus AHCCD Nov-Mar log precipitation ratio, 1979 step: ",
         fmt_pooled(refP79), "; ", fmt_block_wild(refP79), "."),
  paste0("- ERA5 minus AHCCD Nov-Mar temperature, 1977 step: ",
         fmt_pooled(refT77, " C"), "; ", fmt_block_wild(refT77), "."),
  paste0("- ERA5 minus AHCCD Nov-Mar temperature, 1979 step: ",
         fmt_pooled(refT79, " C"), "; ", fmt_block_wild(refT79), "."),
  paste0("- 1979 temperature step with site-specific station slopes: ",
         fmt_pooled(refT79_slopes, " C"), "; ",
         fmt_block_wild(refT79_slopes), "."),
  sprintf("- AHCCD comparison candidates include %d ERA5-domain station locations (%d marked inside the native Nechako basin mask); %d temperature stations and %d precipitation stations enter the respective eligible pooled panels. Snow courses are restricted to the basin.",
          nrow(ahccd_locs), sum(ahccd_locs$in_basin, na.rm = TRUE),
          uniqueN(pooled_sets[["ERA5_minus_AHCCD_temperature_C"]]$station),
          uniqueN(pooled_sets[["ERA5_minus_AHCCD_precip_log_ratio"]]$station)),
  paste0("- 1979 ERA5-minus-AHCCD temperature leave-one-station-out (year-HAC only): ",
         loo_common_statement, ". ", loo_slopes_statement,
         ". All refits appear in ERA5_minus_AHCCD_temperature_1979_leave_one_station_out.csv."),
  sprintf("- Reference fixed-network composites: %s.",
          if (ref_not_confirmed)
            "neither 1977 nor 1979 difference step improves AR(1)-GLS AICc over trend-only" else
              "mixed or inconclusive model comparisons"),
  "",
  "## Directional hypothetical artifact sizes (not pre-registered)",
  sprintf("- Pre-1979 ERA5 basin Nov-Mar P mean = %s mm (context only; NOT used to convert station-level precipitation effects). Station-cell pre-1979 means and panel-design projections are in pre1979_station_cell_precip_margin_QA.csv and the directional sensitivity CSV.",
          fmt(Ppre_basin)),
  margin_line("hypothetical_21mm_wet_pre1979"),
  margin_line("hypothetical_55mm_wet_pre1979"),
  margin_line("hypothetical_0p5degC_cold_pre1979"),
  margin_line("hypothetical_21mm_wet_pre1979", "site_specific_trends"),
  margin_line("hypothetical_55mm_wet_pre1979", "site_specific_trends"),
  margin_line("hypothetical_0p5degC_cold_pre1979", "site_specific_trends"),
  "- A one-sided test rejecting a negative (or positive) hypothesized step is not a two-sided equivalence test. Two-sided TOST uses the 90% CI within symmetric +/- margins and is reported separately.",
  "- The margins (+21/+55 mm, 0.5 C) are sensitivity benchmarks drawn from alternative models or judgment, not validated physical bias thresholds. A mm bias maps nonlinearly to each site's log ratio; the model projection uses a fixed pre-1979 station mean. Conclusions are specific to the period, station sample and regression specification.",
  "",
  "## ERA5-Land basin, seasonality and independent-site review",
  sprintf("- SWE Mar-Apr basin AICc(1979)-AICc(1977) = %s; AICc(1977)-AICc(trend only) = %s (AR1 basin models).",
          fmt(swe_date_delta), fmt(swe_baseline_delta)),
  sprintf("- Across annual basin variables (excluding PET_ERA5_annual), %d/%d prefer 1977 OVER 1979 by AR(1)-GLS AICc, %d by at least 2 units; %d/%d improve on trend-only by at least 2 AICc units with at least one step. These are correlated ERA5-derived variables, not independent experiments.",
          n_annual_pref77, n_annual, n_annual_delta2, n_annual_step_vs_trend, n_annual),
  sprintf("- Trend-only has the lowest AR(1)-GLS AICc for: %s.",
          if (length(annual_trend_only_best)) paste(annual_trend_only_best, collapse = ", ") else "none of the evaluated series"),
  sprintf("- Seven-series pooled panel analysis includes station FE, optional station-specific trends and year-level HAC (Bartlett lag %d). A synchronized %d-year block-wild score sensitivity uses %d replications and seed %d.",
          PANEL_HAC_LAG, PANEL_BLOCK_LENGTH, PANEL_BOOT_REPS, PANEL_SEED),
  sprintf("- Temperature stations with BH-significant ERA5-AHCCD 1979 differences: %d; fixed-neighbour checks in AHCCD_1979_temperature_fixed_neighbour_checks.csv.",
          length(sig_t)),
  "- The 1960-1995 alternative-date scan ranks *descriptive* iid-panel AICc; multiple date searches invalidate treating the smallest unadjusted p-value as confirmatory.",
  "- The file ERA5_minus_snowcourse_monthly_vs_survey_QA_ONLY.csv is NOT a survey-date-matched difference test; it cannot exclude a 1979 SWE artifact. Daily ERA5 SWE on actual visit dates is needed.",
  "- Monthly SWE zero/near-zero QA is saved separately; summer changes may reflect snow presence instead of SWE magnitude.",
  "",
  "## Statistical and data limitations",
  "1. Station-level intercepts remove fixed station-level biases, but unequal time coverage and site-specific trends can still affect pooled estimates. Therefore site-specific slopes and the balanced-network results must also be checked.",
  "2. Regional climatic forcing induces cross-station covariance; year-level HAC and block-wild scores do not render the observed years fully independent and may be sensitive to lag/block length.",
  "3. The pooled iid Gaussian AICc is descriptive under cross-station correlation. Only within-series, identical-row comparisons are displayed; no significance claim is based on its magnitude alone.",
  "4. Fixed-neighbour steps are weak attribution evidence: AHCCD homogenization draws on neighbouring reference stations, and the ERA5-minus-neighbour step reported here adds coefficients from unequal year sets (approximation only). Other BH-flagged 1977 OR 1979 stations are excluded from the neighbour pool. Neighbours are retained only if their joint common-year intersection meets both split-date coverage requirements. Monthly AHCCD coordinate metadata cannot reliably identify historical relocation.",
  "5. Survey dates may deviate by 10 days from nominal Mar/Apr visits, point/grid spatial representativeness differs, and 1977 versus 1979 is intrinsically hard to distinguish from annual data.",
  "6. AHCCD and ERA5 are not fully independent if the same station observations contributed to reanalysis forcing; failure to reject a station-reference difference is not proof of agreement. Leave-one-out refits have year-HAC inference only (no block-wild bootstrap), and the snow-course comparison is not independently date-matched.",
  "",
  sprintf("Observations: %s. Cached observation files were not overwritten.",
          if (length(OBS_DOWNLOADED))
            paste("Downloaded from BC/ECCC:", paste(OBS_DOWNLOADED, collapse = ", ")) else
              "All three public-data snapshots loaded from local cache"),
  "Original Script 11 and other project scripts were not edited."
)
writeLines(summary_text, file.path(OUT, "evidence_summary.md"), useBytes = TRUE)


# Provenance manifest: original observation snapshots (MD5 checksums) and
# derived ERA5 products (path, size, modification time). No file is uploaded.
source_paths <- c(
  find_observation("allmss_archive.csv"),
  find_observation("snow_course_locations.csv"),
  find_observation("ahccd_monthly_bbox.csv"),
  pr_path, t_path, swe_path,
  file.path(BASIN, "swe_homogeneity_annual_series.csv"),
  file.path(ROOT, "trend_analysis_pr_pet", "regional", "regional_avg_monthly.csv"),
  file.path(ROOT, "trend_analysis_pr_pet", "native_basin_mask.rds"))
source_type <- c(rep("independent_observation_snapshot", 3L),
                 rep("Monthly_ERALand_download.R", 3L),
                 rep("11pr_pet_trends.r", 3L))
fi <- file.info(source_paths)
manifest <- data.table(source_type = source_type,
                       file = normalizePath(source_paths, winslash = "/", mustWork = TRUE),
                       bytes = as.numeric(fi$size),
                       modified_utc = format(fi$mtime, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
                       md5 = c(unname(tools::md5sum(source_paths[1:3])),
                               rep(NA_character_, length(source_paths) - 3L)))
write_table(manifest, "input_provenance_manifest.csv")
write_table(data.table::data.table(
  file = names(OBS_SOURCES),
  source_url = unname(OBS_SOURCES),
  retrieval_status = ifelse(names(OBS_SOURCES) %in% OBS_DOWNLOADED,
                            "downloaded_this_run", "existing_local_cache"),
  retrieved_at_utc = ifelse(names(OBS_SOURCES) %in% OBS_DOWNLOADED,
                            format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), NA_character_)),
  "observation_source_provenance.csv")
cat("\n", paste(summary_text, collapse = "\n"), "\n", sep = "")
cat("\nOutputs: ", OUT, "\n", sep = "")

# Complexity: preprocessing O(M log M + S log S) for M AHCCD months and S
# snow-course visits; pooled OLS O(N K^2 + K^3) for N station-years and K
# coefficients, year-HAC O(N K + H Y K^2), and block-wild sensitivity
# O(B (Y K + K^2)) per fitted model for H lags, Y years and B replications.
# The added two slope-structure block-wild margins and station-cell projections
# require O(N K^2 + K^3) time per projection and O(N + K^2) working memory;
# added site-specific 1979 bootstraps use O(B (Y K + K^2)) per station set.
# The LOO sensitivity adds O(J (N K^2 + K^3 + H Y K^2)) time for J station
# omissions (two slope specifications), O(N + Y K + K^2) working memory;
# no RNG use. Fixed-neighbour selection takes O(C N log N) worst-case time
# for C candidate sites and N jointly matching annual observations.
# Memory O(M + S + N K + Y K + K^2). Original site/basin GLS fits are
# unchanged (conservative O(G * I * n^3) time, O(n^2) working memory). No
# imputation, no PDO inference and no causal attribution from placebo ranks.
