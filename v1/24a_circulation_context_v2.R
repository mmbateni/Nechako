################################################################################
# 14a_circulation_context_v2.R
# LARGE-SCALE CIRCULATION CONTEXT OF THE 2022-2025 DROUGHT
# CURRENT MANUSCRIPT RESULTS (Methods/Results 3g/4g)
#
# Inputs produced by Monthly_ERALand_download_v2.R with
# RUN_CIRCULATION_CONTEXT <- TRUE:
#   monthly_data_direct/geopotential_500hPa_monthly.nc
#   monthly_data_direct/ivt_components_monthly.nc
#
# Diagnostics
#   - Z500 calendar-month anomalies [m]
#   - IVT vector-component anomalies
#   - IVT-magnitude anomaly
#   - composites for principal drought episodes
#   - 2023-2024 El Nino-period context composite
#
# IMPORTANT
#   This script provides circulation context, not causal attribution.
#   It stops if the ERA5 circulation files are absent.
#
# MEMORY / PERFORMANCE ARCHITECTURE
#   - No zero-argument terra::mem_info() call.
#   - terra::mem_info(x) is called only after a SpatRaster is opened.
#   - The complete circulation fields are small enough to fit in RAM
#     (~0.31 GB per 912-layer field on the current 161 x 281 grid).
#   - Each circulation field is therefore materialized ONCE with toMemory().
#   - Calendar-month climatology uses terra::tapp() over the in-memory field.
#   - This avoids 12 repeated disk-backed NetCDF reductions per variable.
#   - Intermediate climatology objects are released after each anomaly field.
#
# OUTPUTS -> circulation_context/
#   circulation_domain_monthly_summary.csv
#   circulation_event_composites.csv
#   circulation_event_composite_maps.nc
#   Fig_circulation_event_composites.pdf
#   Fig_circulation_event_composites.png
################################################################################


# ==============================================================================
# 0. WORKING DIRECTORY / SHARED UTILITIES
# ==============================================================================

setwd(
  Sys.getenv(
    "NECHAKO_WD",
    "D:/Nechako_Drought/Nechako/"
  )
)

source("DROUGHT_ANALYSIS_utils_v2.R")

utils_load_packages(
  c(
    "terra",
    "dplyr",
    "data.table",
    "lubridate",
    "ggplot2",
    "patchwork",
    "tidyr",
    "sf"
  )
)

if (!dir.exists(WD_PATH)) {
  stop(
    "Working directory not found: ",
    WD_PATH
  )
}

setwd(WD_PATH)


# ==============================================================================
# 1. OUTPUT / INPUT PATHS
# ==============================================================================

OUT_DIR <- file.path(
  WD_PATH,
  "circulation_context"
)

dir.create(
  OUT_DIR,
  showWarnings = FALSE,
  recursive = TRUE
)

Z_FILE <- file.path(
  WD_PATH,
  "monthly_data_direct",
  "geopotential_500hPa_monthly.nc"
)

IVT_FILE <- file.path(
  WD_PATH,
  "monthly_data_direct",
  "ivt_components_monthly.nc"
)

missing_inputs <- c(
  Z_FILE,
  IVT_FILE
)[
  !file.exists(
    c(
      Z_FILE,
      IVT_FILE
    )
  )
]

if (length(missing_inputs)) {
  stop(
    "Circulation-context input(s) missing:\n  ",
    paste(
      missing_inputs,
      collapse = "\n  "
    ),
    "\nSet RUN_CIRCULATION_CONTEXT <- TRUE in ",
    "Monthly_ERALand_download_v2.R and run its Phase 3 ",
    "before running this script."
  )
}


# ==============================================================================
# 2. AUTHORITATIVE ANALYSIS PERIOD
# ==============================================================================

ANALYSIS_START <- as.Date(
  "1950-01-01"
)

ANALYSIS_END <- as.Date(
  "2025-12-01"
)

EXPECTED_MONTHS <- seq(
  ANALYSIS_START,
  ANALYSIS_END,
  by = "month"
)

EXPECTED_N_MONTHS <- length(
  EXPECTED_MONTHS
)

if (EXPECTED_N_MONTHS != 912L) {
  stop(
    "Internal analysis-period error: expected 912 months, found ",
    EXPECTED_N_MONTHS
  )
}


# ==============================================================================
# 3. CONSOLE PROGRESS HELPERS
# ==============================================================================

stage_start <- function(label) {
  
  cat(
    sprintf(
      "[%s] %-65s",
      format(
        Sys.time(),
        "%H:%M:%S"
      ),
      paste0(
        label,
        "..."
      )
    )
  )
  
  flush.console()
  
  Sys.time()
}


stage_done <- function(t0) {
  
  elapsed <- as.numeric(
    Sys.time() - t0,
    units = "secs"
  )
  
  cat(
    sprintf(
      " done (%.1f s)\n",
      elapsed
    )
  )
  
  flush.console()
  
  invisible(elapsed)
}


stage_message <- function(...) {
  cat(
    sprintf(
      "[%s] ",
      format(Sys.time(), "%H:%M:%S")
    ),
    ...,
    "\n",
    sep = ""
  )
  flush.console()
  invisible(NULL)
}

# Line-based progress for operations that otherwise appear silent for several
# minutes (especially NetCDF reads and large raster reductions).
progress_step <- function(current, total, label, start_time, step = 0.05) {
  if (total <= 0L) return(invisible(NULL))
  pct <- current / total
  if (current == 1L || current == total || pct >= progress_step.last + step) {
    elapsed <- as.numeric(Sys.time() - start_time, units = "secs")
    rate <- if (elapsed > 0) current / elapsed else NA_real_
    eta <- if (is.finite(rate) && rate > 0) (total - current) / rate else NA_real_
    cat(sprintf(
      "[%s] %-34s %5.1f%% (%d/%d) | elapsed %s | ETA %s\n",
      format(Sys.time(), "%H:%M:%S"),
      label,
      100 * pct, current, total,
      format_elapsed_local(elapsed),
      format_elapsed_local(eta)
    ))
    flush.console()
    progress_step.last <<- pct
  }
  invisible(NULL)
}

progress_step.last <- -Inf

format_elapsed_local <- function(seconds) {
  if (!is.finite(seconds)) return("--")
  if (seconds < 60) return(sprintf("%.0fs", seconds))
  sprintf("%dm %.0fs", floor(seconds / 60), seconds %% 60)
}

materialize_with_progress <- function(x, label, target_fraction = 0.05) {
  nr <- terra::nrow(x)
  nc <- terra::ncol(x)
  nl <- terra::nlyr(x)
  if (nr <= 0L || nc <= 0L || nl <= 0L) stop("Cannot materialize an empty raster: ", label)
  
  # Allocate the complete values matrix once, then fill it row blocks. This
  # preserves the existing in-memory architecture but makes the long NetCDF
  # read observable in the console.
  vals <- matrix(NA_real_, nrow = nr * nc, ncol = nl)
  rows_per_chunk <- max(1L, ceiling(nr * target_fraction))
  row_starts <- seq.int(1L, nr, by = rows_per_chunk)
  total_chunks <- length(row_starts)
  
  t0 <- Sys.time()
  progress_step.last <<- -Inf
  terra::readStart(x)
  on.exit(terra::readStop(x), add = TRUE)
  
  for (j in seq_along(row_starts)) {
    r0 <- row_starts[j]
    nrows <- min(rows_per_chunk, nr - r0 + 1L)
    block <- terra::readValues(x, row = r0, nrows = nrows, mat = TRUE)
    if (is.null(dim(block))) block <- matrix(block, ncol = nl)
    idx <- ((r0 - 1L) * nc + 1L):( (r0 - 1L + nrows) * nc )
    vals[idx, ] <- block
    progress_step(j, total_chunks, paste0("Reading ", label), t0)
  }
  
  terra::readStop(x)
  
  out <- terra::rast(x)
  terra::values(out) <- vals
  terra::time(out) <- terra::time(x)
  rm(vals)
  gc(verbose = FALSE)
  out
}


# ==============================================================================
# 4. TIME-AXIS READER
# ==============================================================================

read_dates <- function(
    x,
    label) {
  
  tt <- terra::time(x)
  
  if (
    is.null(tt) ||
    length(tt) != terra::nlyr(x) ||
    anyNA(tt)
  ) {
    stop(
      "Invalid NetCDF time axis: ",
      label
    )
  }
  
  dates <- as.Date(tt)
  
  # Normalize all monthly timestamps to YYYY-MM-01.
  dates <- as.Date(
    format(
      dates,
      "%Y-%m-01"
    )
  )
  
  if (
    length(dates) != terra::nlyr(x) ||
    anyNA(dates)
  ) {
    stop(
      "Failed to normalize NetCDF monthly time axis: ",
      label
    )
  }
  
  dates
}


# ==============================================================================
# 5. MONTHLY RECORD VALIDATION
# ==============================================================================

validate_circulation_dates <- function(
    dates,
    label,
    start_date = ANALYSIS_START,
    end_date = ANALYSIS_END) {
  
  dates <- as.Date(dates)
  
  if (!length(dates)) {
    stop(
      label,
      " contains no dates."
    )
  }
  
  if (anyNA(dates)) {
    stop(
      label,
      " contains NA dates."
    )
  }
  
  if (anyDuplicated(dates)) {
    stop(
      label,
      " contains duplicated months."
    )
  }
  
  if (is.unsorted(dates)) {
    stop(
      label,
      " is not chronologically sorted."
    )
  }
  
  expected <- seq(
    as.Date(start_date),
    as.Date(end_date),
    by = "month"
  )
  
  if (!identical(dates, expected)) {
    
    missing <- setdiff(
      expected,
      dates
    )
    
    extra <- setdiff(
      dates,
      expected
    )
    
    msg <- paste0(
      label,
      " does not exactly match the required monthly record."
    )
    
    if (length(missing)) {
      msg <- paste0(
        msg,
        "\nMissing month(s): ",
        paste(
          missing,
          collapse = ", "
        )
      )
    }
    
    if (length(extra)) {
      msg <- paste0(
        msg,
        "\nUnexpected month(s): ",
        paste(
          extra,
          collapse = ", "
        )
      )
    }
    
    stop(msg)
  }
  
  invisible(TRUE)
}


# ==============================================================================
# 6. CALENDAR-MONTH ANOMALY CALCULATION
#
# Definition:
#
#   anomaly(t) =
#       X(t) - climatology(month(t))
#
# where climatology(month) is calculated independently for each calendar
# month over the complete 1950-2025 record.
#
# PERFORMANCE:
#   The entire field is first materialized once.
#   terra::tapp() then computes all 12 monthly climatologies in one operation.
#
# This is deliberately preferable to:
#
#   12 x terra::mean(x[[ii]])
#
# on a disk-backed NetCDF, which caused the observed ~437-second Z500
# climatology bottleneck.
# ==============================================================================

monthly_anomaly_raster <- function(
    x,
    dates,
    label = "") {
  
  dates <- as.Date(dates)
  
  if (
    length(dates) != terra::nlyr(x)
  ) {
    stop(
      "Date/layer count mismatch in monthly_anomaly_raster(): ",
      label
    )
  }
  
  months <- as.integer(
    format(
      dates,
      "%m"
    )
  )
  
  if (
    anyNA(months) ||
    any(!months %in% 1:12)
  ) {
    stop(
      "Invalid calendar-month values in monthly_anomaly_raster(): ",
      label
    )
  }
  
  tag <- if (nzchar(label)) {
    paste0(
      " [",
      label,
      "]"
    )
  } else {
    ""
  }
  
  
  # --------------------------------------------------------------------------
  # Memory diagnostic
  # --------------------------------------------------------------------------
  
  stage_message(
    "Checking raster memory requirement",
    tag
  )
  
  terra::mem_info(x)
  
  
  # --------------------------------------------------------------------------
  # Materialize the complete field ONCE.
  #
  # Current circulation grid:
  #   161 x 281 x 912
  #   ~0.31 GB according to terra::mem_info()
  #
  # The current machine reports:
  #   available = 16 GB
  #   terra allowed = 8 GB
  #
  # Therefore the field fits comfortably in memory.
  # --------------------------------------------------------------------------
  
  t0 <- stage_start(
    paste0(
      "Materializing circulation field into RAM",
      tag
    )
  )
  
  x_mem <- materialize_with_progress(
    x,
    label = if (nzchar(label)) label else "circulation field"
  )
  
  stage_done(t0)
  
  
  # --------------------------------------------------------------------------
  # Calendar-month climatology.
  #
  # tapp() receives the complete 912-layer in-memory field and groups layers
  # by calendar month. This performs 12 reductions without repeatedly opening
  # and scanning the NetCDF source.
  # --------------------------------------------------------------------------
  
  t0 <- stage_start(
    paste0(
      "Computing calendar-month climatology (tapp)",
      tag
    )
  )
  
  clim <- terra::tapp(
    x_mem,
    index = months,
    fun = "mean",
    na.rm = TRUE,
    cores = max(
      1L,
      parallel::detectCores() - 1L
    )
  )
  
  names(clim) <- sprintf(
    "%02d",
    1:12
  )
  
  stage_done(t0)
  
  
  # --------------------------------------------------------------------------
  # Expand 12 climatological layers to the full 912-month sequence.
  # --------------------------------------------------------------------------
  
  t0 <- stage_start(
    paste0(
      "Constructing time-varying climatology",
      tag
    )
  )
  
  clim_for_time <- clim[[months]]
  
  stage_done(t0)
  
  
  # --------------------------------------------------------------------------
  # Calculate anomalies.
  # --------------------------------------------------------------------------
  
  t0 <- stage_start(
    paste0(
      "Subtracting climatology",
      tag
    )
  )
  
  out <- x_mem -
    clim_for_time
  
  terra::time(out) <- terra::time(
    x_mem
  )
  
  stage_done(t0)
  
  
  # --------------------------------------------------------------------------
  # Release intermediate objects.
  # --------------------------------------------------------------------------
  
  rm(
    x_mem,
    clim,
    clim_for_time
  )
  
  gc(
    verbose = FALSE
  )
  
  out
}


# ==============================================================================
# 7. IVT COMPONENT EXTRACTION
#
# Preferred path:
#   terra::sds()
#
# Fallback:
#   terra::rast()
#
# The returned SpatRasters remain disk-backed until they are passed to
# monthly_anomaly_raster(), which materializes one complete field at a time.
# ==============================================================================

extract_two_ivt_components <- function(
    f) {
  
  s <- try(
    terra::sds(f),
    silent = TRUE
  )
  
  if (
    !inherits(s, "try-error") &&
    length(s) >= 2L
  ) {
    
    nm <- tolower(
      names(s)
    )
    
    east <- grep(
      "eastward|east|u",
      nm
    )
    
    north <- grep(
      "northward|north|v",
      nm
    )
    
    if (
      length(east) &&
      length(north)
    ) {
      
      return(
        list(
          east = s[[east[1]]],
          north = s[[north[1]]]
        )
      )
    }
    
    # Final subdataset fallback: retain downloader ordering.
    return(
      list(
        east = s[[1]],
        north = s[[2]]
      )
    )
  }
  
  
  # --------------------------------------------------------------------------
  # Fallback: ordinary raster interface.
  # --------------------------------------------------------------------------
  
  x <- terra::rast(
    f
  )
  
  nm <- tolower(
    names(x)
  )
  
  east <- grep(
    "eastward|east",
    nm
  )
  
  north <- grep(
    "northward|north",
    nm
  )
  
  if (
    !length(east) ||
    !length(north)
  ) {
    stop(
      "Could not identify eastward/northward IVT components in ",
      f
    )
  }
  
  list(
    east = x[[east[1]]],
    north = x[[north[1]]]
  )
}


# ==============================================================================
# 8. GDAL / TERRA PERFORMANCE OPTIONS
# ==============================================================================

stage_message(
  "Configuring terra/GDAL performance options"
)

terra::gdalCache(
  8192
)

terra::terraOptions(
  progress = 1
)


# ==============================================================================
# 9. OPEN Z500
# ==============================================================================

t0 <- stage_start(
  "Opening Z500"
)

z500 <- terra::rast(
  Z_FILE
)

z_time <- terra::time(
  z500
)

stage_done(t0)


stage_message(
  "Z500 memory diagnostic"
)

# IMPORTANT:
# terra::mem_info() without an argument is invalid in the installed terra
# version. A SpatRaster must be supplied.
terra::mem_info(
  z500
)


# ==============================================================================
# 10. READ / VALIDATE Z500 TIME AXIS
# ==============================================================================

t0 <- stage_start(
  "Reading and validating Z500 monthly record"
)

z_dates <- read_dates(
  z500,
  "Z500"
)

validate_circulation_dates(
  z_dates,
  label = "Z500 monthly record"
)

terra::time(
  z500
) <- z_time

stage_done(t0)

stage_message(
  "Z500:",
  terra::nrow(z500),
  "rows ×",
  terra::ncol(z500),
  "columns ×",
  terra::nlyr(z500),
  "months"
)


# ==============================================================================
# 11. OPEN IVT COMPONENTS
# ==============================================================================

t0 <- stage_start(
  "Opening eastward/northward IVT components"
)

ivt <- extract_two_ivt_components(
  IVT_FILE
)

stage_done(t0)

uivt <- ivt$east
vivt <- ivt$north

rm(ivt)

u_time <- terra::time(
  uivt
)

v_time <- terra::time(
  vivt
)


# ==============================================================================
# 12. READ / VALIDATE IVT TIME AXES
# ==============================================================================

t0 <- stage_start(
  "Reading and validating IVT monthly records"
)

u_dates <- read_dates(
  uivt,
  "eastward IVT"
)

v_dates <- read_dates(
  vivt,
  "northward IVT"
)

validate_circulation_dates(
  u_dates,
  label = "eastward IVT monthly record"
)

validate_circulation_dates(
  v_dates,
  label = "northward IVT monthly record"
)

if (!identical(
  u_dates,
  v_dates
)) {
  stop(
    "Eastward and northward IVT time axes differ."
  )
}

if (!identical(
  z_dates,
  u_dates
)) {
  stop(
    "Z500 and IVT time axes differ."
  )
}

if (
  terra::nlyr(uivt) != EXPECTED_N_MONTHS ||
  terra::nlyr(vivt) != EXPECTED_N_MONTHS
) {
  stop(
    "IVT layer count differs from the expected ",
    EXPECTED_N_MONTHS,
    "-month circulation record."
  )
}

terra::time(
  uivt
) <- u_time

terra::time(
  vivt
) <- v_time

stage_done(t0)

stage_message(
  "Eastward IVT:",
  terra::nrow(uivt),
  "rows ×",
  terra::ncol(uivt),
  "columns ×",
  terra::nlyr(uivt),
  "months"
)

stage_message(
  "Northward IVT:",
  terra::nrow(vivt),
  "rows ×",
  terra::ncol(vivt),
  "columns ×",
  terra::nlyr(vivt),
  "months"
)


# ==============================================================================
# 13. Z500 CALENDAR-MONTH ANOMALY
# ==============================================================================

z_anom <- monthly_anomaly_raster(
  z500,
  z_dates,
  label = "Z500"
)

# Release source Z500 object after anomaly calculation.
rm(z500)

gc(
  verbose = FALSE
)


# ==============================================================================
# 14. EASTWARD IVT CALENDAR-MONTH ANOMALY
# ==============================================================================

u_anom <- monthly_anomaly_raster(
  uivt,
  z_dates,
  label = "eastward IVT"
)

gc(
  verbose = FALSE
)


# ==============================================================================
# 15. NORTHWARD IVT CALENDAR-MONTH ANOMALY
# ==============================================================================

v_anom <- monthly_anomaly_raster(
  vivt,
  z_dates,
  label = "northward IVT"
)

gc(
  verbose = FALSE
)


# ==============================================================================
# 16. IVT MAGNITUDE
#
#   IVT_mag(t,x) = sqrt(U(t,x)^2 + V(t,x)^2)
#
# The U/V component rasters are retained from Sections 14-15 and are therefore
# not reopened. terra raster algebra is used instead of terra::lapp(). This
# avoids the terra 1.9.x lapp() argument-dispatch failure encountered when
# passing function(u, v) to a two-component SpatRaster.
#
# The resulting IVT-magnitude raster remains disk-backed. Section 17 then
# applies the same calendar-month anomaly definition used for the other
# circulation fields.
# ==============================================================================

t0 <- stage_start(
  "Computing IVT magnitude"
)

if (!exists("uivt") || !inherits(uivt, "SpatRaster")) {
  stop(
    "Eastward IVT raster is unavailable for magnitude calculation."
  )
}

if (!exists("vivt") || !inherits(vivt, "SpatRaster")) {
  stop(
    "Northward IVT raster is unavailable for magnitude calculation."
  )
}

if (
  terra::nlyr(uivt) != EXPECTED_N_MONTHS ||
  terra::nlyr(vivt) != EXPECTED_N_MONTHS
) {
  stop(
    "IVT component layer count differs from the expected ",
    EXPECTED_N_MONTHS,
    "-month circulation record."
  )
}

u_mag_dates <- read_dates(
  uivt,
  "eastward IVT for magnitude"
)

v_mag_dates <- read_dates(
  vivt,
  "northward IVT for magnitude"
)

if (!identical(u_mag_dates, z_dates)) {
  stop(
    "Eastward IVT magnitude-input time axis differs from Z500."
  )
}

if (!identical(v_mag_dates, z_dates)) {
  stop(
    "Northward IVT magnitude-input time axis differs from Z500."
  )
}

stage_message(
  "Evaluating sqrt(U^2 + V^2) as a disk-backed raster expression"
)

ivt_mag <- sqrt(
  uivt^2 + vivt^2
)

terra::time(ivt_mag) <- z_dates

stage_done(t0)

# Keep the disk-backed U/V component references alive because the
# publication-map QA below uses the original component fields.  Retaining
# SpatRaster references does not materialize the full NetCDF fields in RAM.


# ==============================================================================
# 17. IVT-MAGNITUDE CALENDAR-MONTH ANOMALY
# ==============================================================================

ivt_mag_anom <- monthly_anomaly_raster(
  ivt_mag,
  z_dates,
  label = "IVT magnitude"
)

rm(ivt_mag)

gc(
  verbose = FALSE
)


# ==============================================================================
# 18. DOMAIN-MEAN MONTHLY SUMMARY
#
# These are means over the downloaded circulation domain.
# They are descriptive circulation-domain summaries, NOT Nechako-basin means.
# ==============================================================================

stage_message(
  "Computing domain-mean monthly circulation diagnostics"
)

mean_domain_with_progress <- function(x, label) {
  n <- terra::nlyr(x)
  out <- numeric(n)
  t0 <- Sys.time()
  progress_step.last <<- -Inf
  for (i in seq_len(n)) {
    out[i] <- terra::global(x[[i]], "mean", na.rm = TRUE)[1, 1]
    progress_step(i, n, paste0("Domain mean: ", label), t0, step = 0.10)
  }
  out
}

t0 <- Sys.time()
z500_domain <- mean_domain_with_progress(z_anom, "Z500")
ivt_mag_domain <- mean_domain_with_progress(ivt_mag_anom, "IVT magnitude")
ivt_east_domain <- mean_domain_with_progress(u_anom, "eastward IVT")
ivt_north_domain <- mean_domain_with_progress(v_anom, "northward IVT")

domain_summary <- data.frame(
  date = z_dates,
  year = as.integer(
    format(
      z_dates,
      "%Y"
    )
  ),
  month = as.integer(
    format(
      z_dates,
      "%m"
    )
  ),
  z500_anom_m = z500_domain,
  ivt_mag_anom = ivt_mag_domain,
  ivt_east_anom = ivt_east_domain,
  ivt_north_anom = ivt_north_domain,
  stringsAsFactors = FALSE
)

stage_done(t0)


DOMAIN_SUMMARY_FILE <- file.path(
  OUT_DIR,
  "circulation_domain_monthly_summary.csv"
)

write.csv(
  domain_summary,
  DOMAIN_SUMMARY_FILE,
  row.names = FALSE
)

stage_message(
  "Wrote:",
  DOMAIN_SUMMARY_FILE
)


# ==============================================================================
# 19. FIXED-WINDOW EVENT DEFINITIONS
#
# Principal drought episodes:
#   2022 May-October
#   2023 May-October
#   2024 May-October
#   2025 May-October
#
# El Nino context:
#   June 2023-May 2024
#
# These are descriptive fixed windows, not dynamically selected circulation
# regimes and not causal-attribution windows.
# ==============================================================================

episodes <- data.frame(
  event = c(
    "2022 drought",
    "2023 drought",
    "2024 drought",
    "2025 drought",
    "2023-2024 El Nino"
  ),
  start = as.Date(
    c(
      "2022-05-01",
      "2023-05-01",
      "2024-05-01",
      "2025-05-01",
      "2023-06-01"
    )
  ),
  end = as.Date(
    c(
      "2022-10-01",
      "2023-10-01",
      "2024-10-01",
      "2025-10-01",
      "2024-05-01"
    )
  ),
  stringsAsFactors = FALSE
)


# ==============================================================================
# 20. COMPOSITE FUNCTION
# ==============================================================================

composite_one <- function(
    start,
    end,
    label) {
  
  ii <- which(
    z_dates >= start &
      z_dates <= end
  )
  
  if (!length(ii)) {
    stop(
      "No circulation months in composite interval: ",
      label
    )
  }
  
  cat(
    sprintf(
      "  [%s] %-25s (%d months)\n",
      format(
        Sys.time(),
        "%H:%M:%S"
      ),
      label,
      length(ii)
    )
  )
  
  flush.console()
  
  list(
    z = terra::mean(
      z_anom[[ii]],
      na.rm = TRUE
    ),
    
    u = terra::mean(
      u_anom[[ii]],
      na.rm = TRUE
    ),
    
    v = terra::mean(
      v_anom[[ii]],
      na.rm = TRUE
    ),
    
    ivt = terra::mean(
      ivt_mag_anom[[ii]],
      na.rm = TRUE
    ),
    
    n = length(ii)
  )
}


# ==============================================================================
# 21. COMPUTE EVENT COMPOSITES
# ==============================================================================

t0 <- stage_start(
  "Computing event composites"
)

comps <- lapply(
  seq_len(
    nrow(episodes)
  ),
  function(i) {
    
    composite_one(
      start = episodes$start[i],
      end = episodes$end[i],
      label = episodes$event[i]
    )
  }
)

stage_done(t0)


# ==============================================================================
# 22. EVENT-COMPOSITE DOMAIN-MEAN SUMMARY
# ==============================================================================

t0 <- stage_start(
  "Computing event-composite domain means"
)

summary_rows <- dplyr::bind_rows(
  lapply(
    seq_along(comps),
    function(i) {
      
      c0 <- comps[[i]]
      
      data.frame(
        event = episodes$event[i],
        start = episodes$start[i],
        end = episodes$end[i],
        n_months = c0$n,
        
        domain_mean_z500_anom_m =
          terra::global(
            c0$z,
            "mean",
            na.rm = TRUE
          )[1, 1],
        
        domain_mean_ivt_mag_anom =
          terra::global(
            c0$ivt,
            "mean",
            na.rm = TRUE
          )[1, 1],
        
        domain_mean_ivt_east_anom =
          terra::global(
            c0$u,
            "mean",
            na.rm = TRUE
          )[1, 1],
        
        domain_mean_ivt_north_anom =
          terra::global(
            c0$v,
            "mean",
            na.rm = TRUE
          )[1, 1],
        
        stringsAsFactors = FALSE
      )
    }
  )
)

stage_done(t0)


EVENT_SUMMARY_FILE <- file.path(
  OUT_DIR,
  "circulation_event_composites.csv"
)

write.csv(
  summary_rows,
  EVENT_SUMMARY_FILE,
  row.names = FALSE
)

stage_message(
  "Wrote:",
  EVENT_SUMMARY_FILE
)


# ==============================================================================
# 23. PERSIST FULL EVENT-COMPOSITE FIELDS AS NETCDF
# ==============================================================================

t0 <- stage_start(
  "Writing event-composite NetCDF"
)

stack_list <- list()
stack_names <- character()

stack_t0 <- Sys.time()
progress_step.last <- -Inf
for (i in seq_along(comps)) {
  
  tag <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    episodes$event[i]
  )
  
  stack_list <- c(
    stack_list,
    list(
      comps[[i]]$z,
      comps[[i]]$ivt,
      comps[[i]]$u,
      comps[[i]]$v
    )
  )
  
  stack_names <- c(
    stack_names,
    paste0(
      tag,
      "_Z500_anom_m"
    ),
    paste0(
      tag,
      "_IVTmag_anom"
    ),
    paste0(
      tag,
      "_IVTeast_anom"
    ),
    paste0(
      tag,
      "_IVTnorth_anom"
    )
  )
  progress_step(i, length(comps), "Assembling composite layers", stack_t0, step = 0.20)
}

comp_stack <- do.call(
  c,
  stack_list
)

names(comp_stack) <- stack_names

COMPOSITE_NETCDF <- file.path(
  OUT_DIR,
  "circulation_event_composite_maps.nc"
)

terra::writeCDF(
  comp_stack,
  COMPOSITE_NETCDF,
  overwrite = TRUE,
  compression = 9
)

stage_done(t0)

stage_message(
  "Wrote:",
  COMPOSITE_NETCDF
)


# ==============================================================================
# 24. PUBLICATION DIAGNOSTIC DATA
# ==============================================================================

t0 <- stage_start(
  "Preparing circulation composite figure"
)

plot_df <- summary_rows |>
  tidyr::pivot_longer(
    cols = c(
      domain_mean_z500_anom_m,
      domain_mean_ivt_mag_anom
    ),
    names_to = "diagnostic",
    values_to = "anomaly"
  ) |>
  dplyr::mutate(
    diagnostic = factor(
      diagnostic,
      levels = c(
        "domain_mean_z500_anom_m",
        "domain_mean_ivt_mag_anom"
      ),
      labels = c(
        "Z500 anomaly (m)",
        "IVT-magnitude anomaly"
      )
    )
  )

p <- ggplot2::ggplot(
  plot_df,
  ggplot2::aes(
    x = event,
    y = anomaly,
    group = diagnostic
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linewidth = 0.35
  ) +
  ggplot2::geom_col() +
  ggplot2::facet_wrap(
    ~ diagnostic,
    scales = "free_y",
    ncol = 1
  ) +
  ggplot2::labs(
    x = NULL,
    y = NULL,
    title =
      "Large-scale circulation context of recent Nechako drought episodes",
    subtitle =
      "Calendar-month anomalies relative to the 1950-2025 climatology"
  ) +
  ggplot2::theme_bw(
    base_size = 10
  ) +
  ggplot2::theme(
    axis.text.x =
      ggplot2::element_text(
        angle = 30,
        hjust = 1
      )
  )

stage_done(t0)


# ==============================================================================
# 25. SAVE FIGURES
# ==============================================================================

FIG_PNG <- file.path(
  OUT_DIR,
  "Fig_circulation_event_composites.png"
)

FIG_PDF <- file.path(
  OUT_DIR,
  "Fig_circulation_event_composites.pdf"
)


t0 <- stage_start(
  "Writing circulation composite PNG"
)

ggplot2::ggsave(
  FIG_PNG,
  p,
  width = 9,
  height = 7,
  dpi = FIG_DPI
)

stage_done(t0)


t0 <- stage_start(
  "Writing circulation composite PDF"
)

ggplot2::ggsave(
  FIG_PDF,
  p,
  width = 9,
  height = 7
)

stage_done(t0)


# ==============================================================================
# 26. PUBLICATION-QUALITY SPATIAL CIRCULATION MAPS
# ==============================================================================
# The manuscript circulation figure uses two spatial scales:
#   A. Regional synoptic context: 170°W–100°W, 25°N–65°N.
#   B. Nechako-focused detail: basin bounding box plus a small geographic buffer.
#
# The IVT anomaly magnitude is the raster field; Z500 anomaly is shown with
# contours; anomalous IVT direction is shown with sparse, explicit arrowheads.
# Vector placement is performed on a regular 2-D grid.  The previous approach
# selected every nth flattened raster row, which produced visually misleading
# diagonal trails.
# ==============================================================================

CIRCULATION_MAP_PNG <- file.path(
  OUT_DIR,
  "Fig_circulation_event_composite_maps.png"
)

CIRCULATION_MAP_PDF <- file.path(
  OUT_DIR,
  "Fig_circulation_event_composite_maps.pdf"
)

BASIN_MAP_FILE <- file.path(
  WD_PATH,
  "Spatial",
  "nechakoBound_dissolve.kmz"
)

if (!file.exists(BASIN_MAP_FILE)) {
  BASIN_MAP_FILE <- file.path(
    WD_PATH,
    "Spatial",
    "nechakoBound_dissolve.shp"
  )
}

if (!file.exists(BASIN_MAP_FILE)) {
  stop(
    "Cannot produce circulation maps: Nechako basin boundary not found. " ,
    "Expected KMZ or SHP in Spatial/."
  )
}

# ------------------------------------------------------------------------------
# Basin boundary: use the shared KMZ-aware reader.
# ------------------------------------------------------------------------------
basin_map_sf <- sf::st_as_sf(
  load_basin_vect(BASIN_MAP_FILE)
)

basin_map_sf <- sf::st_make_valid(
  basin_map_sf
)

basin_map_sf <- sf::st_transform(
  basin_map_sf,
  4326
)

# Automatically derive a readable basin-detail extent.
bb <- sf::st_bbox(
  basin_map_sf
)

basin_xrange <- as.numeric(bb[c("xmin", "xmax")])
basin_yrange <- as.numeric(bb[c("ymin", "ymax")])

basin_dx <- diff(basin_xrange)
basin_dy <- diff(basin_yrange)

basin_buffer_x <- max(1.5, 0.18 * basin_dx)
basin_buffer_y <- max(1.5, 0.18 * basin_dy)

NECHAKO_XLIM <- c(
  max(-170, basin_xrange[1] - basin_buffer_x),
  min(-100, basin_xrange[2] + basin_buffer_x)
)

NECHAKO_YLIM <- c(
  max(25, basin_yrange[1] - basin_buffer_y),
  min(65, basin_yrange[2] + basin_buffer_y)
)

REGIONAL_XLIM <- c(
  -170,
  -100
)

REGIONAL_YLIM <- c(
  25,
  65
)

# ------------------------------------------------------------------------------
# IVT component QA.  The original U/V objects are deliberately retained until
# this test has completed.
# ------------------------------------------------------------------------------
u_mean_abs <- as.numeric(
  terra::global(
    abs(uivt - vivt),
    "mean",
    na.rm = TRUE
  )[1, 1]
)

if (!is.finite(u_mean_abs)) {
  u_mean_abs <- NA_real_
}

u_sd <- as.numeric(
  terra::global(
    uivt,
    "sd",
    na.rm = TRUE
  )[1, 1]
)

v_sd <- as.numeric(
  terra::global(
    vivt,
    "sd",
    na.rm = TRUE
  )[1, 1]
)

DIRECTIONAL_IVT_VALID <-
  is.finite(u_mean_abs) &&
  u_mean_abs > max(.Machine$double.eps * 100, 1e-10) &&
  is.finite(u_sd) &&
  is.finite(v_sd) &&
  (u_sd > 0 || v_sd > 0)

if (DIRECTIONAL_IVT_VALID) {
  stage_message(
    sprintf(
      "IVT component QA passed: mean |U-V| = %.6g; directional arrows enabled.",
      u_mean_abs
    )
  )
} else {
  stage_message(
    sprintf(
      "IVT component QA failed: mean |U-V| = %s; directional arrows disabled.",
      ifelse(
        is.finite(u_mean_abs),
        format(u_mean_abs, scientific = TRUE),
        "NA"
      )
    )
  )
}

# ------------------------------------------------------------------------------
# Common colour scale and Z500 contour interval.
# Using common limits makes the five event composites directly comparable.
# ------------------------------------------------------------------------------
all_ivt_values <- unlist(
  lapply(
    comps,
    function(x) {
      terra::values(x$ivt, mat = FALSE)
    }
  ),
  use.names = FALSE
)

all_ivt_values <- all_ivt_values[
  is.finite(all_ivt_values)
]

IVT_LIMIT <- if (length(all_ivt_values)) {
  max(abs(all_ivt_values), na.rm = TRUE)
} else {
  1
}

if (!is.finite(IVT_LIMIT) || IVT_LIMIT <= 0) {
  IVT_LIMIT <- 1
}

all_z_values <- unlist(
  lapply(
    comps,
    function(x) {
      terra::values(x$z, mat = FALSE)
    }
  ),
  use.names = FALSE
)

all_z_values <- all_z_values[
  is.finite(all_z_values)
]

Z_LIMIT <- if (length(all_z_values)) {
  max(abs(all_z_values), na.rm = TRUE)
} else {
  1
}

if (!is.finite(Z_LIMIT) || Z_LIMIT <= 0) {
  Z_LIMIT <- 1
}

Z_BREAKS <- pretty(
  c(-Z_LIMIT, Z_LIMIT),
  n = 9
)

Z_BREAKS <- Z_BREAKS[
  Z_BREAKS >= -Z_LIMIT &
    Z_BREAKS <= Z_LIMIT
]

Z_BREAKS <- sort(
  unique(
    c(
      Z_BREAKS,
      0
    )
  )
)

# ------------------------------------------------------------------------------
# Helper: select a regular 2-D vector grid.
# ------------------------------------------------------------------------------
thin_vector_grid <- function(
    df,
    nx = 10L,
    ny = 8L) {

  if (!nrow(df)) {
    return(df)
  }

  ux <- sort(unique(df$x))
  uy <- sort(unique(df$y))

  if (!length(ux) || !length(uy)) {
    return(df[0, , drop = FALSE])
  }

  ix <- unique(
    pmax(
      1L,
      pmin(
        length(ux),
        round(
          seq(
            1,
            length(ux),
            length.out = min(nx, length(ux))
          )
        )
      )
    )
  )

  iy <- unique(
    pmax(
      1L,
      pmin(
        length(uy),
        round(
          seq(
            1,
            length(uy),
            length.out = min(ny, length(uy))
          )
        )
      )
    )
  )

  df[
    df$x %in% ux[ix] &
      df$y %in% uy[iy],
    ,
    drop = FALSE
  ]
}

# ------------------------------------------------------------------------------
# Helper: construct visually readable arrows while preserving U/V direction.
# Arrow length is a display parameter only; it does not encode magnitude.
# Longitude increments are adjusted by cos(latitude) so an arrow represents
# the same physical display length in the east-west and north-south directions.
# ------------------------------------------------------------------------------
make_arrow_data <- function(
    d,
    xlim,
    ylim,
    nx,
    ny) {

  ar <- d[
    is.finite(d$u) &
      is.finite(d$v) &
      is.finite(d$x) &
      is.finite(d$y),
    ,
    drop = FALSE
  ]

  ar <- ar[
    ar$x >= xlim[1] &
      ar$x <= xlim[2] &
      ar$y >= ylim[1] &
      ar$y <= ylim[2],
    ,
    drop = FALSE
  ]

  ar <- thin_vector_grid(
    ar,
    nx = nx,
    ny = ny
  )

  if (!nrow(ar)) {
    return(ar)
  }

  theta <- atan2(
    ar$v,
    ar$u
  )

  # Display length in latitude degrees.
  map_height <- diff(ylim)
  arrow_length <- max(0.18, min(1.20, map_height / 38))

  lat_cos <- pmax(
    cos(ar$y * pi / 180),
    0.25
  )

  ar$u_plot <-
    arrow_length *
    cos(theta) /
    lat_cos

  ar$v_plot <-
    arrow_length *
    sin(theta)

  # Keep the complete arrow inside the plotting window where possible.
  ar$xend <- pmax(
    xlim[1],
    pmin(
      xlim[2],
      ar$x + ar$u_plot
    )
  )

  ar$yend <- pmax(
    ylim[1],
    pmin(
      ylim[2],
      ar$y + ar$v_plot
    )
  )

  ar
}

# ------------------------------------------------------------------------------
# Map constructor.
# ------------------------------------------------------------------------------
make_circulation_map <- function(
    i,
    scale = c("regional", "nechako")) {

  scale <- match.arg(scale)

  zdf <- terra::as.data.frame(
    comps[[i]]$z,
    xy = TRUE,
    na.rm = FALSE
  )

  udf <- terra::as.data.frame(
    comps[[i]]$u,
    xy = TRUE,
    na.rm = FALSE
  )

  vdf <- terra::as.data.frame(
    comps[[i]]$v,
    xy = TRUE,
    na.rm = FALSE
  )

  ivtdf <- terra::as.data.frame(
    comps[[i]]$ivt,
    xy = TRUE,
    na.rm = FALSE
  )

  names(zdf)[3] <- "z500"
  names(udf)[3] <- "u"
  names(vdf)[3] <- "v"
  names(ivtdf)[3] <- "ivt_mag"

  d <- zdf |>
    dplyr::left_join(
      udf,
      by = c("x", "y")
    ) |>
    dplyr::left_join(
      vdf,
      by = c("x", "y")
    ) |>
    dplyr::left_join(
      ivtdf,
      by = c("x", "y")
    )

  if (scale == "regional") {
    xlim <- REGIONAL_XLIM
    ylim <- REGIONAL_YLIM
    nx <- 10L
    ny <- 7L
    scale_label <- "Regional circulation context"
  } else {
    xlim <- NECHAKO_XLIM
    ylim <- NECHAKO_YLIM
    nx <- 8L
    ny <- 6L
    scale_label <- "Nechako Basin circulation detail"
  }

  # Retain only the displayed map extent.
  d <- d[
    d$x >= xlim[1] &
      d$x <= xlim[2] &
      d$y >= ylim[1] &
      d$y <= ylim[2],
    ,
    drop = FALSE
  ]

  p0 <- ggplot2::ggplot(
    d,
    ggplot2::aes(
      x = x,
      y = y
    )
  ) +
    ggplot2::geom_raster(
      ggplot2::aes(
        fill = ivt_mag
      )
    ) +
    ggplot2::geom_contour(
      ggplot2::aes(
        z = z500
      ),
      breaks = Z_BREAKS,
      linewidth = 0.22,
      colour = "grey25"
    ) +
    ggplot2::geom_sf(
      data = basin_map_sf,
      inherit.aes = FALSE,
      fill = NA,
      linewidth = if (scale == "nechako") 1.05 else 0.75,
      colour = "black"
    ) +
    ggplot2::coord_sf(
      xlim = xlim,
      ylim = ylim,
      expand = FALSE
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#B2182B",
      mid = "white",
      high = "#2166AC",
      midpoint = 0,
      limits = c(-IVT_LIMIT, IVT_LIMIT),
      oob = scales::squish,
      na.value = "white"
    ) +
    ggplot2::labs(
      title = episodes$event[i],
      subtitle = scale_label,
      fill = "IVT-magnitude
anomaly
(kg m⁻¹ s⁻¹)",
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_bw(
      base_size = if (scale == "nechako") 10 else 9
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = if (scale == "nechako") 11 else 10.5
      ),
      plot.subtitle = ggplot2::element_text(
        size = if (scale == "nechako") 9 else 8.5
      ),
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(
        size = if (scale == "nechako") 8.5 else 8
      ),
      legend.title = ggplot2::element_text(
        size = 8.5
      ),
      legend.text = ggplot2::element_text(
        size = 8
      )
    )

  if (DIRECTIONAL_IVT_VALID) {

    ar <- make_arrow_data(
      d,
      xlim = xlim,
      ylim = ylim,
      nx = nx,
      ny = ny
    )

    if (nrow(ar)) {
      p0 <- p0 +
        ggplot2::geom_segment(
          data = ar,
          ggplot2::aes(
            x = x,
            y = y,
            xend = xend,
            yend = yend
          ),
          inherit.aes = FALSE,
          linewidth = if (scale == "nechako") 0.42 else 0.34,
          colour = "black",
          lineend = "round",
          arrow = grid::arrow(
            length = grid::unit(
              if (scale == "nechako") 0.18 else 0.16,
              "cm"
            ),
            type = "closed"
          )
        )
    }
  }

  p0
}

# ------------------------------------------------------------------------------
# Render the two-scale manuscript figure.
# ------------------------------------------------------------------------------
t0 <- stage_start(
  sprintf(
    "Rendering %d event composites at regional + Nechako scales",
    nrow(episodes)
  )
)

regional_plots <- lapply(
  seq_len(nrow(episodes)),
  function(i) {
    stage_message(
      sprintf(
        "  Regional map %d/%d: %s",
        i,
        nrow(episodes),
        episodes$event[i]
      )
    )
    make_circulation_map(
      i,
      scale = "regional"
    )
  }
)

nechako_plots <- lapply(
  seq_len(nrow(episodes)),
  function(i) {
    stage_message(
      sprintf(
        "  Nechako-detail map %d/%d: %s",
        i,
        nrow(episodes),
        episodes$event[i]
      )
    )
    make_circulation_map(
      i,
      scale = "nechako"
    )
  }
)

# Arrange each event as one row: regional context | Nechako detail.
map_plots <- vector(
  "list",
  length = 2L * nrow(episodes)
)

for (i in seq_len(nrow(episodes))) {
  map_plots[[2L * i - 1L]] <- regional_plots[[i]]
  map_plots[[2L * i]] <- nechako_plots[[i]]
}

p_map <- patchwork::wrap_plots(
  map_plots,
  ncol = 2,
  guides = "collect"
) +
  patchwork::plot_annotation(
    title = "Circulation context of recent Nechako drought episodes",
    subtitle = paste0(
      "Left: regional IVT-magnitude anomalies and Z500 anomalies. ",
      "Right: Nechako Basin detail. ",
      if (DIRECTIONAL_IVT_VALID) {
        "Black arrows show anomalous IVT direction; arrow length is a display parameter only."
      } else {
        "Directional IVT arrows omitted after component QA."
      }
    ),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 15
      ),
      plot.subtitle = ggplot2::element_text(
        size = 9.5
      ),
      legend.position = "right"
    )
  )

stage_message(
  "  Saving two-scale spatial circulation PNG..."
)

ggplot2::ggsave(
  CIRCULATION_MAP_PNG,
  p_map,
  width = 14,
  height = 15,
  dpi = FIG_DPI
)

stage_message(
  "  Saving two-scale spatial circulation PDF..."
)

ggplot2::ggsave(
  CIRCULATION_MAP_PDF,
  p_map,
  width = 14,
  height = 15
)

stage_done(t0)

stage_message(
  "Wrote:",
  CIRCULATION_MAP_PNG
)

stage_message(
  "Wrote:",
  CIRCULATION_MAP_PDF
)

# Separate versions are useful for manuscript main-text/supplement decisions.
regional_only <- patchwork::wrap_plots(
  regional_plots,
  ncol = 2,
  guides = "collect"
) +
  patchwork::plot_annotation(
    title = "Regional circulation context of recent Nechako drought episodes",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 15
      ),
      legend.position = "right"
    )
  )

nechako_only <- patchwork::wrap_plots(
  nechako_plots,
  ncol = 2,
  guides = "collect"
) +
  patchwork::plot_annotation(
    title = "Nechako Basin circulation detail for recent drought episodes",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 15
      ),
      legend.position = "right"
    )
  )

ggplot2::ggsave(
  file.path(
    OUT_DIR,
    "Fig_circulation_event_composite_maps_regional.png"
  ),
  regional_only,
  width = 13,
  height = 10,
  dpi = FIG_DPI
)

ggplot2::ggsave(
  file.path(
    OUT_DIR,
    "Fig_circulation_event_composite_maps_nechako.png"
  ),
  nechako_only,
  width = 13,
  height = 10,
  dpi = FIG_DPI
)

ggplot2::ggsave(
  file.path(
    OUT_DIR,
    "Fig_circulation_event_composite_maps_regional.pdf"
  ),
  regional_only,
  width = 13,
  height = 10
)

ggplot2::ggsave(
  file.path(
    OUT_DIR,
    "Fig_circulation_event_composite_maps_nechako.pdf"
  ),
  nechako_only,
  width = 13,
  height = 10
)

# ==============================================================================
# 27. CLEANUP
# ==============================================================================

rm(
  z_anom,
  u_anom,
  v_anom,
  ivt_mag_anom,
  comps,
  comp_stack
)

gc(
  verbose = FALSE
)


# ==============================================================================
# 29. FINAL REPORT
# ==============================================================================

cat(
  "\n============================================================\n",
  "14circulation_context_v2.R completed successfully.\n",
  "============================================================\n",
  "Analysis period: ",
  format(
    ANALYSIS_START,
    "%Y-%m"
  ),
  " through ",
  format(
    ANALYSIS_END,
    "%Y-%m"
  ),
  " (",
  EXPECTED_N_MONTHS,
  " months)\n",
  sep = ""
)

cat(
  "Z500 input: ",
  Z_FILE,
  "\n",
  sep = ""
)

cat(
  "IVT input:  ",
  IVT_FILE,
  "\n",
  sep = ""
)

cat(
  "Output directory: ",
  OUT_DIR,
  "\n\n",
  sep = ""
)

cat(
  "Outputs:\n",
  "  - circulation_domain_monthly_summary.csv\n",
  "  - circulation_event_composites.csv\n",
  "  - circulation_event_composite_maps.nc\n",
  "  - Fig_circulation_event_composites.png\n",
  "  - Fig_circulation_event_composites.pdf\n",
  "  - Fig_circulation_event_composite_maps.png\n",
  "  - Fig_circulation_event_composite_maps.pdf\n",
  sep = ""
)

cat(
  "\nScientific interpretation:\n",
  "  Circulation diagnostics provide large-scale context only.\n",
  "  They are not causal attribution of drought severity.\n",
  sep = ""
)

cat(
  "\nMemory architecture:\n",
  "  - Each required circulation field is materialized once for anomaly/composite calculations.\n",
  "  - Calendar-month climatologies use terra::tapp().\n",
  "  - No repeated 76-layer disk-backed monthly reductions are used.\n",
  "  - Intermediate fields are released before the next major calculation.\n",
  sep = ""
)

cat(
  "\n============================================================\n"
)

invisible(TRUE)

# ==============================================================================
# 30. BASIN-SCALE IVT/HYDROLOGY ENHANCEMENT
# ==============================================================================
# Adds direction, boundary moisture flux, MFC, precipitation association, and vector maps.
source("14b_circulation_hydrology_v2.R")
