################################################################################
# 14b_circulation_hydrology_v2.R
# BASIN-SCALE IVT DIRECTION, MOISTURE-FLUX, CONVERGENCE, AND PRECIPITATION LINK
#
# Called automatically by 14circulation_context_v2.R after the existing
# circulation-context products have been generated.
#
# Added diagnostics:
#   1. IVT direction from existing U/V components
#   2. Nechako boundary-normal moisture flux by west/east/south/north sectors
#   3. Net atmospheric moisture import into the basin
#   4. Basin-mean moisture-flux convergence (MFC)
#   5. Basin-integrated MFC
#   6. Basin precipitation anomalies
#   7. IVT/transport/MFC versus precipitation association
#   8. Event summaries and IVT vector maps
#
# Explicitly NOT included:
#   IVT directional persistence
#   850/925-hPa moisture flux
#   atmospheric-river detection
#   Lagrangian moisture-source tracking
#   full atmospheric moisture budget
################################################################################

setwd(Sys.getenv("NECHAKO_WD", "D:/Nechako_Drought/Nechako/"))
source("DROUGHT_ANALYSIS_utils_v2.R")
utils_load_packages(c("terra", "sf", "dplyr", "tidyr", "ggplot2", "patchwork"))

if (!dir.exists(WD_PATH)) stop("Working directory not found: ", WD_PATH)
setwd(WD_PATH)

OUT_DIR <- file.path(WD_PATH, "circulation_context")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

Z_FILE <- file.path(WD_PATH, "monthly_data_direct", "geopotential_500hPa_monthly.nc")
IVT_FILE <- file.path(WD_PATH, "monthly_data_direct", "ivt_components_monthly.nc")
P_FILE <- file.path(WD_PATH, "monthly_data_direct", "total_precipitation_monthly.nc")
BASIN_FILE <- file.path(WD_PATH, "Spatial", "nechakoBound_dissolve.kmz")

for (f in c(Z_FILE, IVT_FILE, P_FILE, BASIN_FILE)) {
  if (!file.exists(f)) stop("Required input not found: ", f)
}

START <- as.Date("1950-01-01")
END <- as.Date("2025-12-01")
DATES <- seq(START, END, by = "month")
N <- length(DATES)
if (N != 912L) stop("Expected 912 monthly observations; found ", N)

cat("\n============================================================\n")
cat("14b circulation hydrology diagnostics — progress-enabled run\n")
cat("Analysis period: 1950-01 through 2025-12 (912 months)\n")
cat("The script will print a timestamp before every expensive terra operation.\n")
cat("A long-running operation with no new line is still executing, not waiting for input.\n")
cat("============================================================\n\n")
flush.console()

read_dates <- function(x, label) {
  tt <- terra::time(x)
  if (is.null(tt) || length(tt) != terra::nlyr(x) || anyNA(tt))
    stop("Invalid time axis: ", label)
  d <- as.Date(tt)
  d <- as.Date(format(d, "%Y-%m-01"))
  if (!identical(d, DATES)) stop(label, " does not match 1950-01 to 2025-12.")
  d
}

# -----------------------------------------------------------------------------
# CONSOLE PROGRESS / TIMING HELPERS
# -----------------------------------------------------------------------------
# These diagnostics are deliberately explicit because several terra operations
# below can run for minutes without returning control to the R console.
# A stage message is printed BEFORE each expensive operation, and calendar-month
# climatologies are calculated month-by-month so progress is visible rather than
# appearing as a stalled terra progress bar.
format_elapsed_local <- function(seconds) {
  if (!is.finite(seconds)) return("--")
  if (seconds < 60) return(sprintf("%.0fs", seconds))
  sprintf("%dm %.0fs", floor(seconds / 60), seconds %% 60)
}

stage_message <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
  flush.console()
  invisible(NULL)
}

stage_start <- function(label) {
  stage_message("START: ", label)
  Sys.time()
}

stage_done <- function(t0, label = "DONE") {
  elapsed <- as.numeric(Sys.time() - t0, units = "secs")
  stage_message(label, " | elapsed = ", format_elapsed_local(elapsed))
  invisible(elapsed)
}

progress_line <- function(i, n, label, t0) {
  elapsed <- as.numeric(Sys.time() - t0, units = "secs")
  pct <- 100 * i / n
  rate <- if (elapsed > 0) i / elapsed else NA_real_
  eta <- if (is.finite(rate) && rate > 0) (n - i) / rate else NA_real_
  stage_message(
    sprintf(
      "%s | %5.1f%% (%d/%d) | elapsed %s | ETA %s",
      label, pct, i, n, format_elapsed_local(elapsed), format_elapsed_local(eta)
    )
  )
}

# Calendar-month climatology equivalent to terra::tapp(), but intentionally
# evaluated as 12 visible stages. This changes no scientific definition: each
# grid cell is still standardized against the 1950-2025 mean for its calendar
# month. It makes a long NetCDF-backed reduction observable in the console.
monthly_anomaly <- function(x, dates = DATES, label = "raster") {
  m <- as.integer(format(dates, "%m"))
  if (length(m) != terra::nlyr(x)) {
    stop("Monthly anomaly time axis mismatch for ", label, ".")
  }
  
  t0 <- stage_start(
    sprintf(
      "Calendar-month climatology + anomaly: %s (%d layers; %d x %d grid)",
      label, terra::nlyr(x), terra::nrow(x), terra::ncol(x)
    )
  )
  
  clim_list <- vector("list", 12L)
  for (mm in 1:12) {
    ii <- which(m == mm)
    stage_message(
      sprintf(
        "  Month %02d/12 (%s): averaging %d historical layers...",
        mm, month.abb[mm], length(ii)
      )
    )
    tmm <- Sys.time()
    clim_list[[mm]] <- terra::mean(x[[ii]], na.rm = TRUE)
    stage_message(
      sprintf(
        "  Month %02d/12 (%s): complete | %s",
        mm, month.abb[mm],
        format_elapsed_local(as.numeric(Sys.time() - tmm, units = "secs"))
      )
    )
  }
  
  clim <- do.call(c, clim_list)
  names(clim) <- sprintf("clim_%02d", 1:12)
  
  stage_message("  Applying calendar-month climatologies to all layers...")
  out <- x - clim[[m]]
  terra::time(out) <- terra::time(x)
  
  rm(clim, clim_list)
  gc(verbose = FALSE)
  stage_done(t0, sprintf("COMPLETE: %s anomaly", label))
  out
}

extract_ivt <- function(f) {
  s <- try(terra::sds(f), silent = TRUE)
  if (!inherits(s, "try-error") && length(s) >= 2L) {
    nm <- tolower(names(s))
    ie <- grep("eastward|east", nm)
    in_ <- grep("northward|north", nm)
    if (!length(ie)) ie <- 1L
    if (!length(in_)) in_ <- 2L
    return(list(east = s[[ie[1]]], north = s[[in_[1]]]))
  }
  x <- terra::rast(f)
  nm <- tolower(names(x))
  ie <- grep("eastward|east", nm)
  in_ <- grep("northward|north", nm)
  if (!length(ie) || !length(in_)) stop("Could not identify IVT U/V components in ", f)
  list(east = x[[ie[1]]], north = x[[in_[1]]])
}

mean_raster_ts <- function(x, label) {
  message("Computing: ", label)
  v <- terra::global(x, "mean", na.rm = TRUE)[, 1]
  if (length(v) != N) stop(label, " returned an invalid time series length.")
  v
}

# -----------------------------------------------------------------------------
# 1. LOAD IVT
# -----------------------------------------------------------------------------
t0 <- stage_start("Opening IVT NetCDF and identifying eastward/northward components")
ivt <- extract_ivt(IVT_FILE)
u <- ivt$east
v <- ivt$north
stage_message(sprintf("  IVT grid: %d rows x %d cols x %d months", terra::nrow(u), terra::ncol(u), terra::nlyr(u)))
ud <- read_dates(u, "eastward IVT")
vd <- read_dates(v, "northward IVT")
if (!identical(ud, vd)) stop("Eastward/northward IVT time axes differ.")
stage_done(t0, "COMPLETE: IVT input opened and validated")

u_anom <- monthly_anomaly(u, label = "eastward IVT")
v_anom <- monthly_anomaly(v, label = "northward IVT")

t0 <- stage_start("Computing IVT magnitude sqrt(U^2 + V^2) across 912 layers")
ivt_mag <- sqrt(u^2 + v^2)
terra::time(ivt_mag) <- terra::time(u)
stage_done(t0, "COMPLETE: IVT magnitude")

ivt_mag_anom <- monthly_anomaly(ivt_mag, label = "IVT magnitude")

t0 <- stage_start("Computing 912-month domain means for IVT anomalies")
ivt_u_ts <- mean_raster_ts(u_anom, "domain eastward IVT anomaly")
ivt_v_ts <- mean_raster_ts(v_anom, "domain northward IVT anomaly")
ivt_mag_ts <- mean_raster_ts(ivt_mag_anom, "domain IVT magnitude anomaly")
stage_done(t0, "COMPLETE: IVT domain means")
ivt_dir_ts <- (atan2(ivt_v_ts, ivt_u_ts) * 180 / pi) %% 360
ivt_dir_ts[!is.finite(ivt_dir_ts) | (ivt_u_ts == 0 & ivt_v_ts == 0)] <- NA_real_

# -----------------------------------------------------------------------------
# 2. BASIN GEOMETRY AND BOUNDARY SEGMENTS
# -----------------------------------------------------------------------------
basin <- load_basin_vect(BASIN_FILE)
basin_eq <- terra::project(basin, EQUAL_AREA_CRS)
basin_area_m2 <- sum(terra::expanse(basin_eq, unit = "m"))
if (!is.finite(basin_area_m2) || basin_area_m2 <= 0) stop("Invalid basin area.")

bsf <- sf::st_as_sf(basin_eq)
boundary <- sf::st_boundary(bsf)
boundary <- sf::st_cast(boundary, "LINESTRING")
boundary <- sf::st_segmentize(boundary, dfMaxLength = 10000)
centroid <- sf::st_coordinates(sf::st_centroid(sf::st_union(bsf)))[1, ]

segments <- list()
k <- 0L
for (li in seq_len(nrow(boundary))) {
  xy <- sf::st_coordinates(boundary[li, ])
  if (nrow(xy) < 2L) next
  for (j in seq_len(nrow(xy) - 1L)) {
    p1 <- xy[j, 1:2]
    p2 <- xy[j + 1L, 1:2]
    dx <- p2[1] - p1[1]
    dy <- p2[2] - p1[2]
    L <- sqrt(dx^2 + dy^2)
    if (!is.finite(L) || L <= 0) next
    mid <- c(mean(c(p1[1], p2[1])), mean(c(p1[2], p2[2])))
    n1 <- c(-dy, dx) / L
    n2 <- c(dy, -dx) / L
    radial <- mid - centroid
    nout <- if (sum(n1 * radial) >= sum(n2 * radial)) n1 else n2
    k <- k + 1L
    segments[[k]] <- data.frame(
      segment_id = k, line_id = li,
      x_mid = mid[1], y_mid = mid[2],
      length_m = L, nx_out = nout[1], ny_out = nout[2]
    )
  }
}
segments <- dplyr::bind_rows(segments)
if (!nrow(segments)) stop("No basin boundary segments were generated.")
segments$sector <- ifelse(
  abs(segments$nx_out) >= abs(segments$ny_out),
  ifelse(segments$nx_out < 0, "west", "east"),
  ifelse(segments$ny_out < 0, "south", "north")
)

# Local geographic east/north unit vectors expressed in BC-Albers coordinates.
mid_eq <- sf::st_as_sf(segments, coords = c("x_mid", "y_mid"), crs = EQUAL_AREA_CRS, remove = FALSE)
mid_native <- sf::st_transform(mid_eq, sf::st_crs(u))
mid_native_xy <- sf::st_coordinates(mid_native)
mid_eq_xy <- sf::st_coordinates(mid_eq)
EPS_DEG <- 0.01

make_perturbed <- function(xy, dx = 0, dy = 0) {
  xy2 <- xy
  xy2[, 1] <- xy2[, 1] + dx
  xy2[, 2] <- xy2[, 2] + dy
  sf::st_as_sf(
    data.frame(lon = xy2[, 1], lat = xy2[, 2]),
    coords = c("lon", "lat"),
    crs = sf::st_crs(u)
  )
}

east_eq <- sf::st_coordinates(sf::st_transform(make_perturbed(mid_native_xy, EPS_DEG, 0), sf::st_crs(EQUAL_AREA_CRS)))
north_eq <- sf::st_coordinates(sf::st_transform(make_perturbed(mid_native_xy, 0, EPS_DEG), sf::st_crs(EQUAL_AREA_CRS)))

eb <- east_eq - mid_eq_xy
nb <- north_eq - mid_eq_xy
eb <- eb / sqrt(rowSums(eb^2))
nb <- nb / sqrt(rowSums(nb^2))
segments$east_x <- eb[, 1]
segments$east_y <- eb[, 2]
segments$north_x <- nb[, 1]
segments$north_y <- nb[, 2]

# -----------------------------------------------------------------------------
# 3. BOUNDARY-NORMAL IVT FLUX
# -----------------------------------------------------------------------------
pts <- terra::vect(
  data.frame(x = mid_native_xy[, 1], y = mid_native_xy[, 2]),
  geom = c("x", "y"), crs = terra::crs(u)
)
t0 <- stage_start(sprintf("Extracting U/V at %d basin-boundary segment midpoints across 912 months", nrow(segments)))
U <- as.matrix(terra::extract(u, pts, ID = FALSE))
V <- as.matrix(terra::extract(v, pts, ID = FALSE))
stage_done(t0, "COMPLETE: boundary U/V extraction")
if (nrow(U) != nrow(segments) || nrow(V) != nrow(segments)) stop("Boundary extraction size mismatch.")

qx <- U * matrix(segments$east_x, nrow = nrow(segments), ncol = N) +
  V * matrix(segments$north_x, nrow = nrow(segments), ncol = N)
qy <- U * matrix(segments$east_y, nrow = nrow(segments), ncol = N) +
  V * matrix(segments$north_y, nrow = nrow(segments), ncol = N)
qout <- qx * matrix(segments$nx_out, nrow = nrow(segments), ncol = N) +
  qy * matrix(segments$ny_out, nrow = nrow(segments), ncol = N)
qin <- -qout
flux <- qin * matrix(segments$length_m, nrow = nrow(segments), ncol = N)

sector_flux <- function(sec) {
  ii <- which(segments$sector == sec)
  if (!length(ii)) return(rep(NA_real_, N))
  colSums(flux[ii, , drop = FALSE], na.rm = TRUE)
}

west <- sector_flux("west")
east <- sector_flux("east")
south <- sector_flux("south")
north <- sector_flux("north")
net_import <- rowSums(cbind(west, east, south, north), na.rm = TRUE)

boundary_df <- data.frame(
  date = DATES,
  year = as.integer(format(DATES, "%Y")),
  month = as.integer(format(DATES, "%m")),
  west_inflow_kg_s = west,
  east_inflow_kg_s = east,
  south_inflow_kg_s = south,
  north_inflow_kg_s = north,
  net_moisture_import_kg_s = net_import
)

for (nm in names(boundary_df)[4:8]) {
  xa <- rep(NA_real_, N)
  for (m in 1:12) {
    ii <- which(boundary_df$month == m)
    xa[ii] <- boundary_df[[nm]][ii] - mean(boundary_df[[nm]][ii], na.rm = TRUE)
  }
  boundary_df[[paste0(nm, "_anom")]] <- xa
}

write.csv(segments, file.path(OUT_DIR, "circulation_boundary_geometry.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. MOISTURE-FLUX CONVERGENCE
# -----------------------------------------------------------------------------
R <- 6371000
ex <- terra::ext(u)
nc <- terra::ncol(u)
nr <- terra::nrow(u)
dlon <- abs((ex[2] - ex[1]) / nc) * pi / 180
dlat <- abs((ex[4] - ex[3]) / nr) * pi / 180
lat <- seq(ex[4] - (dlat * 180 / pi) / 2, ex[3] + (dlat * 180 / pi) / 2, length.out = nr)
lat_r <- u[[1]]
terra::values(lat_r) <- rep(lat, each = nc)
coslat <- cos(terra::values(lat_r) * pi / 180)
dx <- R * coslat * dlon
dy <- R * dlat

kx <- matrix(c(-0.5, 0, 0.5), 1, 3)
ky <- matrix(c(-0.5, 0, 0.5), 3, 1)
t0 <- stage_start("Computing zonal IVT derivative dU/dx with terra::focal (912 layers)")
dU <- terra::focal(u, kx, fun = sum, na.rm = FALSE, fillvalue = NA)
stage_done(t0, "COMPLETE: dU/dx numerator")

vcos <- v * cos(lat_r * pi / 180)
t0 <- stage_start("Computing meridional IVT derivative d(V cos(lat))/dy with terra::focal (912 layers)")
dVcos <- terra::focal(vcos, ky, fun = sum, na.rm = FALSE, fillvalue = NA)
stage_done(t0, "COMPLETE: dV/dy numerator")

t0 <- stage_start("Combining derivatives to compute moisture-flux convergence (MFC)")
dUdx <- dU / dx
dVdy <- dVcos / dy / pmax(coslat, .Machine$double.eps)
mfc <- -(dUdx + dVdy)
terra::time(mfc) <- terra::time(u)
stage_done(t0, "COMPLETE: MFC field")

mfc_anom <- monthly_anomaly(mfc, label = "MFC")

t0 <- stage_start("Projecting MFC to BC Albers and constructing basin mask/cell areas")
mfc_template <- terra::project(u[[1]], EQUAL_AREA_CRS, method = "near")
basin_mask <- terra::rasterize(basin_eq, mfc_template, field = 1, background = NA)
cell_area <- terra::cellSize(basin_mask, unit = "m")
mfc_eq <- terra::project(mfc, basin_mask, method = "bilinear")
mfc_anom_eq <- terra::project(mfc_anom, basin_mask, method = "bilinear")
stage_done(t0, "COMPLETE: MFC projected to basin grid")

t0 <- stage_start("Computing basin MFC time series and basin-integrated MFC")
mfc_mean <- mean_raster_ts(terra::mask(mfc_eq, basin_mask), "basin-mean MFC")
mfc_mean_anom <- mean_raster_ts(terra::mask(mfc_anom_eq, basin_mask), "basin-mean MFC anomaly")
mfc_int <- mean_raster_ts(terra::mask(mfc_eq * cell_area, basin_mask), "basin-integrated MFC")
mfc_int_anom <- mean_raster_ts(terra::mask(mfc_anom_eq * cell_area, basin_mask), "basin-integrated MFC anomaly")
stage_done(t0, "COMPLETE: basin MFC summaries")

# -----------------------------------------------------------------------------
# 5. BASIN PRECIPITATION ANOMALY
# -----------------------------------------------------------------------------
t0 <- stage_start("Opening monthly precipitation NetCDF and validating time axis")
p <- terra::rast(P_FILE)
pdates <- read_dates(p, "precipitation")
stage_done(t0, "COMPLETE: precipitation input")

t0 <- stage_start("Masking precipitation to Nechako Basin and computing 912 monthly basin means")
p_basin <- terra::mask(p, terra::project(basin, terra::crs(p)), touches = TRUE)
p_raw <- terra::global(p_basin, "mean", na.rm = TRUE)[, 1]
stage_done(t0, "COMPLETE: basin precipitation means")
p_mm_day <- p_raw
if (is.finite(median(p_raw, na.rm = TRUE)) && median(p_raw, na.rm = TRUE) < 0.1) p_mm_day <- p_raw * 1000

p_anom <- rep(NA_real_, N)
for (m in 1:12) {
  ii <- which(as.integer(format(DATES, "%m")) == m)
  p_anom[ii] <- p_mm_day[ii] - mean(p_mm_day[ii], na.rm = TRUE)
}

# -----------------------------------------------------------------------------
# 6. MONTHLY HYDROLOGY/CIRCULATION PRODUCT
# -----------------------------------------------------------------------------
hydro <- data.frame(
  date = DATES,
  year = as.integer(format(DATES, "%Y")),
  month = as.integer(format(DATES, "%m")),
  ivt_mag_anom = ivt_mag_ts,
  ivt_east_anom = ivt_u_ts,
  ivt_north_anom = ivt_v_ts,
  ivt_direction_to_deg = ivt_dir_ts,
  west_inflow_kg_s = west,
  east_inflow_kg_s = east,
  south_inflow_kg_s = south,
  north_inflow_kg_s = north,
  net_moisture_import_kg_s = net_import,
  west_inflow_kg_s_anom = boundary_df$west_inflow_kg_s_anom,
  east_inflow_kg_s_anom = boundary_df$east_inflow_kg_s_anom,
  south_inflow_kg_s_anom = boundary_df$south_inflow_kg_s_anom,
  north_inflow_kg_s_anom = boundary_df$north_inflow_kg_s_anom,
  net_moisture_import_kg_s_anom = boundary_df$net_moisture_import_kg_s_anom,
  mfc_basin_mean_kg_m2_s = mfc_mean,
  mfc_basin_mean_anom_kg_m2_s = mfc_mean_anom,
  mfc_basin_integrated_kg_s = mfc_int,
  mfc_basin_integrated_anom_kg_s = mfc_int_anom,
  precip_basin_mm_day = p_mm_day,
  precip_basin_anom_mm_day = p_anom
)

write.csv(hydro, file.path(OUT_DIR, "circulation_hydrology_monthly.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 7. IVT/PRECIPITATION ASSOCIATION
# -----------------------------------------------------------------------------
assoc_periods <- list(
  "Full record" = hydro,
  "May-October" = hydro[hydro$month %in% 5:10, , drop = FALSE],
  "August-January" = hydro[hydro$month %in% c(8:12, 1), , drop = FALSE]
)
predictors <- c("ivt_mag_anom", "net_moisture_import_kg_s_anom", "mfc_basin_mean_anom_kg_m2_s")
rows <- list(); k <- 0L
for (pn in names(assoc_periods)) {
  d <- assoc_periods[[pn]]
  for (pr in predictors) {
    ok <- is.finite(d[[pr]]) & is.finite(d$precip_basin_anom_mm_day)
    x <- d[[pr]][ok]; y <- d$precip_basin_anom_mm_day[ok]
    k <- k + 1L
    rows[[k]] <- data.frame(
      period = pn,
      predictor = pr,
      n = length(x),
      pearson_r = if (length(x) >= 3 && sd(x) > 0 && sd(y) > 0) cor(x, y, method = "pearson") else NA_real_,
      spearman_rho = if (length(x) >= 3 && sd(x) > 0 && sd(y) > 0) cor(x, y, method = "spearman") else NA_real_
    )
  }
}
association <- dplyr::bind_rows(rows)
write.csv(association, file.path(OUT_DIR, "circulation_ivt_precip_association.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 8. EVENT SUMMARY
# -----------------------------------------------------------------------------
stage_message("Building predefined drought/event windows and El Niño comparison window")
events <- data.frame(
  event = c("2022 drought", "2023 drought", "2024 drought", "2025 drought", "2023-2024 El Nino"),
  start = as.Date(c("2022-05-01", "2023-05-01", "2024-05-01", "2025-05-01", "2023-06-01")),
  end = as.Date(c("2022-10-01", "2023-10-01", "2024-10-01", "2025-10-01", "2024-05-01"))
)

summary_rows <- dplyr::bind_rows(lapply(seq_len(nrow(events)), function(i) {
  ii <- which(hydro$date >= events$start[i] & hydro$date <= events$end[i])
  ubar <- mean(hydro$ivt_east_anom[ii], na.rm = TRUE)
  vbar <- mean(hydro$ivt_north_anom[ii], na.rm = TRUE)
  data.frame(
    event = events$event[i],
    start = events$start[i],
    end = events$end[i],
    n_months = length(ii),
    domain_mean_ivt_mag_anom = mean(hydro$ivt_mag_anom[ii], na.rm = TRUE),
    domain_mean_ivt_east_anom = ubar,
    domain_mean_ivt_north_anom = vbar,
    domain_ivt_direction_to_deg = (atan2(vbar, ubar) * 180 / pi) %% 360,
    basin_west_inflow_anom_kg_s = mean(hydro$west_inflow_kg_s_anom[ii], na.rm = TRUE),
    basin_east_inflow_anom_kg_s = mean(hydro$east_inflow_kg_s_anom[ii], na.rm = TRUE),
    basin_south_inflow_anom_kg_s = mean(hydro$south_inflow_kg_s_anom[ii], na.rm = TRUE),
    basin_north_inflow_anom_kg_s = mean(hydro$north_inflow_kg_s_anom[ii], na.rm = TRUE),
    basin_net_moisture_import_anom_kg_s = mean(hydro$net_moisture_import_kg_s_anom[ii], na.rm = TRUE),
    basin_mfc_anom_kg_m2_s = mean(hydro$mfc_basin_mean_anom_kg_m2_s[ii], na.rm = TRUE),
    basin_integrated_mfc_anom_kg_s = mean(hydro$mfc_basin_integrated_anom_kg_s[ii], na.rm = TRUE),
    basin_precip_anom_mm_day = mean(hydro$precip_basin_anom_mm_day[ii], na.rm = TRUE)
  )
}))

# Preserve existing 14 output and append the new hydrologic columns.
old_file <- file.path(OUT_DIR, "circulation_event_composites.csv")
if (file.exists(old_file)) {
  old <- tryCatch(read.csv(old_file, stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(old) && nrow(old) == nrow(summary_rows)) {
    old$start <- as.Date(old$start)
    old$end <- as.Date(old$end)
    old$n_months <- as.integer(old$n_months)
    summary_rows <- dplyr::left_join(old, summary_rows, by = c("event", "start", "end", "n_months"))
  }
}
write.csv(summary_rows, old_file, row.names = FALSE)

# -----------------------------------------------------------------------------
# 9. VECTOR MAPS
# -----------------------------------------------------------------------------
t0 <- stage_start("Opening Z500 NetCDF and computing Z500 calendar-month anomalies for event maps")
z <- terra::rast(Z_FILE)
z_anom <- monthly_anomaly(z, label = "Z500")
stage_done(t0, "COMPLETE: Z500 anomaly field")

z_event <- lapply(seq_len(nrow(events)), function(i) {
  stage_message(sprintf("  Event composite %d/%d: %s", i, nrow(events), events$event[i]))
  ii <- which(DATES >= events$start[i] & DATES <= events$end[i])
  list(
    z = terra::mean(z_anom[[ii]], na.rm = TRUE),
    u = terra::mean(u_anom[[ii]], na.rm = TRUE),
    v = terra::mean(v_anom[[ii]], na.rm = TRUE),
    ivt = terra::mean(ivt_mag_anom[[ii]], na.rm = TRUE),
    mfc = terra::mean(mfc_anom[[ii]], na.rm = TRUE)
  )
})

# Persist event fields, including MFC, without changing the existing NetCDF
# variable names for the original four fields.
stack_list <- list()
stack_names <- character()
for (i in seq_along(z_event)) {
  tag <- gsub("[^A-Za-z0-9]+", "_", events$event[i])
  stack_list <- c(
    stack_list,
    list(
      z_event[[i]]$z,
      z_event[[i]]$ivt,
      z_event[[i]]$u,
      z_event[[i]]$v,
      z_event[[i]]$mfc
    )
  )
  stack_names <- c(
    stack_names,
    paste0(tag, "_Z500_anom_m"),
    paste0(tag, "_IVTmag_anom"),
    paste0(tag, "_IVTeast_anom"),
    paste0(tag, "_IVTnorth_anom"),
    paste0(tag, "_MFC_anom")
  )
}
stage_message(sprintf("Writing %d event-map fields to circulation_event_composite_maps.nc (NetCDF compression=9)", length(stack_list)))
t0 <- Sys.time()
stack <- do.call(c, stack_list)
names(stack) <- stack_names
terra::writeCDF(stack, file.path(OUT_DIR, "circulation_event_composite_maps.nc"), overwrite = TRUE, compression = 9)
stage_done(t0, "COMPLETE: event composite NetCDF written")

basin_plot <- sf::st_transform(
  sf::st_as_sf(basin),
  sf::st_crs(u)
)

# -----------------------------------------------------------------------------
# Manuscript-quality basin-scale IVT/Z500 maps
# -----------------------------------------------------------------------------
# 14a provides the regional circulation context.  This 14b figure is therefore
# deliberately basin-focused: it preserves the IVT anomaly field and Z500
# contours but uses a sparse, regular 2-D arrow grid so individual arrowheads
# remain visible at publication size.
# -----------------------------------------------------------------------------

# Component-distinctness QA.  Directional arrows are shown only when the two
# retained IVT components contain non-degenerate spatial information.
u_mean_abs <- as.numeric(
  terra::global(
    abs(u - v),
    "mean",
    na.rm = TRUE
  )[1, 1]
)

if (!is.finite(u_mean_abs)) {
  u_mean_abs <- NA_real_
}

u_sd <- as.numeric(
  terra::global(
    u,
    "sd",
    na.rm = TRUE
  )[1, 1]
)

v_sd <- as.numeric(
  terra::global(
    v,
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
      "IVT component QA passed for basin maps: mean |U-V| = %.6g; directional arrows enabled.",
      u_mean_abs
    )
  )
} else {
  stage_message(
    sprintf(
      "IVT component QA failed for basin maps: mean |U-V| = %s; directional arrows disabled.",
      ifelse(
        is.finite(u_mean_abs),
        format(u_mean_abs, scientific = TRUE),
        "NA"
      )
    )
  )
}

basin_bb <- sf::st_bbox(
  basin_plot
)

basin_xrange <- as.numeric(
  basin_bb[c("xmin", "xmax")]
)

basin_yrange <- as.numeric(
  basin_bb[c("ymin", "ymax")]
)

basin_dx <- diff(basin_xrange)
basin_dy <- diff(basin_yrange)

BASIN_XLIM <- c(
  basin_xrange[1] - max(1.5, 0.18 * basin_dx),
  basin_xrange[2] + max(1.5, 0.18 * basin_dx)
)

BASIN_YLIM <- c(
  basin_yrange[1] - max(1.5, 0.18 * basin_dy),
  basin_yrange[2] + max(1.5, 0.18 * basin_dy)
)

all_ivt_values <- unlist(
  lapply(
    z_event,
    function(x) terra::values(x$ivt, mat = FALSE)
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
    z_event,
    function(x) terra::values(x$z, mat = FALSE)
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
  n = 8
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

thin_vector_grid <- function(
    df,
    nx = 8L,
    ny = 6L) {

  if (!nrow(df)) {
    return(df)
  }

  ux <- sort(unique(df$x))
  uy <- sort(unique(df$y))

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

make_basin_arrow_data <- function(d) {

  ar <- d[
    is.finite(d$u) &
      is.finite(d$v) &
      is.finite(d$x) &
      is.finite(d$y),
    ,
    drop = FALSE
  ]

  ar <- ar[
    ar$x >= BASIN_XLIM[1] &
      ar$x <= BASIN_XLIM[2] &
      ar$y >= BASIN_YLIM[1] &
      ar$y <= BASIN_YLIM[2],
    ,
    drop = FALSE
  ]

  ar <- thin_vector_grid(
    ar,
    nx = 8L,
    ny = 6L
  )

  if (!nrow(ar)) {
    return(ar)
  }

  theta <- atan2(
    ar$v,
    ar$u
  )

  # Geographic display length.  The direction is determined only by U/V;
  # this length is deliberately independent of IVT magnitude.
  arrow_length <- max(
    0.15,
    min(
      0.65,
      diff(BASIN_YLIM) / 30
    )
  )

  lat_cos <- pmax(
    cos(ar$y * pi / 180),
    0.25
  )

  ar$xend <- ar$x + arrow_length * cos(theta) / lat_cos
  ar$yend <- ar$y + arrow_length * sin(theta)

  ar$xend <- pmax(
    BASIN_XLIM[1],
    pmin(BASIN_XLIM[2], ar$xend)
  )

  ar$yend <- pmax(
    BASIN_YLIM[1],
    pmin(BASIN_YLIM[2], ar$yend)
  )

  ar
}

make_map <- function(i) {

  a <- terra::as.data.frame(
    z_event[[i]]$z,
    xy = TRUE,
    na.rm = FALSE
  )
  names(a)[3] <- "z500"

  b <- terra::as.data.frame(
    z_event[[i]]$u,
    xy = TRUE,
    na.rm = FALSE
  )
  names(b)[3] <- "u"

  cdf <- terra::as.data.frame(
    z_event[[i]]$v,
    xy = TRUE,
    na.rm = FALSE
  )
  names(cdf)[3] <- "v"

  idf <- terra::as.data.frame(
    z_event[[i]]$ivt,
    xy = TRUE,
    na.rm = FALSE
  )
  names(idf)[3] <- "ivt"

  d <- a |>
    dplyr::left_join(
      b,
      by = c("x", "y")
    ) |>
    dplyr::left_join(
      cdf,
      by = c("x", "y")
    ) |>
    dplyr::left_join(
      idf,
      by = c("x", "y")
    )

  d <- d[
    d$x >= BASIN_XLIM[1] &
      d$x <= BASIN_XLIM[2] &
      d$y >= BASIN_YLIM[1] &
      d$y <= BASIN_YLIM[2],
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
        fill = ivt
      )
    ) +
    ggplot2::geom_contour(
      ggplot2::aes(
        z = z500
      ),
      breaks = Z_BREAKS,
      linewidth = 0.24,
      colour = "grey25"
    ) +
    ggplot2::geom_sf(
      data = basin_plot,
      inherit.aes = FALSE,
      fill = NA,
      colour = "black",
      linewidth = 1.05
    ) +
    ggplot2::coord_sf(
      xlim = BASIN_XLIM,
      ylim = BASIN_YLIM,
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
      title = events$event[i],
      subtitle = "Nechako Basin detail: IVT-magnitude anomaly + Z500 anomaly contours",
      fill = "IVT-magnitude\nanomaly\n(kg m⁻¹ s⁻¹)",
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_bw(
      base_size = 10
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 11
      ),
      plot.subtitle = ggplot2::element_text(
        size = 8.5
      ),
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(
        size = 8.5
      ),
      legend.title = ggplot2::element_text(
        size = 8.5
      ),
      legend.text = ggplot2::element_text(
        size = 8
      )
    )

  if (DIRECTIONAL_IVT_VALID) {
    ar <- make_basin_arrow_data(d)

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
          colour = "black",
          linewidth = 0.42,
          lineend = "round",
          arrow = grid::arrow(
            length = grid::unit(
              0.18,
              "cm"
            ),
            type = "closed"
          )
        )
    }
  }

  p0
}

t0 <- stage_start(
  sprintf(
    "Rendering %d basin-scale IVT/Z500 map panels",
    nrow(events)
  )
)

plots <- lapply(
  seq_len(nrow(events)),
  function(i) {
    stage_message(
      sprintf(
        "  Rendering basin map %d/%d: %s",
        i,
        nrow(events),
        events$event[i]
      )
    )
    make_map(i)
  }
)

pvec <- patchwork::wrap_plots(
  plots,
  ncol = 2,
  guides = "collect"
) +
  patchwork::plot_annotation(
    title = "Nechako Basin IVT and Z500 circulation context",
    subtitle = if (DIRECTIONAL_IVT_VALID) {
      "Black arrows show anomalous IVT direction; arrow length is a display parameter only."
    } else {
      "Directional IVT arrows omitted after component QA."
    },
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

stage_message("  Saving PNG basin-scale circulation-map figure...")
ggplot2::ggsave(
  file.path(
    OUT_DIR,
    "Fig_circulation_IVT_Z500_maps.png"
  ),
  pvec,
  width = 13,
  height = 10,
  dpi = FIG_DPI
)

stage_message("  Saving PDF basin-scale circulation-map figure...")
ggplot2::ggsave(
  file.path(
    OUT_DIR,
    "Fig_circulation_IVT_Z500_maps.pdf"
  ),
  pvec,
  width = 13,
  height = 10
)

stage_done(
  t0,
  "COMPLETE: manuscript-quality basin-scale IVT/Z500 maps"
)

message("14b circulation hydrology diagnostics completed.")
message("Basin area: ", round(basin_area_m2 / 1e6, 1), " km2")
message("Boundary segments: ", nrow(segments))
message("Outputs written to: ", OUT_DIR)

cat("\n============================================================\n")
cat("14b progress-enabled run completed successfully.\n")
cat("Outputs: ", OUT_DIR, "\n", sep = "")
cat("============================================================\n")
