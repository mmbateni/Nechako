################################################################################
# 14circulation_context.R
################################################################################
################################################################################
# BASIN-BOUNDARY-NORMAL INTEGRATED VAPOUR TRANSPORT (IVT)
# NECHAKO DROUGHT FRAMEWORK
#
# Purpose
#   Calculate the monthly atmospheric-moisture flux normal to the authoritative
#   Nechako basin boundary and integrate it around the complete basin perimeter.
#
# Scientific definition
#   IVT = (U, V), where U and V are the eastward and northward vertically
#   integrated moisture-transport components [kg m-1 s-1].
#
#   For each boundary segment i with outward unit normal n_i and length L_i:
#       q_out,i = IVT_i . n_i
#       F_in,i  = -q_out,i L_i
#
#   The basin-integrated inward moisture transport is approximated by
#       F_net,in = - integral_boundary(IVT . n_out ds)
#                  ~= sum_i F_in,i                              [kg s-1]
#
# Sign convention
#   Positive values = net atmospheric moisture import into the Nechako basin.
#   Negative values = net atmospheric moisture export from the basin.
#
# Spatial / data contract
#   - ERA5-Land/ERA5 IVT U and V are used on their native grid.
#   - The raster is never reprojected or resampled for boundary extraction.
#   - Only the basin vector is projected to EQUAL_AREA_CRS for boundary
#     geometry and segment lengths.
#   - IVT is sampled at boundary-segment midpoints using bilinear interpolation
#     where supported by terra::extract(); nearest-cell extraction is retained
#     as a deterministic fallback.
#   - Boundary segments are <= 10 km long.
#   - Local geographic east/north basis vectors are transformed into the same
#     projected coordinate system as the boundary normals before the dot product.
#
# Outputs -> circulation_context/basin_boundary_ivt/
#   basin_boundary_ivt_monthly.csv
#   basin_boundary_ivt_boundary_geometry.csv
#   basin_boundary_ivt_summary.csv
#   basin_boundary_ivt_event_summary.csv
#   basin_boundary_ivt_qa.csv
#
# Section 14 (circulation-hydrology extension) -> circulation_context/
#   circulation_hydrology_monthly.csv       IVT direction/magnitude anomalies,
#                                           boundary transport, MFC, precipitation
#   circulation_ivt_precip_association.csv  full record, May-Oct, Aug-Jan
#   circulation_boundary_geometry.csv
#   circulation_hydrology_qa.csv
#   circulation_event_composites.csv        new event columns (incl. MFC) joined in
#   circulation_event_composite_maps.nc     overwritten; 5 fields/event incl. MFC
#   Fig_circulation_IVT_Z500_maps.pdf/.png
#
# This script does NOT perform atmospheric-river detection, IVT directional
# persistence, 850/925-hPa moisture flux, moisture-source tracking, or a full
# atmospheric moisture budget.
################################################################################

# ==============================================================================
# 0. WORKING DIRECTORY / SHARED UTILITIES
# ==============================================================================

setwd(Sys.getenv("NECHAKO_WD", "D:/Nechako_Drought/Nechako/"))
source("DROUGHT_ANALYSIS_utils.R")
utils_load_packages(c("terra", "sf", "dplyr", "data.table", "lubridate", "ggplot2", "patchwork"))

if (!dir.exists(WD_PATH)) stop("Working directory not found: ", WD_PATH)
setwd(WD_PATH)

# ==============================================================================
# 1. INPUTS / OUTPUTS / ANALYSIS PERIOD
# ==============================================================================

OUT_DIR <- file.path(WD_PATH, "circulation_context", "basin_boundary_ivt")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

IVT_FILE <- file.path(
  WD_PATH, "monthly_data_direct", "ivt_components_monthly.nc"
)
if (!file.exists(IVT_FILE)) stop("Required IVT input not found: ", IVT_FILE)

# Section 14 inputs and output folder (shared with the 14a circulation products).
HYDRO_OUT_DIR <- file.path(WD_PATH, "circulation_context")
dir.create(HYDRO_OUT_DIR, showWarnings = FALSE, recursive = TRUE)
Z_FILE <- file.path(WD_PATH, "monthly_data_direct", "geopotential_500hPa_monthly.nc")
P_FILE <- file.path(WD_PATH, "monthly_data_direct", "total_precipitation_monthly.nc")
for (f in c(Z_FILE, P_FILE)) {
  if (!file.exists(f)) stop("Required Section 14 input not found: ", f)
}

# The shared utility defines BASIN_SHP as the authoritative Nechako basin
# Shapefile directory. Do not use the legacy KMZ input here.
if (!exists("BASIN_SHP", inherits = TRUE) ||
    !nzchar(BASIN_SHP) || !dir.exists(BASIN_SHP)) {
  stop(
    "Authoritative Nechako basin Shapefile directory not found: ",
    if (exists("BASIN_SHP", inherits = TRUE)) BASIN_SHP else "<BASIN_SHP undefined>"
  )
}

ANALYSIS_START <- as.Date("1950-01-01")
ANALYSIS_END   <- as.Date("2025-12-01")
DATES <- seq(ANALYSIS_START, ANALYSIS_END, by = "month")
N <- length(DATES)
if (N != 912L) stop("Expected 912 monthly observations; found ", N)

SEGMENT_MAX_M <- 10000
EPS_DEG <- 0.01

# ==============================================================================
# 2. PROGRESS HELPERS
# ==============================================================================

stage_start <- function(label) {
  cat(sprintf("[%s] %-70s", format(Sys.time(), "%H:%M:%S"), paste0(label, "...")))
  flush.console()
  Sys.time()
}

stage_done <- function(t0, suffix = "done") {
  elapsed <- as.numeric(Sys.time() - t0, units = "secs")
  cat(sprintf(" %s (%.1f s)\n", suffix, elapsed))
  flush.console()
  invisible(elapsed)
}

# ==============================================================================
# 3. DATE / IVT COMPONENT HELPERS
# ==============================================================================

read_dates <- function(x, label) {
  tt <- terra::time(x)
  if (is.null(tt) || length(tt) != terra::nlyr(x) || anyNA(tt)) {
    stop("Invalid NetCDF time axis: ", label)
  }
  d <- as.Date(tt)
  d <- as.Date(format(d, "%Y-%m-01"))
  if (!identical(d, DATES)) {
    stop(
      label, " does not match the authoritative monthly record 1950-01 to 2025-12."
    )
  }
  d
}

extract_ivt <- function(f) {
  s <- try(terra::sds(f), silent = TRUE)
  if (!inherits(s, "try-error") && length(s) >= 2L) {
    nm <- tolower(names(s))
    east <- grep("eastward|east|^u$|_u$", nm)
    north <- grep("northward|north|^v$|_v$", nm)
    if (length(east) && length(north)) {
      return(list(east = s[[east[1L]]], north = s[[north[1L]]]))
    }
    return(list(east = s[[1L]], north = s[[2L]]))
  }
  
  x <- terra::rast(f)
  nm <- tolower(names(x))
  east <- grep("eastward|east|^u$|_u$", nm)
  north <- grep("northward|north|^v$|_v$", nm)
  if (!length(east) || !length(north)) {
    stop("Could not identify eastward/northward IVT components in ", f)
  }
  list(east = x[[east[1L]]], north = x[[north[1L]]])
}

# ==============================================================================
# 4. LOAD AND VALIDATE IVT
# ==============================================================================

cat("\n============================================================\n")
cat("14c basin-boundary-normal IVT diagnostics\n")
cat("Analysis period: 1950-01 through 2025-12 (912 months)\n")
cat("Positive flux = atmospheric moisture import into the basin\n")
cat("============================================================\n\n")

terra::gdalCache(8192)
terra::terraOptions(progress = 1)

t0 <- stage_start("Opening IVT NetCDF and identifying U/V components")
ivt <- extract_ivt(IVT_FILE)
u <- ivt$east
v <- ivt$north
stage_done(t0, "COMPLETE")

if (!inherits(u, "SpatRaster") || !inherits(v, "SpatRaster")) {
  stop("IVT components must be SpatRaster objects.")
}
if (terra::nlyr(u) != N || terra::nlyr(v) != N) {
  stop("IVT layer count does not equal the expected 912 months.")
}

u_dates <- read_dates(u, "eastward IVT")
v_dates <- read_dates(v, "northward IVT")
if (!identical(u_dates, v_dates)) stop("Eastward and northward IVT time axes differ.")

if (is.na(terra::crs(u)) || !nzchar(terra::crs(u))) stop("Eastward IVT has no CRS.")
if (is.na(terra::crs(v)) || !nzchar(terra::crs(v))) stop("Northward IVT has no CRS.")
if (!terra::same.crs(u, v)) stop("Eastward and northward IVT CRSs differ.")
if (!terra::compareGeom(u[[1]], v[[1]], stopOnError = FALSE,
                        crs = TRUE, ext = TRUE, rowcol = TRUE, res = TRUE)) {
  stop("Eastward and northward IVT spatial geometries differ.")
}

cat(sprintf(
  "IVT grid: %d rows x %d cols x %d months | CRS=%s\n",
  terra::nrow(u), terra::ncol(u), terra::nlyr(u), terra::crs(u)
))

# ==============================================================================
# 5. BASIN GEOMETRY AND BOUNDARY SEGMENTATION
# ==============================================================================

t0 <- stage_start("Loading authoritative Nechako basin and projecting boundary geometry")
basin <- load_basin_vect(BASIN_SHP)
if (!inherits(basin, "SpatVector")) stop("Basin input did not produce a SpatVector.")

basin_eq <- terra::project(basin, EQUAL_AREA_CRS)
basin_area_m2 <- sum(terra::expanse(basin_eq, unit = "m"), na.rm = TRUE)
if (!is.finite(basin_area_m2) || basin_area_m2 <= 0) stop("Invalid basin area.")

bsf <- sf::st_as_sf(basin_eq)
bsf <- sf::st_make_valid(bsf)
bsf_union <- sf::st_union(bsf)

boundary <- sf::st_boundary(bsf_union)
boundary <- sf::st_cast(boundary, "LINESTRING")
boundary <- sf::st_segmentize(boundary, dfMaxLength = SEGMENT_MAX_M)

centroid <- sf::st_coordinates(sf::st_centroid(bsf_union))[1, ]
if (length(centroid) != 2L || any(!is.finite(centroid))) {
  stop("Could not determine a finite basin centroid for outward-normal QA.")
}
stage_done(t0, "COMPLETE")

# Build one row per straight boundary segment. The normal is determined from
# the local segment geometry and then oriented outward using the basin centroid.
segments <- list()
k <- 0L
for (li in seq_len(length(boundary))) {
  xy <- sf::st_coordinates(boundary[li])
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
    
    score1 <- sum(n1 * radial)
    score2 <- sum(n2 * radial)
    if (!is.finite(score1) || !is.finite(score2)) {
      stop("Non-finite outward-normal orientation score at boundary segment.")
    }
    nout <- if (score1 >= score2) n1 else n2
    
    k <- k + 1L
    segments[[k]] <- data.frame(
      segment_id = k,
      line_id = li,
      x_mid_eq = mid[1],
      y_mid_eq = mid[2],
      length_m = L,
      nx_out = nout[1],
      ny_out = nout[2],
      normal_radial_score = max(score1, score2),
      stringsAsFactors = FALSE
    )
  }
}

segments <- dplyr::bind_rows(segments)
if (!nrow(segments)) stop("No basin boundary segments were generated.")

normal_norm <- sqrt(segments$nx_out^2 + segments$ny_out^2)
if (any(!is.finite(normal_norm)) || any(abs(normal_norm - 1) > 1e-8)) {
  stop("Boundary outward normals are not unit vectors.")
}

# Cardinal sector is retained only as a descriptive diagnostic. It is based on
# the dominant component of the local outward normal; it is NOT a replacement
# for the full boundary integral.
segments$sector <- ifelse(
  abs(segments$nx_out) >= abs(segments$ny_out),
  ifelse(segments$nx_out < 0, "west", "east"),
  ifelse(segments$ny_out < 0, "south", "north")
)

# ==============================================================================
# 6. TRANSFORM BOUNDARY MIDPOINTS AND LOCAL EAST/NORTH BASIS
# ==============================================================================

# IVT U/V are geographic eastward/northward components. The boundary normal is
# defined in EQUAL_AREA_CRS. We therefore construct local geographic east and
# north unit vectors at every midpoint and express them in EQUAL_AREA_CRS.

t0 <- stage_start("Constructing local east/north basis vectors at boundary midpoints")

mid_eq <- sf::st_as_sf(
  segments,
  coords = c("x_mid_eq", "y_mid_eq"),
  crs = sf::st_crs(EQUAL_AREA_CRS),
  remove = FALSE
)

mid_native <- sf::st_transform(mid_eq, sf::st_crs(terra::crs(u)))
mid_native_xy <- sf::st_coordinates(mid_native)
mid_eq_xy <- sf::st_coordinates(mid_eq)

if (nrow(mid_native_xy) != nrow(segments)) stop("Boundary midpoint transformation failed.")

make_perturbed <- function(xy, dx = 0, dy = 0) {
  xy2 <- xy
  xy2[, 1] <- xy2[, 1] + dx
  xy2[, 2] <- xy2[, 2] + dy
  sf::st_as_sf(
    data.frame(x = xy2[, 1], y = xy2[, 2]),
    coords = c("x", "y"),
    crs = sf::st_crs(terra::crs(u))
  )
}

east_eq <- sf::st_coordinates(sf::st_transform(
  make_perturbed(mid_native_xy, EPS_DEG, 0),
  sf::st_crs(EQUAL_AREA_CRS)
))
north_eq <- sf::st_coordinates(sf::st_transform(
  make_perturbed(mid_native_xy, 0, EPS_DEG),
  sf::st_crs(EQUAL_AREA_CRS)
))

eb <- east_eq - mid_eq_xy
nb <- north_eq - mid_eq_xy

eb_norm <- sqrt(rowSums(eb^2))
nb_norm <- sqrt(rowSums(nb^2))
if (any(!is.finite(eb_norm)) || any(!is.finite(nb_norm)) ||
    any(eb_norm <= 0) || any(nb_norm <= 0)) {
  stop("Could not construct finite local east/north basis vectors.")
}

eb <- eb / eb_norm
nb <- nb / nb_norm

segments$east_x <- eb[, 1]
segments$east_y <- eb[, 2]
segments$north_x <- nb[, 1]
segments$north_y <- nb[, 2]

# Basis orthogonality is a useful projection QA check. The transformed
# geographic east/north vectors should be close to orthogonal locally.
basis_dot <- segments$east_x * segments$north_x +
  segments$east_y * segments$north_y
segments$east_north_dot <- basis_dot

if (any(!is.finite(basis_dot)) || any(abs(basis_dot) > 0.01)) {
  stop(
    "Local east/north basis vectors are not sufficiently orthogonal; " ,
    "check the IVT CRS and geographic-coordinate interpretation."
  )
}
stage_done(t0, "COMPLETE")

# ==============================================================================
# 7. EXTRACT IVT AT BOUNDARY MIDPOINTS
# ==============================================================================

t0 <- stage_start(sprintf(
  "Extracting U/V at %d boundary-segment midpoints across 912 months",
  nrow(segments)
))

pts <- terra::vect(
  data.frame(x = mid_native_xy[, 1], y = mid_native_xy[, 2]),
  geom = c("x", "y"),
  crs = terra::crs(u)
)

# Bilinear interpolation is preferable for a point located between grid cells.
# terra versions that do not accept method="bilinear" for extract() are handled
# by a deterministic nearest-cell fallback.
extract_component <- function(r, points) {
  out <- try(terra::extract(r, points, ID = FALSE, method = "bilinear"), silent = TRUE)
  if (inherits(out, "try-error")) {
    cat("  bilinear extraction unavailable; using nearest-cell extraction.\n")
    out <- terra::extract(r, points, ID = FALSE)
  }
  as.matrix(out)
}

U <- extract_component(u, pts)
V <- extract_component(v, pts)

if (nrow(U) != nrow(segments) || nrow(V) != nrow(segments)) {
  stop("Boundary extraction size mismatch.")
}
if (ncol(U) != N || ncol(V) != N) {
  stop("Boundary extraction did not return all 912 monthly layers.")
}

stage_done(t0, "COMPLETE")

# ==============================================================================
# 8. BOUNDARY-NORMAL IVT FLUX
# ==============================================================================

# Convert geographic U/V into the EQUAL_AREA_CRS basis, then take the dot
# product with the outward normal. The resulting normal IVT has units
# kg m-1 s-1. Multiplication by segment length gives kg s-1.

t0 <- stage_start("Computing boundary-normal IVT and integrated monthly flux")

EAST_X <- matrix(segments$east_x, nrow = nrow(segments), ncol = N)
EAST_Y <- matrix(segments$east_y, nrow = nrow(segments), ncol = N)
NORTH_X <- matrix(segments$north_x, nrow = nrow(segments), ncol = N)
NORTH_Y <- matrix(segments$north_y, nrow = nrow(segments), ncol = N)
NX_OUT <- matrix(segments$nx_out, nrow = nrow(segments), ncol = N)
NY_OUT <- matrix(segments$ny_out, nrow = nrow(segments), ncol = N)
LENGTH_M <- matrix(segments$length_m, nrow = nrow(segments), ncol = N)

qx <- U * EAST_X + V * NORTH_X
qy <- U * EAST_Y + V * NORTH_Y

qout <- qx * NX_OUT + qy * NY_OUT
qin <- -qout

flux <- qin * LENGTH_M

# Segment-level diagnostics.
qout_abs_weighted <- abs(qout) * LENGTH_M

sector_flux <- function(sec) {
  ii <- which(segments$sector == sec)
  if (!length(ii)) return(rep(NA_real_, N))
  colSums(flux[ii, , drop = FALSE], na.rm = TRUE)
}

west  <- sector_flux("west")
east  <- sector_flux("east")
south <- sector_flux("south")
north <- sector_flux("north")
net_import <- rowSums(cbind(west, east, south, north), na.rm = TRUE)

# Full boundary integral calculated directly from every segment. This is the
# authoritative result; sector totals are only a decomposition of it.
net_import_direct <- colSums(flux, na.rm = TRUE)

if (any(abs(net_import - net_import_direct) > 1e-8 *
        pmax(1, abs(net_import_direct)), na.rm = TRUE)) {
  stop("Sector-summed net flux does not reproduce the direct boundary integral.")
}

stage_done(t0, "COMPLETE")

# ==============================================================================
# 9. MONTHLY OUTPUT
# ==============================================================================

boundary_df <- data.frame(
  date = DATES,
  year = as.integer(format(DATES, "%Y")),
  month = as.integer(format(DATES, "%m")),
  west_inflow_kg_s = west,
  east_inflow_kg_s = east,
  south_inflow_kg_s = south,
  north_inflow_kg_s = north,
  net_moisture_import_kg_s = net_import,
  net_moisture_import_direct_kg_s = net_import_direct,
  stringsAsFactors = FALSE
)

# Calendar-month anomalies relative to the complete 1950-2025 record.
for (nm in c(
  "west_inflow_kg_s", "east_inflow_kg_s", "south_inflow_kg_s",
  "north_inflow_kg_s", "net_moisture_import_kg_s"
)) {
  xa <- rep(NA_real_, N)
  for (m in 1:12) {
    ii <- which(boundary_df$month == m)
    clim <- mean(boundary_df[[nm]][ii], na.rm = TRUE)
    if (is.finite(clim)) xa[ii] <- boundary_df[[nm]][ii] - clim
  }
  boundary_df[[paste0(nm, "_anom")]] <- xa
}

# Useful conversion from kg/s to mean daily atmospheric-water volume is not
# performed here because the scientifically primary quantity is the transport
# rate kg/s. No additional precipitation/unit conversion is introduced.

data.table::fwrite(
  boundary_df,
  file.path(OUT_DIR, "basin_boundary_ivt_monthly.csv")
)

# ==============================================================================
# 10. BOUNDARY GEOMETRY OUTPUT
# ==============================================================================

boundary_geometry <- segments |>
  dplyr::mutate(
    basin_area_m2 = basin_area_m2,
    segment_fraction_of_perimeter = length_m / sum(length_m),
    outward_normal_azimuth_deg = (atan2(nx_out, ny_out) * 180 / pi) %% 360
  )

data.table::fwrite(
  boundary_geometry,
  file.path(OUT_DIR, "basin_boundary_ivt_boundary_geometry.csv")
)

# ==============================================================================
# 11. SUMMARY / QA
# ==============================================================================

summary_df <- data.frame(
  metric = c(
    "analysis_start",
    "analysis_end",
    "n_months",
    "n_boundary_segments",
    "segment_max_m",
    "boundary_perimeter_m",
    "basin_area_m2",
    "mean_net_import_kg_s",
    "sd_net_import_kg_s",
    "min_net_import_kg_s",
    "max_net_import_kg_s",
    "mean_abs_net_import_kg_s",
    "mean_west_inflow_kg_s",
    "mean_east_inflow_kg_s",
    "mean_south_inflow_kg_s",
    "mean_north_inflow_kg_s"
  ),
  value = c(
    as.character(ANALYSIS_START),
    as.character(ANALYSIS_END),
    N,
    nrow(segments),
    SEGMENT_MAX_M,
    sum(segments$length_m),
    basin_area_m2,
    mean(net_import, na.rm = TRUE),
    sd(net_import, na.rm = TRUE),
    min(net_import, na.rm = TRUE),
    max(net_import, na.rm = TRUE),
    mean(abs(net_import), na.rm = TRUE),
    mean(west, na.rm = TRUE),
    mean(east, na.rm = TRUE),
    mean(south, na.rm = TRUE),
    mean(north, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

data.table::fwrite(
  summary_df,
  file.path(OUT_DIR, "basin_boundary_ivt_summary.csv")
)

# QA table: geometry, extraction, and numerical-integral checks. These are
# deliberately explicit so a run can be audited without inspecting the RDS.
qa_df <- data.frame(
  check = c(
    "IVT layer count",
    "IVT U/V time axes identical",
    "IVT U/V geometry identical",
    "Boundary segments generated",
    "Boundary normals unit length",
    "Local east/north basis orthogonality",
    "Boundary midpoint extraction dimensions",
    "Sector sum equals direct boundary integral",
    "All monthly net-import values finite"
  ),
  result = c(
    N == 912L,
    identical(u_dates, v_dates),
    isTRUE(terra::compareGeom(u[[1]], v[[1]], stopOnError = FALSE,
                              crs = TRUE, ext = TRUE, rowcol = TRUE, res = TRUE)),
    nrow(segments) > 0L,
    all(abs(normal_norm - 1) <= 1e-8),
    all(abs(basis_dot) <= 0.01),
    nrow(U) == nrow(segments) && nrow(V) == nrow(segments) &&
      ncol(U) == N && ncol(V) == N,
    all(abs(net_import - net_import_direct) <=
          1e-8 * pmax(1, abs(net_import_direct)), na.rm = TRUE),
    all(is.finite(net_import))
  ),
  details = c(
    sprintf("%d months", N),
    sprintf("%d dates", length(u_dates)),
    sprintf("%d x %d grid geometry", terra::nrow(u), terra::ncol(u)),
    sprintf("%d segments", nrow(segments)),
    sprintf("max |norm-1| = %.3e", max(abs(normal_norm - 1))),
    sprintf("max |east dot north| = %.3e", max(abs(basis_dot))),
    sprintf("U=%d x %d; V=%d x %d", nrow(U), ncol(U), nrow(V), ncol(V)),
    sprintf("max absolute difference = %.6e kg/s", max(abs(net_import - net_import_direct), na.rm = TRUE)),
    sprintf("finite = %d/%d", sum(is.finite(net_import)), length(net_import))
  ),
  stringsAsFactors = FALSE
)

if (!all(qa_df$result)) {
  stop("One or more basin-boundary IVT QA checks failed.")
}

data.table::fwrite(
  qa_df,
  file.path(OUT_DIR, "basin_boundary_ivt_qa.csv")
)

# ==============================================================================
# 12. EVENT SUMMARIES CONSISTENT WITH CIRCULATION-CONTEXT WORK
# ==============================================================================

episodes <- data.frame(
  event = c(
    "2022 drought", "2023 drought", "2024 drought",
    "2025 drought", "2023-2024 El Nino"
  ),
  start = as.Date(c(
    "2022-05-01", "2023-05-01", "2024-05-01",
    "2025-05-01", "2023-06-01"
  )),
  end = as.Date(c(
    "2022-10-01", "2023-10-01", "2024-10-01",
    "2025-10-01", "2024-05-01"
  )),
  stringsAsFactors = FALSE
)

event_rows <- lapply(seq_len(nrow(episodes)), function(i) {
  ii <- which(boundary_df$date >= episodes$start[i] &
                boundary_df$date <= episodes$end[i])
  if (!length(ii)) stop("No boundary-IVT months in event interval: ", episodes$event[i])
  
  data.frame(
    event = episodes$event[i],
    start = episodes$start[i],
    end = episodes$end[i],
    n_months = length(ii),
    mean_net_import_kg_s = mean(boundary_df$net_moisture_import_kg_s[ii], na.rm = TRUE),
    mean_net_import_anom_kg_s = mean(boundary_df$net_moisture_import_kg_s_anom[ii], na.rm = TRUE),
    mean_west_inflow_kg_s = mean(boundary_df$west_inflow_kg_s[ii], na.rm = TRUE),
    mean_east_inflow_kg_s = mean(boundary_df$east_inflow_kg_s[ii], na.rm = TRUE),
    mean_south_inflow_kg_s = mean(boundary_df$south_inflow_kg_s[ii], na.rm = TRUE),
    mean_north_inflow_kg_s = mean(boundary_df$north_inflow_kg_s[ii], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
})

event_summary <- dplyr::bind_rows(event_rows)
data.table::fwrite(
  event_summary,
  file.path(OUT_DIR, "basin_boundary_ivt_event_summary.csv")
)


# ==============================================================================
# 13. REGIONAL EXTENSION — WHOLE BASIN + 7 SUBBASINS
# ==============================================================================
# PURPOSE
#   Extend the authoritative basin-boundary-normal IVT calculation to the
#   seven supplied Nechako subbasins using the same:
#     - native IVT U/V fields,
#     - boundary-midpoint extraction,
#     - local geographic east/north basis transformation,
#     - outward-normal convention,
#     - <=10 km boundary segmentation, and
#     - midpoint boundary integral.
#
#   The whole-basin monthly series remains the authoritative basin result from
#   Sections 5–12 above. This section calculates the seven subbasin series and
#   assembles an 8-panel manuscript-style figure:
#       (a) whole Nechako Basin
#       (b)–(h) seven supplied subbasins
#
# OUTPUTS -> circulation_context/basin_boundary_ivt/
#   regional_basin_boundary_ivt_monthly.csv
#   regional_basin_boundary_ivt_summary.csv
#   Fig_MS_basin_boundary_normal_IVT_8panel.pdf
#   Fig_MS_basin_boundary_normal_IVT_8panel.png
#
# IMPORTANT
#   Positive values = net atmospheric moisture import into the region.
#   The figure uses the raw integrated boundary-normal IVT [kg s-1], matching
#   the existing whole-basin figure. A common y-axis is used across all panels.
# ==============================================================================

t0 <- stage_start(
  "Loading authoritative Nechako regional geometry for 7-subbasin IVT"
)

regional_regions_14 <- load_nechako_regions()

# The canonical utility returns a named list of SpatVector geometries (one
# element per region), not a single multi-feature SpatVector. Accept that
# contract explicitly and validate every returned regional geometry.
if (!is.list(regional_regions_14) || inherits(regional_regions_14, "SpatVector")) {
  stop(
    "load_nechako_regions() returned an unexpected object of class: ",
    paste(class(regional_regions_14), collapse = "/"),
    ". Expected a named list of SpatVector regional geometries."
  )
}

regional_region_ids_14 <- names(regional_regions_14)
if (is.null(regional_region_ids_14) ||
    any(!nzchar(regional_region_ids_14)) ||
    anyDuplicated(regional_region_ids_14)) {
  stop("load_nechako_regions() must return a non-empty uniquely named regional list.")
}

regional_spatvector_check_14 <- vapply(
  regional_regions_14,
  inherits,
  logical(1),
  what = "SpatVector"
)

if (any(!regional_spatvector_check_14)) {
  stop(
    "One or more objects returned by load_nechako_regions() are not SpatVector ",
    "geometries: ",
    paste(names(regional_regions_14)[!regional_spatvector_check_14], collapse = ", ")
  )
}

if (!length(regional_region_ids_14)) {
  stop("No regions were returned by load_nechako_regions().")
}

whole_region_id_14 <- if (exists("NECHAKO_REGION_BASIN_ID", inherits = TRUE)) {
  NECHAKO_REGION_BASIN_ID
} else {
  "Nechako_Basin"
}

if (!whole_region_id_14 %in% regional_region_ids_14) {
  stop(
    "Authoritative whole-basin region '",
    whole_region_id_14,
    "' was not found in load_nechako_regions()."
  )
}

subbasin_ids_14 <- setdiff(regional_region_ids_14, whole_region_id_14)

if (length(subbasin_ids_14) != 7L) {
  stop(
    "Expected exactly 7 supplied subbasins after excluding the whole-basin ",
    "union; found ", length(subbasin_ids_14), "."
  )
}

stage_done(
  t0,
  sprintf("COMPLETE (%d whole-basin + %d subbasin regions)", 1L, length(subbasin_ids_14))
)

# Reuse the same boundary-integral method as the whole-basin calculation above.
# The helper deliberately returns only the final monthly integrated flux, which
# keeps the regional extension memory-efficient relative to retaining all
# segment x month matrices after each region has been processed.
calculate_region_boundary_ivt_14 <- function(region_vect, region_id, region_name) {
  
  if (!inherits(region_vect, "SpatVector")) {
    stop("Region ", region_id, " is not a SpatVector.")
  }
  
  region_eq <- terra::project(region_vect, EQUAL_AREA_CRS)
  region_sf <- sf::st_as_sf(region_eq)
  region_sf <- sf::st_make_valid(region_sf)
  region_union <- sf::st_union(region_sf)
  
  if (length(region_union) == 0L || sf::st_is_empty(region_union)) {
    stop("Region ", region_id, " has empty geometry after union.")
  }
  
  region_area_m2 <- sum(
    sf::st_area(region_union),
    na.rm = TRUE
  )
  region_area_m2 <- as.numeric(region_area_m2)
  
  if (!is.finite(region_area_m2) || region_area_m2 <= 0) {
    stop("Region ", region_id, " has invalid area.")
  }
  
  region_boundary <- sf::st_boundary(region_union)
  region_boundary <- sf::st_cast(region_boundary, "LINESTRING")
  region_boundary <- sf::st_segmentize(
    region_boundary,
    dfMaxLength = SEGMENT_MAX_M
  )
  
  region_centroid <- sf::st_coordinates(
    sf::st_centroid(region_union)
  )[1, ]
  
  if (length(region_centroid) != 2L ||
      any(!is.finite(region_centroid))) {
    stop("Could not determine a finite centroid for region ", region_id, ".")
  }
  
  region_segments <- list()
  kk <- 0L
  
  for (li in seq_len(length(region_boundary))) {
    xy <- sf::st_coordinates(region_boundary[li])
    if (nrow(xy) < 2L) next
    
    for (j in seq_len(nrow(xy) - 1L)) {
      p1 <- xy[j, 1:2]
      p2 <- xy[j + 1L, 1:2]
      
      dx <- p2[1] - p1[1]
      dy <- p2[2] - p1[2]
      L <- sqrt(dx^2 + dy^2)
      
      if (!is.finite(L) || L <= 0) next
      
      mid <- c(
        mean(c(p1[1], p2[1])),
        mean(c(p1[2], p2[2]))
      )
      
      n1 <- c(-dy, dx) / L
      n2 <- c(dy, -dx) / L
      radial <- mid - region_centroid
      
      score1 <- sum(n1 * radial)
      score2 <- sum(n2 * radial)
      
      if (!is.finite(score1) || !is.finite(score2)) {
        stop(
          "Non-finite outward-normal orientation score in region ",
          region_id, "."
        )
      }
      
      nout <- if (score1 >= score2) n1 else n2
      
      kk <- kk + 1L
      region_segments[[kk]] <- data.frame(
        segment_id = kk,
        x_mid_eq = mid[1],
        y_mid_eq = mid[2],
        length_m = L,
        nx_out = nout[1],
        ny_out = nout[2],
        stringsAsFactors = FALSE
      )
    }
  }
  
  seg <- dplyr::bind_rows(region_segments)
  
  if (!nrow(seg)) {
    stop("No boundary segments were generated for region ", region_id, ".")
  }
  
  normal_norm <- sqrt(seg$nx_out^2 + seg$ny_out^2)
  
  if (any(!is.finite(normal_norm)) ||
      any(abs(normal_norm - 1) > 1e-8)) {
    stop("Boundary outward normals failed QA for region ", region_id, ".")
  }
  
  mid_eq_reg <- sf::st_as_sf(
    seg,
    coords = c("x_mid_eq", "y_mid_eq"),
    crs = sf::st_crs(EQUAL_AREA_CRS),
    remove = FALSE
  )
  
  mid_native_reg <- sf::st_transform(
    mid_eq_reg,
    sf::st_crs(terra::crs(u))
  )
  mid_native_xy_reg <- sf::st_coordinates(mid_native_reg)
  mid_eq_xy_reg <- sf::st_coordinates(mid_eq_reg)
  
  if (nrow(mid_native_xy_reg) != nrow(seg)) {
    stop("Boundary midpoint transformation failed for region ", region_id, ".")
  }
  
  make_perturbed_reg <- function(xy, dx = 0, dy = 0) {
    xy2 <- xy
    xy2[, 1] <- xy2[, 1] + dx
    xy2[, 2] <- xy2[, 2] + dy
    
    sf::st_as_sf(
      data.frame(x = xy2[, 1], y = xy2[, 2]),
      coords = c("x", "y"),
      crs = sf::st_crs(terra::crs(u))
    )
  }
  
  east_eq_reg <- sf::st_coordinates(sf::st_transform(
    make_perturbed_reg(mid_native_xy_reg, EPS_DEG, 0),
    sf::st_crs(EQUAL_AREA_CRS)
  ))
  
  north_eq_reg <- sf::st_coordinates(sf::st_transform(
    make_perturbed_reg(mid_native_xy_reg, 0, EPS_DEG),
    sf::st_crs(EQUAL_AREA_CRS)
  ))
  
  eb_reg <- east_eq_reg - mid_eq_xy_reg
  nb_reg <- north_eq_reg - mid_eq_xy_reg
  
  eb_norm_reg <- sqrt(rowSums(eb_reg^2))
  nb_norm_reg <- sqrt(rowSums(nb_reg^2))
  
  if (any(!is.finite(eb_norm_reg)) ||
      any(!is.finite(nb_norm_reg)) ||
      any(eb_norm_reg <= 0) ||
      any(nb_norm_reg <= 0)) {
    stop("Invalid local east/north basis for region ", region_id, ".")
  }
  
  eb_reg <- eb_reg / eb_norm_reg
  nb_reg <- nb_reg / nb_norm_reg
  
  basis_dot_reg <-
    eb_reg[, 1] * nb_reg[, 1] +
    eb_reg[, 2] * nb_reg[, 2]
  
  if (any(!is.finite(basis_dot_reg)) ||
      any(abs(basis_dot_reg) > 0.01)) {
    stop(
      "Local east/north basis vectors failed orthogonality QA for region ",
      region_id, "."
    )
  }
  
  pts_reg <- terra::vect(
    data.frame(
      x = mid_native_xy_reg[, 1],
      y = mid_native_xy_reg[, 2]
    ),
    geom = c("x", "y"),
    crs = terra::crs(u)
  )
  
  U_reg <- extract_component(u, pts_reg)
  V_reg <- extract_component(v, pts_reg)
  
  if (nrow(U_reg) != nrow(seg) ||
      nrow(V_reg) != nrow(seg) ||
      ncol(U_reg) != N ||
      ncol(V_reg) != N) {
    stop("Boundary extraction dimensions failed for region ", region_id, ".")
  }
  
  EAST_X_reg <- matrix(
    eb_reg[, 1], nrow = nrow(seg), ncol = N
  )
  EAST_Y_reg <- matrix(
    eb_reg[, 2], nrow = nrow(seg), ncol = N
  )
  NORTH_X_reg <- matrix(
    nb_reg[, 1], nrow = nrow(seg), ncol = N
  )
  NORTH_Y_reg <- matrix(
    nb_reg[, 2], nrow = nrow(seg), ncol = N
  )
  NX_OUT_reg <- matrix(
    seg$nx_out, nrow = nrow(seg), ncol = N
  )
  NY_OUT_reg <- matrix(
    seg$ny_out, nrow = nrow(seg), ncol = N
  )
  LENGTH_reg <- matrix(
    seg$length_m, nrow = nrow(seg), ncol = N
  )
  
  qx_reg <- U_reg * EAST_X_reg + V_reg * NORTH_X_reg
  qy_reg <- U_reg * EAST_Y_reg + V_reg * NORTH_Y_reg
  
  qout_reg <- qx_reg * NX_OUT_reg + qy_reg * NY_OUT_reg
  qin_reg <- -qout_reg
  flux_reg <- qin_reg * LENGTH_reg
  
  net_import_reg <- colSums(flux_reg, na.rm = TRUE)
  
  if (any(!is.finite(net_import_reg))) {
    stop("Non-finite monthly integrated IVT values for region ", region_id, ".")
  }
  
  data.frame(
    region_id = region_id,
    region_name = region_name,
    date = DATES,
    year = as.integer(format(DATES, "%Y")),
    month = as.integer(format(DATES, "%m")),
    net_moisture_import_kg_s = net_import_reg,
    n_boundary_segments = nrow(seg),
    boundary_perimeter_m = sum(seg$length_m),
    area_m2 = region_area_m2,
    stringsAsFactors = FALSE
  )
}

regional_ivt_list_14 <- list()

for (rid in subbasin_ids_14) {
  cat(sprintf("  Calculating basin-boundary-normal IVT: %s ... ", rid))
  flush.console()
  
  regional_ivt_list_14[[rid]] <- calculate_region_boundary_ivt_14(
    region_vect = regional_regions_14[[rid]],
    region_id = rid,
    region_name = gsub("_", " ", rid, fixed = TRUE)
  )
  
  cat("complete\n")
}

regional_ivt_14 <- dplyr::bind_rows(regional_ivt_list_14)

if (nrow(regional_ivt_14) != 7L * N) {
  stop(
    "Regional IVT output should contain exactly ",
    7L * N,
    " rows; found ",
    nrow(regional_ivt_14), "."
  )
}

# Retain a compact machine-readable regional monthly product.
data.table::fwrite(
  regional_ivt_14,
  file.path(
    OUT_DIR,
    "regional_basin_boundary_ivt_monthly.csv"
  )
)

regional_ivt_summary_14 <- regional_ivt_14 |>
  dplyr::group_by(region_id, region_name) |>
  dplyr::summarise(
    n_months = dplyr::n(),
    n_boundary_segments = dplyr::first(n_boundary_segments),
    boundary_perimeter_m = dplyr::first(boundary_perimeter_m),
    area_m2 = dplyr::first(area_m2),
    mean_net_import_kg_s = mean(
      net_moisture_import_kg_s,
      na.rm = TRUE
    ),
    sd_net_import_kg_s = sd(
      net_moisture_import_kg_s,
      na.rm = TRUE
    ),
    min_net_import_kg_s = min(
      net_moisture_import_kg_s,
      na.rm = TRUE
    ),
    max_net_import_kg_s = max(
      net_moisture_import_kg_s,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# Add the already-computed whole-basin result to the same summary table.
whole_summary_row_14 <- data.frame(
  region_id = whole_region_id_14,
  region_name = "Nechako Basin",
  n_months = N,
  n_boundary_segments = nrow(segments),
  boundary_perimeter_m = sum(segments$length_m),
  area_m2 = basin_area_m2,
  mean_net_import_kg_s = mean(
    boundary_df$net_moisture_import_kg_s,
    na.rm = TRUE
  ),
  sd_net_import_kg_s = sd(
    boundary_df$net_moisture_import_kg_s,
    na.rm = TRUE
  ),
  min_net_import_kg_s = min(
    boundary_df$net_moisture_import_kg_s,
    na.rm = TRUE
  ),
  max_net_import_kg_s = max(
    boundary_df$net_moisture_import_kg_s,
    na.rm = TRUE
  ),
  stringsAsFactors = FALSE
)

regional_ivt_summary_14 <- dplyr::bind_rows(
  whole_summary_row_14,
  regional_ivt_summary_14
)

data.table::fwrite(
  regional_ivt_summary_14,
  file.path(
    OUT_DIR,
    "regional_basin_boundary_ivt_summary.csv"
  )
)

# ------------------------------------------------------------------------------
# 13.1 Assemble the 8-panel 2022–2025 figure.
#
# The whole-basin boundary-integrated moisture transport is intrinsically much
# larger in magnitude than the transport across individual subbasin
# boundaries. A single fixed y-axis range therefore compresses the subbasin
# signals when the whole basin and subbasins are plotted together. The units
# remain identical (kg s^-1), but the whole-basin panel and the seven-subbasin
# panels use separate y-axis ranges. The seven subbasins retain one common
# y-axis range so that their magnitudes remain directly comparable.
# ------------------------------------------------------------------------------

whole_plot_14 <- boundary_df |>
  dplyr::filter(
    date >= as.Date("2022-01-01"),
    date <= as.Date("2025-12-01")
  ) |>
  dplyr::transmute(
    region_id = whole_region_id_14,
    region_name = "Nechako Basin",
    date = date,
    net_moisture_import_kg_s = net_moisture_import_kg_s
  )

sub_plot_14 <- regional_ivt_14 |>
  dplyr::filter(
    date >= as.Date("2022-01-01"),
    date <= as.Date("2025-12-01")
  ) |>
  dplyr::select(
    region_id,
    region_name,
    date,
    net_moisture_import_kg_s
  )

plot_ivt_8panel_14 <- dplyr::bind_rows(
  whole_plot_14,
  sub_plot_14
)

expected_plot_dates_14 <- seq(
  as.Date("2022-01-01"),
  as.Date("2025-12-01"),
  by = "month"
)

plot_coverage_14 <- plot_ivt_8panel_14 |>
  dplyr::group_by(region_id) |>
  dplyr::summarise(
    n_months = dplyr::n(),
    n_expected = sum(date %in% expected_plot_dates_14),
    n_finite = sum(is.finite(net_moisture_import_kg_s)),
    .groups = "drop"
  )

if (nrow(plot_coverage_14) != 8L ||
    any(
      plot_coverage_14$n_months != length(expected_plot_dates_14) |
      plot_coverage_14$n_expected != length(expected_plot_dates_14) |
      plot_coverage_14$n_finite != length(expected_plot_dates_14)
    )) {
  stop(
    "The 8-panel IVT figure requires exactly 48 finite monthly values ",
    "for each of the whole basin and seven subbasins."
  )
}

# Preserve the authoritative regional order: whole basin first, followed by
# the seven supplied subbasins in the order returned by load_nechako_regions().
panel_ids_14 <- c(
  whole_region_id_14,
  subbasin_ids_14
)

panel_names_14 <- c(
  "Nechako Basin",
  gsub("_", " ", subbasin_ids_14, fixed = TRUE)
)

panel_labels_14 <- paste0(
  "(", letters[seq_along(panel_ids_14)], ") ",
  panel_names_14
)

plot_ivt_8panel_14$region_label <- factor(
  panel_labels_14[
    match(plot_ivt_8panel_14$region_id, panel_ids_14)
  ],
  levels = panel_labels_14
)

# Use one y-axis range for the whole basin. The subbasin panels normally use
# one common range so their magnitudes remain directly comparable. The first
# subbasin panel (panel 2 of the figure) is intentionally given its own tighter
# range because its smaller month-to-month variability is otherwise compressed
# by the other subbasins.
whole_plot_values_14 <- whole_plot_14$net_moisture_import_kg_s
sub_plot_values_14 <- sub_plot_14$net_moisture_import_kg_s

ivt_whole_plot_lim_14 <- range(
  pretty(
    range(whole_plot_values_14, finite = TRUE),
    n = 7
  ),
  finite = TRUE
)

ivt_sub_plot_lim_14 <- range(
  pretty(
    range(sub_plot_values_14, finite = TRUE),
    n = 7
  ),
  finite = TRUE
)

# Panel 2 is the first supplied subbasin. Derive its limits from that panel's
# own 2022-2025 monthly values so the actual temporal volatility is visible.
# The same pretty() logic is retained; only the data used to determine the
# limits changes. This does not alter any IVT calculations or source data.
second_panel_region_id_14 <- subbasin_ids_14[1L]
second_panel_values_14 <- sub_plot_14$net_moisture_import_kg_s[
  sub_plot_14$region_id == second_panel_region_id_14
]

ivt_second_panel_lim_14 <- range(
  pretty(
    range(second_panel_values_14, finite = TRUE),
    n = 7
  ),
  finite = TRUE
)

# Keep the date axis consistent across all panels.
ivt_date_breaks_14 <- seq(
  as.Date("2022-01-01"),
  as.Date("2025-07-01"),
  by = "6 months"
)

# Common panel theme. The y-axis range is supplied separately for the whole
# basin and the subbasins below.
ivt_panel_theme_14 <- ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_line(
      colour = "grey90",
      linewidth = 0.3
    ),
    panel.grid.minor = ggplot2::element_line(
      colour = "grey95",
      linewidth = 0.2
    ),
    panel.border = ggplot2::element_rect(
      colour = "grey55",
      fill = NA,
      linewidth = 0.4
    ),
    axis.text.x = ggplot2::element_text(
      angle = 35,
      hjust = 1,
      size = 7.2
    ),
    axis.text.y = ggplot2::element_text(size = 7.2),
    axis.title.x = ggplot2::element_text(size = 9.5),
    axis.title.y = ggplot2::element_text(size = 9.5),
    axis.ticks = ggplot2::element_line(linewidth = 0.3),
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 10,
      hjust = 0
    ),
    plot.margin = ggplot2::margin(3, 5, 3, 5)
  )

make_ivt_panel_14 <- function(df, panel_label, y_lim) {
  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = date,
      y = net_moisture_import_kg_s,
      group = 1
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.45
    ) +
    ggplot2::geom_line(
      linewidth = 0.65,
      colour = "black"
    ) +
    ggplot2::geom_point(
      size = 1.35,
      colour = "black"
    ) +
    ggplot2::scale_x_date(
      breaks = ivt_date_breaks_14,
      date_labels = "%Y-%m",
      limits = c(
        as.Date("2022-01-01"),
        as.Date("2025-12-01")
      ),
      expand = ggplot2::expansion(
        mult = c(0.005, 0.005)
      )
    ) +
    ggplot2::scale_y_continuous(
      limits = y_lim,
      expand = ggplot2::expansion(
        mult = c(0.02, 0.02)
      ),
      labels = function(x) {
        format(
          x,
          scientific = FALSE,
          trim = TRUE,
          big.mark = ""
        )
      }
    ) +
    ggplot2::labs(
      title = panel_label,
      x = "Date",
      y = expression("Net atmospheric moisture import (kg s"^-1*")")
    ) +
    ivt_panel_theme_14
}

# Build the whole-basin panel with its own y-axis range.
p_ivt_whole_14 <- make_ivt_panel_14(
  whole_plot_14,
  panel_labels_14[1],
  ivt_whole_plot_lim_14
)

# Build each subbasin separately. Panel 2 receives its own tighter y-axis;
# the remaining six subbasin panels retain the common subbasin range.
p_ivt_sub_14 <- lapply(
  seq_along(subbasin_ids_14),
  function(i) {
    rid <- subbasin_ids_14[i]
    df_i <- sub_plot_14 |>
      dplyr::filter(region_id == rid)
    
    y_lim_i <- if (i == 1L) {
      ivt_second_panel_lim_14
    } else {
      ivt_sub_plot_lim_14
    }
    
    make_ivt_panel_14(
      df_i,
      panel_labels_14[i + 1L],
      y_lim_i
    )
  }
)

p_ivt_8panel_14 <- patchwork::wrap_plots(
  c(list(p_ivt_whole_14), p_ivt_sub_14),
  ncol = 2
) +
  patchwork::plot_annotation(
    title = "Basin-boundary-normal atmospheric moisture transport",
    subtitle = paste0(
      "Nechako Basin and seven subbasins, January 2022–December 2025; ",
      "positive values indicate net import"
    ),
    caption = paste0(
      "Computed from the boundary integral of IVT normal to each regional ",
      "boundary; Script 14. The whole-basin panel uses a separate y-axis ",
      "range from the seven-subbasin panels; all panels retain kg s^-1 units."
    ),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 13,
        hjust = 0
      ),
      plot.subtitle = ggplot2::element_text(
        size = 9,
        colour = "grey35",
        hjust = 0
      ),
      plot.caption = ggplot2::element_text(
        size = 7.5,
        colour = "grey40",
        hjust = 0
      ),
      plot.margin = ggplot2::margin(5, 8, 5, 8)
    )
  )

for (ext in c("pdf", "png")) {
  out_ivt_fig_14 <- file.path(
    OUT_DIR,
    paste0(
      "Fig_MS_basin_boundary_normal_IVT_8panel.",
      ext
    )
  )
  
  ggplot2::ggsave(
    out_ivt_fig_14,
    p_ivt_8panel_14,
    width = 13.5,
    height = 14.5,
    units = "in",
    dpi = if (ext == "png") 300 else "print",
    device = ext,
    limitsize = FALSE
  )
  
  cat(
    sprintf(
      "  ✓ 8-panel regional IVT figure written: %s\n",
      out_ivt_fig_14
    )
  )
}

cat(
  "  ✓ Regional IVT extension complete: whole basin + ",
  length(subbasin_ids_14),
  " subbasins; 8-panel figure generated.\n",
  sep = ""
)


# ==============================================================================
# 14. CIRCULATION-HYDROLOGY EXTENSION (WHOLE BASIN)
# ==============================================================================
# Adds, on top of the authoritative whole-basin boundary flux of Sections 5-9:
#   14.1 IVT magnitude and direction anomalies (domain means)
#   14.2 Moisture-flux convergence (MFC): basin-mean and basin-integrated
#   14.3 Basin precipitation anomalies
#   14.4 Monthly hydrology/circulation product
#   14.5 IVT / boundary transport / MFC versus precipitation association
#   14.6 Event summaries (including MFC), joined into
#        circulation_event_composites.csv
#   14.7 Event-composite fields (including MFC) written to
#        circulation_event_composite_maps.nc
#   14.8 Basin-scale IVT/Z500 maps
#
# The boundary flux used here is the QA-checked result from Sections 5-9
# (authoritative BASIN_SHP geometry, bilinear midpoint extraction). It is not
# recomputed.
#
# OUTPUTS -> circulation_context/  (same folder as the 14a circulation products)
#   circulation_hydrology_monthly.csv
#   circulation_ivt_precip_association.csv
#   circulation_boundary_geometry.csv
#   circulation_hydrology_qa.csv
#   circulation_event_composites.csv        (new columns joined in / created)
#   circulation_event_composite_maps.nc     (overwritten: 5 fields per event)
#   Fig_circulation_IVT_Z500_maps.pdf/.png
#
# SHARED-FILE NOTE
#   circulation_event_composite_maps.nc is overwritten with five fields per
#   event (Z500, IVT magnitude, IVT east, IVT north, MFC anomalies), using the
#   same layer names as 14a plus "<event>_MFC_anom". circulation_event_
#   composites.csv keeps any existing 14a columns; columns with the same name
#   as a newly computed column are replaced, so re-running is idempotent.
#   If 14a is run AFTER this script, it rewrites both files without the new
#   fields.
# ==============================================================================

cat("\n============================================================\n")
cat("14. Circulation-hydrology extension (whole basin)\n")
cat("============================================================\n")

HYDRO_FIG_DPI <- if (exists("FIG_DPI", inherits = TRUE)) FIG_DPI else 300
MONTH_INDEX <- as.integer(format(DATES, "%m"))

# ------------------------------------------------------------------------------
# 14.0 Helpers
# ------------------------------------------------------------------------------

# Calendar-month anomaly of a 912-layer raster relative to its 1950-2025
# calendar-month mean. Equivalent to terra::tapp() by month.
monthly_anomaly_14 <- function(x, label) {
  if (terra::nlyr(x) != N) {
    stop("Monthly anomaly: ", label, " does not have ", N, " layers.")
  }
  t0 <- stage_start(sprintf("Calendar-month anomaly: %s", label))
  clim <- vector("list", 12L)
  for (mm in 1:12) {
    clim[[mm]] <- terra::mean(x[[which(MONTH_INDEX == mm)]], na.rm = TRUE)
  }
  clim <- do.call(c, clim)
  out <- x - clim[[MONTH_INDEX]]
  terra::time(out) <- DATES
  stage_done(t0, "COMPLETE")
  out
}

# Calendar-month anomaly of a monthly numeric series.
series_anomaly_14 <- function(x) {
  out <- rep(NA_real_, length(x))
  for (mm in 1:12) {
    ii <- which(MONTH_INDEX == mm)
    clim <- mean(x[ii], na.rm = TRUE)
    if (is.finite(clim)) out[ii] <- x[ii] - clim
  }
  out
}

# Unweighted domain mean per layer (same definition as 14a's domain means).
domain_mean_ts_14 <- function(x, label) {
  out <- terra::global(x, "mean", na.rm = TRUE)[, 1]
  if (length(out) != N) stop(label, " returned an invalid series length.")
  out
}

# Direction the vector points TOWARD, in compass degrees clockwise from north
# (0 = northward, 90 = eastward, 180 = southward, 270 = westward).
compass_to_deg_14 <- function(u_comp, v_comp) {
  d <- (atan2(u_comp, v_comp) * 180 / pi) %% 360
  d[!is.finite(d) | (u_comp == 0 & v_comp == 0)] <- NA_real_
  d
}

safe_cor_14 <- function(x, y, method) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) {
    return(NA_real_)
  }
  stats::cor(x[ok], y[ok], method = method)
}

# ------------------------------------------------------------------------------
# 14.1 IVT MAGNITUDE AND DIRECTION ANOMALIES
# ------------------------------------------------------------------------------
# IVT magnitude is sqrt(U^2 + V^2) of the monthly-mean components (as in 14a).
# Domain means are unweighted grid means over the native IVT grid, matching the
# domain-mean definition used in 14a.

u_anom_14 <- monthly_anomaly_14(u, "eastward IVT")
v_anom_14 <- monthly_anomaly_14(v, "northward IVT")

t0 <- stage_start("Computing IVT magnitude sqrt(U^2 + V^2) for 912 months")
ivt_mag_14 <- sqrt(u^2 + v^2)
terra::time(ivt_mag_14) <- DATES
stage_done(t0, "COMPLETE")

ivt_mag_anom_14 <- monthly_anomaly_14(ivt_mag_14, "IVT magnitude")

t0 <- stage_start("Computing domain means of IVT and IVT anomalies")
ivt_east_ts_14       <- domain_mean_ts_14(u, "domain eastward IVT")
ivt_north_ts_14      <- domain_mean_ts_14(v, "domain northward IVT")
ivt_east_anom_ts_14  <- domain_mean_ts_14(u_anom_14, "domain eastward IVT anomaly")
ivt_north_anom_ts_14 <- domain_mean_ts_14(v_anom_14, "domain northward IVT anomaly")
ivt_mag_anom_ts_14   <- domain_mean_ts_14(ivt_mag_anom_14, "domain IVT magnitude anomaly")
stage_done(t0, "COMPLETE")

ivt_direction_ts_14      <- compass_to_deg_14(ivt_east_ts_14, ivt_north_ts_14)
ivt_anom_direction_ts_14 <- compass_to_deg_14(ivt_east_anom_ts_14, ivt_north_anom_ts_14)

# ------------------------------------------------------------------------------
# 14.2 MOISTURE-FLUX CONVERGENCE (MFC)
# ------------------------------------------------------------------------------
# On a regular lon/lat grid:
#   div(IVT) = 1/(R cos(phi)) * [ dU/dlambda + d(V cos(phi))/dphi ]
#   MFC      = -div(IVT)                                   [kg m-2 s-1]
#
# Centred differences via terra::focal(). terra applies the weight matrix to
# the neighbourhood as displayed: columns run west -> east and rows run
# north -> south. Therefore:
#   W_X = [-0.5, 0, +0.5]           gives (U_east  - U_west)  / 2
#   W_Y = [+0.5, 0, -0.5] (column)  gives (Vc_north - Vc_south) / 2
#
# The grid spacing varies with latitude, so dx and cos(phi) are single-layer
# RASTERS (not numeric vectors): terra applies a single-layer raster cell by
# cell to every layer, whereas a numeric vector would be applied layer by layer.
#
# Positive MFC = moisture-flux convergence. By the divergence theorem the
# basin-integrated MFC should track the boundary net import from Section 8;
# this is reported as a QA diagnostic in 14.4.

if (!terra::is.lonlat(u)) {
  stop("MFC calculation requires IVT on a geographic lon/lat grid.")
}

R_EARTH_M <- 6371000
res_deg_14 <- terra::res(u)
dlon_rad_14 <- res_deg_14[1] * pi / 180
dlat_rad_14 <- res_deg_14[2] * pi / 180

lat_r_14 <- terra::init(u[[1]], "y")
coslat_r_14 <- cos(lat_r_14 * pi / 180)
dx_r_14 <- R_EARTH_M * coslat_r_14 * dlon_rad_14   # single-layer raster [m]
dy_m_14 <- R_EARTH_M * dlat_rad_14                 # scalar [m]

W_X_14 <- matrix(c(-0.5, 0, 0.5), nrow = 1)
W_Y_14 <- matrix(c(0.5, 0, -0.5), ncol = 1)

t0 <- stage_start("MFC: zonal derivative dU/dx (terra::focal, 912 layers)")
dU_half_14 <- terra::focal(u, w = W_X_14, fun = "sum", na.rm = FALSE)
dUdx_14 <- dU_half_14 / dx_r_14
stage_done(t0, "COMPLETE")

t0 <- stage_start("MFC: meridional derivative d(V cos(lat))/dy (912 layers)")
dVc_half_14 <- terra::focal(v * coslat_r_14, w = W_Y_14, fun = "sum", na.rm = FALSE)
dVdy_14 <- dVc_half_14 / (dy_m_14 * coslat_r_14)
stage_done(t0, "COMPLETE")

t0 <- stage_start("MFC: combining derivatives")
mfc_14 <- -(dUdx_14 + dVdy_14)
terra::time(mfc_14) <- DATES
names(mfc_14) <- paste0("mfc_", format(DATES, "%Y_%m"))
stage_done(t0, "COMPLETE")
rm(dU_half_14, dVc_half_14, dUdx_14, dVdy_14)

# Area weights on the NATIVE IVT grid (the raster is not reprojected; only the
# basin vector is). Weight = fraction of cell inside basin x true cell area.
basin_native_14 <- terra::aggregate(terra::project(basin, terra::crs(u)))
cover_u_14 <- terra::rasterize(basin_native_14, u[[1]], cover = TRUE)
w_area_u_14 <- cover_u_14 * terra::cellSize(u[[1]], unit = "m")
basin_area_grid_m2_14 <- terra::global(w_area_u_14, "sum", na.rm = TRUE)[1, 1]

w_area_u_14 <- terra::mask(w_area_u_14, mfc_14[[1]])
mfc_area_m2_14 <- terra::global(w_area_u_14, "sum", na.rm = TRUE)[1, 1]
if (!is.finite(mfc_area_m2_14) || mfc_area_m2_14 <= 0) {
  stop("No valid MFC cells inside the basin.")
}
if (mfc_area_m2_14 < 0.99 * basin_area_grid_m2_14) {
  warning(sprintf(
    "MFC is undefined over %.1f%% of the gridded basin area (domain-edge cells).",
    100 * (1 - mfc_area_m2_14 / basin_area_grid_m2_14)
  ), call. = FALSE)
}

t0 <- stage_start("MFC: basin-integrated [kg s-1] and basin-mean [kg m-2 s-1]")
mfc_int_ts_14 <- terra::global(mfc_14 * w_area_u_14, "sum", na.rm = TRUE)[, 1]
stage_done(t0, "COMPLETE")
if (length(mfc_int_ts_14) != N) stop("Basin-integrated MFC series has invalid length.")

mfc_mean_ts_14 <- mfc_int_ts_14 / mfc_area_m2_14
mfc_int_anom_ts_14 <- series_anomaly_14(mfc_int_ts_14)
mfc_mean_anom_ts_14 <- series_anomaly_14(mfc_mean_ts_14)

# Full MFC anomaly field, needed for the event-composite maps.
mfc_anom_14 <- monthly_anomaly_14(mfc_14, "MFC")

# ------------------------------------------------------------------------------
# 14.3 BASIN PRECIPITATION ANOMALY
# ------------------------------------------------------------------------------

t0 <- stage_start("Opening precipitation NetCDF and validating time axis")
p_14 <- terra::rast(P_FILE)
if (terra::nlyr(p_14) != N) stop("Precipitation layer count is not ", N, ".")
p_dates_14 <- read_dates(p_14, "precipitation")
stage_done(t0, "COMPLETE")

basin_p_14 <- terra::aggregate(terra::project(basin, terra::crs(p_14)))
cover_p_14 <- terra::rasterize(basin_p_14, p_14[[1]], cover = TRUE)
w_area_p_14 <- cover_p_14 * terra::cellSize(p_14[[1]], unit = "m")
w_area_p_14 <- terra::mask(w_area_p_14, p_14[[1]])
p_area_m2_14 <- terra::global(w_area_p_14, "sum", na.rm = TRUE)[1, 1]
if (!is.finite(p_area_m2_14) || p_area_m2_14 <= 0) {
  stop("No valid precipitation cells inside the basin.")
}

t0 <- stage_start("Area-weighted basin-mean precipitation for 912 months")
p_raw_14 <- terra::global(p_14 * w_area_p_14, "sum", na.rm = TRUE)[, 1] / p_area_m2_14
stage_done(t0, "COMPLETE")

# ERA5-Land monthly-mean total precipitation is stored in metres of water per
# day. Convert to mm/day when the file units or magnitudes indicate metres.
p_units_14 <- unique(tolower(trimws(terra::units(p_14))))
p_units_14 <- paste(p_units_14[nzchar(p_units_14)], collapse = ";")
p_median_14 <- stats::median(p_raw_14, na.rm = TRUE)
if (grepl("^m($| )", p_units_14) || (is.finite(p_median_14) && p_median_14 < 0.1)) {
  p_mm_day_14 <- p_raw_14 * 1000
  p_conversion_14 <- "x1000 (m/day -> mm/day)"
} else {
  p_mm_day_14 <- p_raw_14
  p_conversion_14 <- "none (assumed mm/day)"
}
cat(sprintf("  Precipitation units in file: '%s'; conversion applied: %s\n",
            p_units_14, p_conversion_14))
p_anom_14 <- series_anomaly_14(p_mm_day_14)

# ------------------------------------------------------------------------------
# 14.4 MONTHLY HYDROLOGY/CIRCULATION PRODUCT + QA
# ------------------------------------------------------------------------------

if (!identical(boundary_df$date, DATES)) {
  stop("boundary_df dates do not match the authoritative monthly record.")
}

hydro_14 <- data.frame(
  date = DATES,
  year = as.integer(format(DATES, "%Y")),
  month = MONTH_INDEX,
  ivt_east_domain_mean = ivt_east_ts_14,
  ivt_north_domain_mean = ivt_north_ts_14,
  ivt_direction_to_deg = ivt_direction_ts_14,
  ivt_mag_anom = ivt_mag_anom_ts_14,
  ivt_east_anom = ivt_east_anom_ts_14,
  ivt_north_anom = ivt_north_anom_ts_14,
  ivt_anom_direction_to_deg = ivt_anom_direction_ts_14,
  west_inflow_kg_s = boundary_df$west_inflow_kg_s,
  east_inflow_kg_s = boundary_df$east_inflow_kg_s,
  south_inflow_kg_s = boundary_df$south_inflow_kg_s,
  north_inflow_kg_s = boundary_df$north_inflow_kg_s,
  net_moisture_import_kg_s = boundary_df$net_moisture_import_kg_s,
  west_inflow_kg_s_anom = boundary_df$west_inflow_kg_s_anom,
  east_inflow_kg_s_anom = boundary_df$east_inflow_kg_s_anom,
  south_inflow_kg_s_anom = boundary_df$south_inflow_kg_s_anom,
  north_inflow_kg_s_anom = boundary_df$north_inflow_kg_s_anom,
  net_moisture_import_kg_s_anom = boundary_df$net_moisture_import_kg_s_anom,
  mfc_basin_mean_kg_m2_s = mfc_mean_ts_14,
  mfc_basin_mean_anom_kg_m2_s = mfc_mean_anom_ts_14,
  mfc_basin_integrated_kg_s = mfc_int_ts_14,
  mfc_basin_integrated_anom_kg_s = mfc_int_anom_ts_14,
  precip_basin_mm_day = p_mm_day_14,
  precip_basin_anom_mm_day = p_anom_14,
  stringsAsFactors = FALSE
)

data.table::fwrite(
  hydro_14,
  file.path(HYDRO_OUT_DIR, "circulation_hydrology_monthly.csv")
)

data.table::fwrite(
  boundary_geometry,
  file.path(HYDRO_OUT_DIR, "circulation_boundary_geometry.csv")
)

# Divergence-theorem consistency: basin-integrated MFC (area integral on the
# IVT grid) versus boundary net import (line integral). They differ through
# discretization, but should be strongly and positively related. A negative
# correlation indicates a sign error.
mfc_vs_boundary_r_14 <- safe_cor_14(mfc_int_ts_14, net_import, "pearson")
mfc_vs_boundary_slope_14 <- tryCatch(
  unname(stats::coef(stats::lm(mfc_int_ts_14 ~ net_import))[2]),
  error = function(e) NA_real_
)
if (is.finite(mfc_vs_boundary_r_14) && mfc_vs_boundary_r_14 < 0) {
  warning("Basin-integrated MFC is negatively correlated with boundary net import; ",
          "check MFC sign conventions.", call. = FALSE)
}

hydro_qa_14 <- data.frame(
  metric = c(
    "mfc_valid_area_fraction_of_gridded_basin",
    "basin_area_gridded_ivt_m2",
    "basin_area_vector_m2",
    "mfc_integrated_vs_boundary_import_pearson_r",
    "mfc_integrated_vs_boundary_import_ols_slope",
    "mean_mfc_integrated_kg_s",
    "mean_boundary_net_import_kg_s",
    "precip_units_in_file",
    "precip_conversion_applied",
    "precip_valid_area_m2"
  ),
  value = c(
    mfc_area_m2_14 / basin_area_grid_m2_14,
    basin_area_grid_m2_14,
    basin_area_m2,
    mfc_vs_boundary_r_14,
    mfc_vs_boundary_slope_14,
    mean(mfc_int_ts_14, na.rm = TRUE),
    mean(net_import, na.rm = TRUE),
    p_units_14,
    p_conversion_14,
    p_area_m2_14
  ),
  stringsAsFactors = FALSE
)

data.table::fwrite(
  hydro_qa_14,
  file.path(HYDRO_OUT_DIR, "circulation_hydrology_qa.csv")
)

# ------------------------------------------------------------------------------
# 14.5 IVT / TRANSPORT / MFC VERSUS PRECIPITATION ASSOCIATION
# ------------------------------------------------------------------------------
# Descriptive correlations only. Monthly anomalies are serially correlated, so
# no naive significance test is reported.

assoc_periods_14 <- list(
  "Full record"    = hydro_14,
  "May-October"    = hydro_14[hydro_14$month %in% 5:10, , drop = FALSE],
  "August-January" = hydro_14[hydro_14$month %in% c(8:12, 1L), , drop = FALSE]
)

assoc_predictors_14 <- c(
  "ivt_mag_anom",
  "net_moisture_import_kg_s_anom",
  "mfc_basin_mean_anom_kg_m2_s"
)

association_14 <- dplyr::bind_rows(lapply(names(assoc_periods_14), function(pn) {
  d <- assoc_periods_14[[pn]]
  dplyr::bind_rows(lapply(assoc_predictors_14, function(pr) {
    data.frame(
      period = pn,
      predictor = pr,
      response = "precip_basin_anom_mm_day",
      n = sum(is.finite(d[[pr]]) & is.finite(d$precip_basin_anom_mm_day)),
      pearson_r = safe_cor_14(d[[pr]], d$precip_basin_anom_mm_day, "pearson"),
      spearman_rho = safe_cor_14(d[[pr]], d$precip_basin_anom_mm_day, "spearman"),
      stringsAsFactors = FALSE
    )
  }))
}))

data.table::fwrite(
  association_14,
  file.path(HYDRO_OUT_DIR, "circulation_ivt_precip_association.csv")
)

# ------------------------------------------------------------------------------
# 14.6 EVENT SUMMARY (INCLUDING MFC) -> circulation_event_composites.csv
# ------------------------------------------------------------------------------
# Uses the same `episodes` windows defined in Section 12.

event_hydro_14 <- dplyr::bind_rows(lapply(seq_len(nrow(episodes)), function(i) {
  ii <- which(hydro_14$date >= episodes$start[i] & hydro_14$date <= episodes$end[i])
  if (!length(ii)) stop("No months in event interval: ", episodes$event[i])
  ubar <- mean(hydro_14$ivt_east_anom[ii], na.rm = TRUE)
  vbar <- mean(hydro_14$ivt_north_anom[ii], na.rm = TRUE)
  data.frame(
    event = episodes$event[i],
    start = episodes$start[i],
    end = episodes$end[i],
    n_months = length(ii),
    domain_mean_ivt_mag_anom = mean(hydro_14$ivt_mag_anom[ii], na.rm = TRUE),
    domain_mean_ivt_east_anom = ubar,
    domain_mean_ivt_north_anom = vbar,
    domain_ivt_anom_direction_to_deg = compass_to_deg_14(ubar, vbar),
    basin_west_inflow_anom_kg_s = mean(hydro_14$west_inflow_kg_s_anom[ii], na.rm = TRUE),
    basin_east_inflow_anom_kg_s = mean(hydro_14$east_inflow_kg_s_anom[ii], na.rm = TRUE),
    basin_south_inflow_anom_kg_s = mean(hydro_14$south_inflow_kg_s_anom[ii], na.rm = TRUE),
    basin_north_inflow_anom_kg_s = mean(hydro_14$north_inflow_kg_s_anom[ii], na.rm = TRUE),
    basin_net_moisture_import_anom_kg_s = mean(hydro_14$net_moisture_import_kg_s_anom[ii], na.rm = TRUE),
    basin_mfc_mean_anom_kg_m2_s = mean(hydro_14$mfc_basin_mean_anom_kg_m2_s[ii], na.rm = TRUE),
    basin_mfc_integrated_anom_kg_s = mean(hydro_14$mfc_basin_integrated_anom_kg_s[ii], na.rm = TRUE),
    basin_precip_anom_mm_day = mean(hydro_14$precip_basin_anom_mm_day[ii], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

EVENT_COMPOSITES_FILE_14 <- file.path(HYDRO_OUT_DIR, "circulation_event_composites.csv")
event_keys_14 <- c("event", "start", "end")

if (file.exists(EVENT_COMPOSITES_FILE_14)) {
  old_events_14 <- data.table::fread(EVENT_COMPOSITES_FILE_14, data.table = FALSE)
  if (!all(event_keys_14 %in% names(old_events_14))) {
    stop("Existing circulation_event_composites.csv lacks event/start/end columns; ",
         "refusing to join.")
  }
  old_events_14$start <- as.Date(old_events_14$start)
  old_events_14$end <- as.Date(old_events_14$end)
  
  old_key_14 <- paste(old_events_14$event, old_events_14$start, old_events_14$end)
  new_key_14 <- paste(event_hydro_14$event, event_hydro_14$start, event_hydro_14$end)
  if (!setequal(old_key_14, new_key_14)) {
    warning("Event windows in the existing circulation_event_composites.csv differ ",
            "from Section 12; unmatched events are kept with NA columns.", call. = FALSE)
  }
  
  # Replace same-named columns with the newly computed values (idempotent).
  replaced_14 <- setdiff(intersect(names(old_events_14), names(event_hydro_14)),
                         event_keys_14)
  old_events_14 <- old_events_14[, setdiff(names(old_events_14), replaced_14), drop = FALSE]
  
  event_out_14 <- dplyr::full_join(old_events_14, event_hydro_14, by = event_keys_14)
  cat(sprintf("  Joined %d new/updated column(s) into existing %s\n",
              ncol(event_hydro_14) - length(event_keys_14),
              basename(EVENT_COMPOSITES_FILE_14)))
} else {
  event_out_14 <- event_hydro_14
  cat(sprintf("  %s not found; created it with the hydrology event columns.\n",
              basename(EVENT_COMPOSITES_FILE_14)))
}

lead_cols_14 <- c(event_keys_14, "n_months")
event_out_14 <- event_out_14[, c(lead_cols_14, setdiff(names(event_out_14), lead_cols_14)), drop = FALSE]

data.table::fwrite(event_out_14, EVENT_COMPOSITES_FILE_14)

# ------------------------------------------------------------------------------
# 14.7 EVENT-COMPOSITE FIELDS -> circulation_event_composite_maps.nc
# ------------------------------------------------------------------------------

t0 <- stage_start("Opening Z500 NetCDF and validating time axis")
z_14 <- terra::rast(Z_FILE)
if (terra::nlyr(z_14) != N) stop("Z500 layer count is not ", N, ".")
z_dates_14 <- read_dates(z_14, "Z500")
stage_done(t0, "COMPLETE")

z_median_14 <- stats::median(terra::values(z_14[[1]], mat = FALSE), na.rm = TRUE)
if (is.finite(z_median_14) && z_median_14 > 20000) {
  warning("Z500 values look like geopotential (m2 s-2), not height (m); ",
          "layer names keep the '_m' suffix used by 14a.", call. = FALSE)
}

z_anom_14 <- monthly_anomaly_14(z_14, "Z500")
z_same_grid_14 <- isTRUE(terra::compareGeom(
  z_14[[1]], u[[1]], stopOnError = FALSE,
  crs = TRUE, ext = TRUE, rowcol = TRUE, res = TRUE
))
if (!z_same_grid_14) {
  cat("  Z500 grid differs from IVT grid; Z500 event composites are bilinearly\n",
      "  resampled to the IVT grid for stacking and mapping only.\n", sep = "")
}

t0 <- stage_start(sprintf("Event composites for %d windows (Z500, IVT, MFC)", nrow(episodes)))
event_fields_14 <- lapply(seq_len(nrow(episodes)), function(i) {
  ii <- which(DATES >= episodes$start[i] & DATES <= episodes$end[i])
  zc <- terra::mean(z_anom_14[[ii]], na.rm = TRUE)
  if (!z_same_grid_14) zc <- terra::resample(zc, u[[1]], method = "bilinear")
  list(
    z   = zc,
    ivt = terra::mean(ivt_mag_anom_14[[ii]], na.rm = TRUE),
    u   = terra::mean(u_anom_14[[ii]], na.rm = TRUE),
    v   = terra::mean(v_anom_14[[ii]], na.rm = TRUE),
    mfc = terra::mean(mfc_anom_14[[ii]], na.rm = TRUE)
  )
})
stage_done(t0, "COMPLETE")

stack_list_14 <- list()
stack_names_14 <- character()
for (i in seq_along(event_fields_14)) {
  tag <- gsub("[^A-Za-z0-9]+", "_", episodes$event[i])
  ef <- event_fields_14[[i]]
  stack_list_14 <- c(stack_list_14, list(ef$z, ef$ivt, ef$u, ef$v, ef$mfc))
  stack_names_14 <- c(
    stack_names_14,
    paste0(tag, "_Z500_anom_m"),
    paste0(tag, "_IVTmag_anom"),
    paste0(tag, "_IVTeast_anom"),
    paste0(tag, "_IVTnorth_anom"),
    paste0(tag, "_MFC_anom")
  )
}

t0 <- stage_start("Writing circulation_event_composite_maps.nc (overwrite)")
# One named NetCDF variable per field and event, so the layer names survive a
# round trip (a single multi-layer variable does not keep terra layer names).
comp_stack_14 <- terra::sds(stack_list_14)
names(comp_stack_14) <- stack_names_14
terra::writeCDF(
  comp_stack_14,
  file.path(HYDRO_OUT_DIR, "circulation_event_composite_maps.nc"),
  overwrite = TRUE,
  compression = 9
)
stage_done(t0, "COMPLETE")

# ------------------------------------------------------------------------------
# 14.8 BASIN-SCALE IVT/Z500 MAPS -> Fig_circulation_IVT_Z500_maps.*
# ------------------------------------------------------------------------------
# Raster = IVT-magnitude anomaly; contours = Z500 anomaly; arrows = anomalous
# IVT direction on a sparse regular 2-D grid (arrow length is display-only).

basin_plot_14 <- sf::st_as_sf(basin_native_14)

# Component QA: arrows are drawn only if the U and V composites are distinct.
uv_diff_14 <- mean(vapply(event_fields_14, function(ef) {
  terra::global(abs(ef$u - ef$v), "mean", na.rm = TRUE)[1, 1]
}, numeric(1)), na.rm = TRUE)
DIRECTIONAL_IVT_VALID_14 <- is.finite(uv_diff_14) && uv_diff_14 > 1e-10
cat(sprintf("  IVT component QA (mean |U-V| of composites = %s): arrows %s\n",
            format(uv_diff_14, digits = 4),
            if (DIRECTIONAL_IVT_VALID_14) "enabled" else "disabled"))

bb_14 <- sf::st_bbox(basin_plot_14)
bx_14 <- as.numeric(bb_14[c("xmin", "xmax")])
by_14 <- as.numeric(bb_14[c("ymin", "ymax")])
BASIN_XLIM_14 <- bx_14 + c(-1, 1) * max(1.5, 0.18 * diff(bx_14))
BASIN_YLIM_14 <- by_14 + c(-1, 1) * max(1.5, 0.18 * diff(by_14))

symmetric_limit_14 <- function(field) {
  vals <- unlist(lapply(event_fields_14, function(ef) terra::values(ef[[field]], mat = FALSE)),
                 use.names = FALSE)
  vals <- vals[is.finite(vals)]
  lim <- if (length(vals)) max(abs(vals)) else 1
  if (!is.finite(lim) || lim <= 0) 1 else lim
}
IVT_LIMIT_14 <- symmetric_limit_14("ivt")
Z_LIMIT_14 <- symmetric_limit_14("z")
Z_BREAKS_14 <- pretty(c(-Z_LIMIT_14, Z_LIMIT_14), n = 8)
Z_BREAKS_14 <- sort(unique(c(Z_BREAKS_14[abs(Z_BREAKS_14) <= Z_LIMIT_14], 0)))

thin_vector_grid_14 <- function(df, nx = 8L, ny = 6L) {
  if (!nrow(df)) return(df)
  ux <- sort(unique(df$x))
  uy <- sort(unique(df$y))
  ix <- unique(round(seq(1, length(ux), length.out = min(nx, length(ux)))))
  iy <- unique(round(seq(1, length(uy), length.out = min(ny, length(uy)))))
  df[df$x %in% ux[ix] & df$y %in% uy[iy], , drop = FALSE]
}

make_arrow_data_14 <- function(d) {
  ar <- d[is.finite(d$u) & is.finite(d$v), , drop = FALSE]
  ar <- thin_vector_grid_14(ar)
  if (!nrow(ar)) return(ar)
  theta <- atan2(ar$v, ar$u)
  arrow_len <- max(0.15, min(0.65, diff(BASIN_YLIM_14) / 30))
  lat_cos <- pmax(cos(ar$y * pi / 180), 0.25)
  ar$xend <- pmin(BASIN_XLIM_14[2], pmax(BASIN_XLIM_14[1], ar$x + arrow_len * cos(theta) / lat_cos))
  ar$yend <- pmin(BASIN_YLIM_14[2], pmax(BASIN_YLIM_14[1], ar$y + arrow_len * sin(theta)))
  ar
}

make_map_14 <- function(i) {
  ef <- event_fields_14[[i]]
  d <- terra::as.data.frame(c(ef$z, ef$u, ef$v, ef$ivt), xy = TRUE, na.rm = FALSE)
  names(d) <- c("x", "y", "z500", "u", "v", "ivt")
  d <- d[d$x >= BASIN_XLIM_14[1] & d$x <= BASIN_XLIM_14[2] &
           d$y >= BASIN_YLIM_14[1] & d$y <= BASIN_YLIM_14[2], , drop = FALSE]
  
  p0 <- ggplot2::ggplot(d, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_raster(ggplot2::aes(fill = ivt)) +
    ggplot2::geom_contour(ggplot2::aes(z = z500), breaks = Z_BREAKS_14,
                          linewidth = 0.24, colour = "grey25") +
    ggplot2::geom_sf(data = basin_plot_14, inherit.aes = FALSE,
                     fill = NA, colour = "black", linewidth = 1.05) +
    ggplot2::coord_sf(xlim = BASIN_XLIM_14, ylim = BASIN_YLIM_14, expand = FALSE) +
    ggplot2::scale_fill_gradient2(
      low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 0,
      limits = c(-IVT_LIMIT_14, IVT_LIMIT_14), oob = scales::squish,
      na.value = "white"
    ) +
    ggplot2::labs(
      title = episodes$event[i],
      subtitle = "IVT-magnitude anomaly + Z500 anomaly contours",
      fill = "IVT-magnitude\nanomaly\n(kg m\u207b\u00b9 s\u207b\u00b9)",
      x = NULL, y = NULL
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 11),
      plot.subtitle = ggplot2::element_text(size = 8.5),
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 8.5),
      legend.title = ggplot2::element_text(size = 8.5),
      legend.text = ggplot2::element_text(size = 8)
    )
  
  if (DIRECTIONAL_IVT_VALID_14) {
    ar <- make_arrow_data_14(d)
    if (nrow(ar)) {
      p0 <- p0 + ggplot2::geom_segment(
        data = ar,
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
        inherit.aes = FALSE, colour = "black", linewidth = 0.42, lineend = "round",
        arrow = grid::arrow(length = grid::unit(0.18, "cm"), type = "closed")
      )
    }
  }
  p0
}

t0 <- stage_start(sprintf("Rendering %d basin-scale IVT/Z500 map panels", nrow(episodes)))
map_panels_14 <- lapply(seq_len(nrow(episodes)), make_map_14)

p_maps_14 <- patchwork::wrap_plots(map_panels_14, ncol = 2, guides = "collect") +
  patchwork::plot_annotation(
    title = "Nechako Basin IVT and Z500 circulation context",
    subtitle = if (DIRECTIONAL_IVT_VALID_14) {
      "Black arrows show anomalous IVT direction; arrow length is a display parameter only."
    } else {
      "Directional IVT arrows omitted after component QA."
    },
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 15),
      plot.subtitle = ggplot2::element_text(size = 9.5),
      legend.position = "right"
    )
  )

ggplot2::ggsave(file.path(HYDRO_OUT_DIR, "Fig_circulation_IVT_Z500_maps.png"),
                p_maps_14, width = 13, height = 10, dpi = HYDRO_FIG_DPI)
ggplot2::ggsave(file.path(HYDRO_OUT_DIR, "Fig_circulation_IVT_Z500_maps.pdf"),
                p_maps_14, width = 13, height = 10)
stage_done(t0, "COMPLETE")

# Release the large 912-layer intermediates.
rm(u_anom_14, v_anom_14, ivt_mag_14, ivt_mag_anom_14, mfc_14, mfc_anom_14,
   z_14, z_anom_14, p_14)
invisible(gc(verbose = FALSE))

cat("  \u2713 Circulation-hydrology extension complete.\n")


# ==============================================================================
# 15. FINAL REPORT
# ==============================================================================

cat("\n============================================================\n")
cat("14c basin-boundary-normal IVT completed successfully\n")
cat(sprintf("  Period:             %s to %s\n", ANALYSIS_START, ANALYSIS_END))
cat(sprintf("  Boundary segments:  %d\n", nrow(segments)))
cat(sprintf("  Boundary perimeter: %.1f km\n", sum(segments$length_m) / 1000))
cat(sprintf("  Basin area:         %.3f km2\n", basin_area_m2 / 1e6))
cat(sprintf("  Mean net import:    %.3f kg/s\n", mean(net_import, na.rm = TRUE)))
cat(sprintf("  SD net import:      %.3f kg/s\n", sd(net_import, na.rm = TRUE)))
cat(sprintf("  Output directory:   %s\n", normalizePath(OUT_DIR, winslash = "/", mustWork = FALSE)))
cat("  Regional 8-panel:   Fig_MS_basin_boundary_normal_IVT_8panel.pdf/.png\n")
cat(sprintf("  MFC vs boundary import (Pearson r): %.3f\n", mfc_vs_boundary_r_14))
cat(sprintf("  Hydrology outputs:  %s\n", normalizePath(HYDRO_OUT_DIR, winslash = "/", mustWork = FALSE)))
cat("============================================================\n")

invisible(TRUE)

################################################################################
# COMPLEXITY / ARCHITECTURAL NOTES
#
# Time complexity:
#   O(S * T) for boundary extraction and flux calculation, where S is the
#   number of boundary segments and T=912 monthly layers. Raster extraction is
#   additionally dependent on the native IVT grid size and terra's internal
#   block processing.
#
# Space complexity:
#   O(S*T) because the boundary U/V matrices and transformed flux matrices are
#   held in memory. This is substantially smaller than materializing the full
#   912-layer IVT raster when S << number of raster cells.
#
# Numerical method:
#   The boundary integral is a midpoint-rule discretization of
#   -integral(IVT . n_out ds). Segmentization is capped at 10 km. The result is
#   therefore a numerical line integral, not an analytic continuous integral.
#
# Important edge cases:
#   - Invalid/missing IVT time axes stop the run.
#   - U/V geometry or CRS mismatch stops the run.
#   - Empty/invalid basin geometry stops the run.
#   - Non-unit outward normals stop the run.
#   - Failure of local east/north basis construction stops the run.
#   - Bilinear point extraction is attempted first; nearest-cell extraction is a
#     deterministic compatibility fallback.
#   - A concave or multipart basin is handled after polygon union; holes, if
#     present, remain genuine polygon boundaries and are therefore included in
#     the mathematical boundary integral.
#
# Architectural trade-offs:
#   - Boundary vectors are retained in an equal-area projected CRS for stable
#     segment geometry, while IVT components remain native geographic U/V.
#     Local east/north basis transformation reconciles these representations.
#   - The full boundary integral is authoritative; west/east/south/north sector
#     totals are descriptive summaries only.
#   - Section 14 adds MFC and precipitation diagnostics without changing the
#     boundary-flux calculation; the boundary flux is reused, not recomputed.
#     No moisture-budget closure is imposed; basin-integrated MFC vs boundary
#     import is reported only as a consistency check.
################################################################################