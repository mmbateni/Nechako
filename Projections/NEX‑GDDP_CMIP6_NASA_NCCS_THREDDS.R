
# ======= Nechako Basin NEX-GDDP CMIP6 (NASA THREDDS NetCDFSubset) =======
# Packages
library(httr)
library(terra)
library(dplyr)
library(lubridate)
library(glue)

out_dir <- file.path(getwd(), "nex_gddp_canesm5_nechako")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

area_for_request <- c(56.5, -128.5, 51.9, -122)  # [North, West, South, East]

# NASA NCCS THREDDS NCSS base (global downscaled, daily, 0.25°)
# Docs + example netcdfSubset usage: https://www.nccs.nasa.gov/services/data-collections/land-based-products/nex-gddp-cmip6
ncss_base <- "https://ds.nccs.nasa.gov/thredds/ncss/AMES/NEX/GDDP-CMIP6"

# Variables typically available in NEX-GDDP CMIP6
nex_vars <- c("tas","tasmin","tasmax","pr","hurs","huss","rlds","rsds","sfcWind")

# Scenarios (adjust if your application prefers a subset)
scenarios <- c("ssp245","ssp370","ssp585")

# Years to request (chunk by year—NEX files are often organized per year)
years <- 2015:2100

# 1) Build dataset path; NEX-GDDP uses per-variable/per-year files:
# Path pattern: /CanESM5/<scenario>/r1i1p1f1/<var>/<var>_day_CanESM5_<scenario>_r1i1p1f1_gn_<YYYY>.nc
nex_dataset_path <- function(var, scenario, year) {
  fname <- glue("{var}_day_CanESM5_{scenario}_r1i1p1f1_gn_{year}.nc")
  glue("{ncss_base}/CanESM5/{scenario}/r1i1p1f1/{var}/{fname}")
}

# 2) Build NetCDFSubset URL for a bbox (lat/lon) and full year
build_ncss_bbox_url <- function(dataset_path, north, west, south, east) {
  query <- list(
    var = sub("_day.*$", "", basename(dataset_path)),    # 'tas','pr', etc.
    north = north, south = south, east = east, west = west,
    disableProjSubset = "on", horizStride = 1,
    # Full year; NCSS ignores time if not specified for full file, but include for safety:
    time_start = paste0(format(as.Date(paste0(substr(basename(dataset_path), nchar(basename(dataset_path))-7, nchar(basename(dataset_path))-4), "-01-01")), "%Y-%m-%d"), "T00:00:00Z"),
    time_end   = paste0(format(as.Date(paste0(substr(basename(dataset_path), nchar(basename(dataset_path))-7, nchar(basename(dataset_path))-4), "-12-31")), "%Y-%m-%d"), "T23:59:59Z"),
    accept = "netcdf"
  )
  paste0(dataset_path, "?", URLencode(paste(names(query), query, sep="=", collapse="&")))
}

# 3) Downloader
download_nex_subset <- function(var, scenario, yr) {
  dataset_path <- nex_dataset_path(var, scenario, yr)
  ncss_url <- build_ncss_bbox_url(dataset_path,
                                  north = area_for_request[1], west = area_for_request[2],
                                  south = area_for_request[3], east = area_for_request[4])
  out_nc <- file.path(out_dir, glue("NEX_{var}_{scenario}_{yr}_Nechako.nc"))
  message("Requesting: ", ncss_url)
  resp <- httr::GET(ncss_url)
  if (httr::status_code(resp) != 200) {
    warning("Failed: ", httr::status_code(resp), " for ", ncss_url)
    return(invisible(NULL))
  }
  writeBin(httr::content(resp, "raw"), out_nc)
  message("Saved: ", out_nc)
}

# 4) Run downloads (filter variables if needed)
for (scen in scenarios) {
  for (v in nex_vars) {
    for (yr in years) {
      download_nex_subset(var = v, scenario = scen, yr = yr)
    }
  }
}

# 5) Example: load one file and inspect
example_nc <- list.files(out_dir, pattern = "NEX_pr_ssp585_2090", full.names = TRUE)[1]
if (!is.na(example_nc)) {
  r <- terra::rast(example_nc)
  if (is.na(terra::crs(r))) terra::crs(r) <- "EPSG:4326"
  print(r)
  plot(r[[1]], main = "NEX-GDDP CMIP6 pr (first day, Nechako bbox)")
}
