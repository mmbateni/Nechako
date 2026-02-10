
# ======= Nechako Basin CanDCS-M6 (PCIC) Retrieval via THREDDS/NCSS =======
# Packages
library(httr)
library(terra)
library(dplyr)
library(lubridate)
library(glue)

# 0) Paths, region, scenarios
out_dir <- file.path(getwd(), "candcs_m6_nechako")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Same bounding box as your ERA5 code: CDS/NCSS expects [North, West, South, East]
area_for_request <- c(56.5, -128.5, 51.9, -122)  # [N, W, S, E]

# THREDDS NCSS base (PCIC). If your site differs, update base_url accordingly:
# Tip: You can locate dataset names in the PCIC M6 map portal and copy the full THREDDS path.
# e.g., https://data.pacificclimate.org/thredds/ncss/downscaled_cmip6_multi/<dataset>.nc
base_url <- "https://data.pacificclimate.org/thredds/ncss"

# Variables available in CanDCS-M6
m6_vars <- c("tasmin", "tasmax", "pr")

# Scenarios to pull (M6 includes SSP3-7.0 for a subset of models)
scenarios <- c("ssp126","ssp245","ssp370","ssp585")

# Years to request (chunked by decade for manageable downloads)
year_chunks <- split(2015:2100, ceiling(seq_along(2015:2100)/10))

# 1) Helper to construct canonical dataset names used by PCIC M6 portal
# NOTE: Confirm the exact dataset name in the portal and adjust if needed.
m6_dataset_name <- function(var, scenario) {
  glue("{var}_day_MBCn_CanESM5_historical-{scenario}_r1i1p2f1_19500101-21001231_Canada.nc")
}

# 2) Function to build an NCSS (NetCDFSubset) URL for a bbox + time range
build_ncss_url <- function(dataset_id, var, north, west, south, east, t_start, t_end) {
  # NCSS query: ?var=...&north=...&south=...&east=...&west=...&disableProjSubset=on&horizStride=1
  #             &time_start=YYYY-MM-DDT00:00:00Z&time_end=YYYY-MM-DDT23:59:59Z&accept=netcdf
  query <- list(
    var = var,
    north = north, south = south, east = east, west = west,
    disableProjSubset = "on", horizStride = 1,
    time_start = paste0(format(t_start, "%Y-%m-%d"), "T00:00:00Z"),
    time_end   = paste0(format(t_end,   "%Y-%m-%d"), "T23:59:59Z"),
    accept = "netcdf"
  )
  paste0(base_url, "/downscaled_cmip6_multi/", dataset_id, "?", URLencode(paste(names(query), query, sep="=", collapse="&")))
}

# 3) Downloader
download_m6_subset <- function(var, scenario, years_vector) {
  dataset_id <- m6_dataset_name(var, scenario)
  for (yr in years_vector) {
    t_start <- as.Date(sprintf("%d-01-01", yr))
    t_end   <- as.Date(sprintf("%d-12-31", yr))
    ncss_url <- build_ncss_url(dataset_id, var,
                               north = area_for_request[1], west = area_for_request[2],
                               south = area_for_request[3], east = area_for_request[4],
                               t_start = t_start, t_end = t_end)
    out_nc <- file.path(out_dir, glue("M6_{var}_{scenario}_{yr}_Nechako.nc"))
    message("Requesting: ", ncss_url)
    resp <- httr::GET(ncss_url)
    stopifnot(httr::status_code(resp) == 200)
    writeBin(httr::content(resp, "raw"), out_nc)
    message("Saved: ", out_nc)
  }
}

# 4) Run downloads (chunk by decade to limit request size)
for (scen in scenarios) {
  for (v in m6_vars) {
    for (chunk in year_chunks) {
      download_m6_subset(var = v, scenario = scen, years_vector = chunk)
    }
  }
}

# 5) Example: load one file and inspect
example_nc <- list.files(out_dir, pattern = "M6_tasmax_ssp585_2090", full.names = TRUE)[1]
if (!is.na(example_nc)) {
  r <- terra::rast(example_nc)
  if (is.na(terra::crs(r))) terra::crs(r) <- "EPSG:4326"
  print(r)
  plot(r[[1]], main = "CanDCS-M6 tasmax (first day, Nechako bbox)")
}
``
