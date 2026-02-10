
# =============================================================================
# Fast ERA5-Land 'swvl4' via Earth Data Hub (DestinE) using PAT.
# Streams requested bbox & months (1950–2024); writes monthly NetCDF.
# Uses NumPy chunk manager (no Dask dependency).
# =============================================================================

suppressPackageStartupMessages({
  library(reticulate)
  library(lubridate)
})

# ---- 1) Python environment ---------------------------------------------------
# Create/activate a conda env and install required Python packages (Windows-safe)
setup_py_env <- function(env_name = "era5") {
  conda_bin <- reticulate::conda_binary()
  if (is.na(conda_bin) || !nzchar(conda_bin)) {
    reticulate::install_miniconda()  # one-time
  }
  envs <- reticulate::conda_list()
  if (!(env_name %in% envs$name)) {
    reticulate::conda_create(env_name, packages = "python=3.11")
  }
  reticulate::use_condaenv(env_name, required = TRUE)
  
  py <- reticulate::conda_python(env_name)
  # Upgrade pip & install scientific + HTTP/Zarr stack (no Dask)
  system2(py, c("-m", "pip", "install", "--upgrade", "pip"))
  system2(py, c("-m", "pip", "install",
                "numpy>=1.26",
                "xarray==2024.7.0",
                "zarr==2.17.0",
                "numcodecs>=0.12.1",
                "fsspec>=2024.2.0",
                "aiohttp", "requests"))
  # Sanity checks
  stopifnot(py_module_available("xarray"))
  stopifnot(py_module_available("zarr"))
  stopifnot(py_module_available("fsspec"))
  stopifnot(py_module_available("aiohttp"))
  stopifnot(py_module_available("requests"))
}

setup_py_env("era5")
message("✅ Python bound to: ", reticulate::py_config()$python)

# ---- 2) User settings --------------------------------------------------------
years   <- 1994:2024
months  <- 1:12

# BBox (latitude N/S, longitude W/E):
latN <- 56.5; latS <- 51.9
lonW <- -128.5; lonE <- -122.0

out_dir <- "D:/Nechako_Drought/ERA5-Land/swvl4_1950_2024"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- 3) EDH dataset + authentication ----------------------------------------
# EDH ERA5-Land hourly Zarr (0.1°). Variables include 'swvl4'. (EDH dataset page)
# PAT method documented by EDH: https://edh:<PAT>@host/path (EDH "Getting started")
# Sources: EDH getting-started & dataset page; ERA5-Land documentation.               [1][2][3]
EDH_PAT <- "edh_pat_f628161a6b72e1620777ce8cb0b656975623f0c71b160cb28cd422b420e067ac95da82ec9cdaa7ef43353ffa295854b1"  # << your PAT
EDH_URL <- sprintf("https://edh:%s@data.earthdatahub.destine.eu/era5/reanalysis-era5-land-no-antartica-v0.zarr",
                   EDH_PAT)

# ---- 4) Import Python modules -----------------------------------------------
xr       <- reticulate::import("xarray")
builtins <- reticulate::import_builtins()
slice    <- builtins$slice

# ---- 5) Open dataset (NumPy chunk manager; chunks disabled) -----------------
open_ds <- function() {
  xr$open_dataset(
    EDH_URL,
    engine = "zarr",
    # Disable chunking -> avoid Dask. xarray will load NumPy arrays when accessed.  [docs]
    chunks = NULL,
    storage_options = dict(client_kwargs = dict(trust_env = TRUE)),
    # Tell the backend to use NumPy chunk manager (not Dask)
    backend_kwargs  = dict(consolidated = TRUE, chunked_array_type = "numpy")
  )
}

ds <- tryCatch(open_ds(), error = function(e) {
  stop("Failed to open EDH dataset. Check your PAT/network.\nDetails: ", e$message)
})

# Confirm variable present
vars <- names(py_to_r(ds$data_vars))
if (!"swvl4" %in% vars) {
  ds$close()
  stop("Variable 'swvl4' not found. Available (first 50): ",
       paste(utils::head(vars, 50), collapse = ", "),
       "\n(ERA5-Land variables list includes 'swvl4'; confirm endpoint and PAT.)")
}

# ---- 6) Monthly loop (subset bbox+time; write NetCDF) -----------------------
min_bytes  <- 10L * 1024L
retry_wait <- function(k) { Sys.sleep(min(60, 3 * 2^k)) }

for (y in years) {
  message(sprintf("📅 Year %d ...", y))
  for (m in months) {
    first <- make_date(y, m, 1)
    last  <- ceiling_date(first, unit = "month") - days(1)
    
    out_nc <- file.path(out_dir, sprintf("era5land_swvl4_%d_%02d.nc", y, m))
    if (file.exists(out_nc) && is.finite(file.size(out_nc)) && file.size(out_nc) >= min_bytes) {
      message(sprintf("⏭️  %d-%02d already present (%.2f MB)", y, m, file.size(out_nc)/(1024^2)))
      next
    }
    
    attempt <- 0; ok <- FALSE
    while (attempt < 3 && !ok) {
      attempt <- attempt + 1
      message(sprintf("⬇️  %d-%02d (attempt %d/3)", y, m, attempt))
      res <- try({
        swvl4 <- ds$get_item("swvl4")$
          sel(time      = slice(as.character(first), as.character(last)))$
          sel(latitude  = slice(latN, latS),
              longitude = slice(lonW, lonE))
        swvl4$to_netcdf(out_nc)
        swvl4$close()
        ok <- file.exists(out_nc) && file.size(out_nc) >= min_bytes
      }, silent = TRUE)
      
      if (!ok && attempt < 3) {
        message("⚠️  transient issue; backing off ...")
        retry_wait(attempt)
      }
    }
    
    if (ok) {
      message(sprintf("💾 Wrote %s (%.2f MB)", basename(out_nc), file.size(out_nc)/(1024^2)))
    } else {
      warning(sprintf("❌ Failed to save %d-%02d after retries.", y, m))
    }
  }
}

ds$close()
message("🎉 Done. Files saved in: ", out_dir)
