# ---- 1. Setup & Libraries ----
library(ecmwfr)    # Pure R interface to CDS (No Python setup needed)
library(terra)     # Memory-safe raster handling
library(lubridate) # Robust date handling (Handles leap years correctly)

# ---- 1. Credentials ----
# Prefer environment variables:
user_uid <- "bateni2000@gmail.com"
user_key <- "47081be6-c313-4832-9275-c2c5deeed751"

if (is.na(user_uid) || is.na(user_key) || !nzchar(user_uid) || !nzchar(user_key)) {
  stop("CDS credentials not found. Set CDS_USER and CDS_KEY env vars, or hardcode temporarily.")
}

# Register the key (PAT) with ecmwfr
ecmwfr::wf_set_key(user = user_uid, key = user_key)

# ---- 2. Settings ----
dataset        <- "reanalysis-era5-land"
variable_name  <- "snowfall"

base_download_directory <- "D:/Nechako_Drought/ERA5-Land download-N and B"
output_directory        <- file.path(base_download_directory, variable_name)
if (!dir.exists(output_directory)) dir.create(output_directory, recursive = TRUE)

start_year <- 1950
end_year   <- 2024
months_of_interest <- 1:12
times <- sprintf("%02d:00", 0:23)             # "00:00" to "23:00"
area_bbox <- c(56.5, -128.5, 51.9, -122.0)    # [North, West, South, East]

# Thresholds
presence_min_bytes <- 10L * 1024L          # Treat file as "present" if >= 10 KB
report_small_bytes <- 1L  * 1024L * 1024L  # Report months smaller than 1 MB

# Post-download handling
keep_archives <- FALSE   # TRUE to keep ZIPs after extraction

# Parallel batch tuning (CDS typically allows up to ~20 parallel jobs)
max_workers <- 12L       # adjust to your preference ≤ 20  (CDS queue limit)  [1] [1](https://bluegreen-labs.github.io/ecmwfr/articles/advanced_vignette.html)[2](https://cloud.r-project.org/web/packages/ecmwfr/vignettes/advanced_vignette.html)
batch_retry_secs <- 60L  # polling interval for batch jobs (default is 5–30s)
batch_timeout    <- 3600 # how long to wait for download to start (per job)

# ---- 3. Helper utilities (ZIP-aware existence & extraction) ----

# Canonical monthly base name (no extension)
monthly_basename <- function(year, month) {
  sprintf("era5_%s_%d_%02d", variable_name, year, month)
}
monthly_nc_path <- function(year, month) file.path(output_directory, paste0(monthly_basename(year, month), ".nc"))
monthly_zip_path <- function(year, month) file.path(output_directory, paste0(monthly_basename(year, month), ".zip"))
monthly_any_pattern <- function(year, month) sprintf("^%s(\\..*)?$", monthly_basename(year, month))

# Detect ZIP by magic bytes (works even if extension is absent/misleading)
is_zip_file <- function(fp) {
  if (!file.exists(fp)) return(FALSE)
  con <- file(fp, "rb")
  on.exit(close(con), add = TRUE)
  header <- tryCatch(readBin(con, "raw", n = 4), error = function(e) raw(0))
  length(header) >= 2 && rawToChar(header[1:2]) == "PK"
}

# Extract a ZIP into .nc and rename to canonical name
ensure_nc_from_archive <- function(zip_fp, desired_nc_fp) {
  if (!file.exists(zip_fp)) return(FALSE)
  listing <- tryCatch(utils::unzip(zip_fp, list = TRUE), error = function(e) NULL)
  if (is.null(listing) || nrow(listing) == 0) return(FALSE)
  
  inside_nc <- listing$Name[grepl("\\.nc$", listing$Name, ignore.case = TRUE)]
  if (length(inside_nc) == 0) {
    utils::unzip(zip_fp, exdir = output_directory, overwrite = TRUE)
    candidates <- list.files(output_directory, pattern = "\\.nc$", full.names = TRUE)
    if (length(candidates) == 0) return(FALSE)
    info <- file.info(candidates)
    inside_fp <- candidates[which.max(info$mtime)]
  } else {
    utils::unzip(zip_fp, files = inside_nc, exdir = output_directory, overwrite = TRUE)
    inside_fp <- file.path(output_directory, inside_nc[1])
  }
  
  if (!file.exists(inside_fp)) return(FALSE)
  if (normalizePath(inside_fp, winslash = "/", mustWork = FALSE) !=
      normalizePath(desired_nc_fp, winslash = "/", mustWork = FALSE)) {
    ok <- tryCatch(file.rename(inside_fp, desired_nc_fp), error = function(e) FALSE)
    if (!ok) {
      ok <- tryCatch(file.copy(inside_fp, desired_nc_fp, overwrite = TRUE), error = function(e) FALSE)
      if (!ok) return(FALSE)
      unlink(inside_fp)
    }
  }
  
  if (!keep_archives) try(unlink(zip_fp), silent = TRUE)
  file.exists(desired_nc_fp) && file.size(desired_nc_fp) >= presence_min_bytes
}

# Ensure local .nc exists (extract from ZIP or “no-extension ZIP” if needed)
ensure_local_nc <- function(year, month) {
  nc_fp  <- monthly_nc_path(year, month)
  if (file.exists(nc_fp) && file.size(nc_fp) >= presence_min_bytes) return(TRUE)
  
  pattern <- monthly_any_pattern(year, month)
  artifacts <- list.files(output_directory, pattern = pattern, full.names = TRUE)
  if (length(artifacts) == 0) return(FALSE)
  
  for (fp in artifacts) {
    if (is_zip_file(fp)) {
      if (ensure_nc_from_archive(fp, nc_fp)) return(TRUE)
    }
  }
  
  file.exists(nc_fp) && file.size(nc_fp) >= presence_min_bytes
}

# High-level check
has_month_locally <- function(year, month) {
  ensure_local_nc(year, month)
}

# Return basenames of existing monthly files for a given year (robust pattern)
existing_month_files_for_year <- function(year) {
  pattern <- sprintf("^era5_%s_%d_\\d{2}\\.nc$", variable_name, year)
  basename(list.files(output_directory, pattern = pattern, full.names = FALSE))
}

# Months that are missing (by finalized .nc presence)
missing_months_for_year <- function(year, months_vec = months_of_interest) {
  expected <- sapply(months_vec, \(m) sprintf("era5_%s_%d_%02d.nc", variable_name, year, m), USE.NAMES = FALSE)
  existing <- existing_month_files_for_year(year)
  missing_basenames <- setdiff(expected, existing)
  if (length(missing_basenames) == 0) return(integer(0))
  as.integer(sub("^.*_(\\d{4})_(\\d{2})\\.nc$", "\\2", missing_basenames))
}

# Months with small final .nc (< 1 MB) for reporting
small_months_for_year <- function(year, threshold_bytes = report_small_bytes) {
  small <- integer(0)
  for (m in months_of_interest) {
    nc_fp <- monthly_nc_path(year, m)
    ensure_local_nc(year, m) # normalize/extract before checking
    if (file.exists(nc_fp)) {
      sz <- file.size(nc_fp)
      if (!is.na(sz) && sz < threshold_bytes) small <- c(small, m)
    }
  }
  small
}

# ---- 4. Build monthly requests ONLY for months not available locally ----

# Build a single monthly request (existence check occurs BEFORE building request)
build_month_request <- function(year, month) {
  # Existence check first (don’t build a CDS request if present locally)
  if (has_month_locally(year, month)) return(NULL)
  
  current_date <- lubridate::make_date(year, month, 1)
  days_list <- sprintf("%02d", 1:lubridate::days_in_month(current_date))
  target_name <- sprintf("era5_%s_%d_%02d.nc", variable_name, year, month)
  
  list(
    product_type       = "reanalysis",
    format             = "netcdf",
    variable           = variable_name,
    year               = as.character(year),
    month              = sprintf("%02d", month),
    day                = days_list,
    time               = times,
    area               = area_bbox,
    dataset_short_name = dataset,
    target             = target_name,
    # If supported by the dataset, this avoids ZIP packaging; safe to include
    download_format    = "unarchived"   # dataset-dependent; ignored if unsupported  [2] [5](https://training.pages.data.meditwin-project.eu/meditwin-summer-school-2025/chapter1/6_download_cds_api.html)
  )
}

# Collect all missing-month requests
requests <- list()
for (year in start_year:end_year) {
  # Status & reporting
  present <- sapply(months_of_interest, function(m) has_month_locally(year, m))
  missing <- months_of_interest[!present]
  small   <- small_months_for_year(year, threshold_bytes = report_small_bytes)
  
  if (length(missing) == 0) {
    if (length(small) == 0) {
      message(sprintf("⏭️  %d — all 12 months present.", year))
    } else {
      message(sprintf("⏭️  %d — all months present. ⚠️ Small-size (<1MB): %s",
                      year, paste(sprintf("%02d", sort(small)), collapse = ", ")))
    }
    next
  } else {
    present_count <- 12 - length(missing)
    msg <- sprintf("📅 %d — %d/%d months present; missing: %s",
                   year, present_count, 12, paste(sprintf("%02d", sort(missing)), collapse = ", "))
    if (length(small) > 0) {
      msg <- paste0(msg, sprintf(" | ⚠️ Small-size (<1MB): %s",
                                 paste(sprintf("%02d", sort(small)), collapse = ", ")))
    }
    message(msg)
  }
  
  # Build requests for missing months only
  for (m in missing) {
    req <- build_month_request(year, m)
    if (!is.null(req)) {
      requests[[length(requests) + 1]] <- req
    }
  }
}

message(sprintf("🧾 Prepared %d monthly requests needing download.", length(requests)))

# ---- 5. Submit requests in PARALLEL batch (up to ~20 workers) ----
if (length(requests) > 0) {
  workers <- min(max_workers, length(requests))
  message(sprintf("🚀 Submitting %d parallel requests (workers = %d) ...", length(requests), workers))
  # Per ecmwfr docs: workers ≤ 20; retry controls polling frequency; time_out controls start wait. [1] [3](https://bluegreen-labs.github.io/ecmwfr/reference/wf_request.html)[4](https://rdrr.io/cran/ecmwfr/man/wf_request.html)
  ecmwfr::wf_request_batch(
    request_list  = requests,
    workers       = workers,
    user          = user_uid,
    path          = output_directory,
    time_out      = batch_timeout,
    retry         = batch_retry_secs
    # total_timeout has a default; you can set if desired:
    # total_timeout = length(requests) * batch_timeout / workers
  )
} else {
  message("✅ Nothing to download — all requested months are already available locally.")
}

# After batch completes, normalize any leftover archives or no-extension files
message("🔎 Post-batch normalization (ZIP extraction & size checks)...")
for (year in start_year:end_year) {
  for (m in months_of_interest) {
    # Ensure final .nc file exists & is minimally sized
    if (!has_month_locally(year, m)) {
      # Try to extract from any artifacts if present
      ensure_local_nc(year, m)
    }
  }
}

message("✅ Downloads finished. Proceeding to merge ...")

# ---- 6. Smart Merge (CDO if available; Terra fallback) ----

# Collect only finalized .nc files
nc_files <- list.files(
  output_directory,
  pattern = sprintf("^era5_%s_\\d{4}_\\d{2}\\.nc$", variable_name),
  full.names = TRUE
)

if (length(nc_files) == 0) stop("No NetCDF files found to merge.")

# Sort lexicographically by year-month
nc_files <- nc_files[order(nc_files)]

merged_filename <- file.path(output_directory, sprintf("era5_land_%s_merged.nc", variable_name))
if (file.exists(merged_filename)) {
  message(sprintf("🧹 Removing old merged file: %s", merged_filename))
  unlink(merged_filename)
}

# Option A: CDO (fastest, if available)
cdo_path <- Sys.which("cdo")
merge_done <- FALSE

if (nzchar(cdo_path)) {
  message("🚀 CDO found — using high-speed merge (mergetime) ...")
  args <- c("mergetime", shQuote(nc_files), shQuote(merged_filename))
  res <- tryCatch(system2(cdo_path, args = args, stdout = TRUE, stderr = TRUE), error = function(e) e)
  
  if (file.exists(merged_filename) && isTRUE(file.size(merged_filename) > 1024)) {
    message(sprintf("✅ Merge successful via CDO: %s", merged_filename))
    merge_done <- TRUE
  } else {
    message("⚠️  CDO merge did not produce a valid file. Falling back to Terra.")
    if (inherits(res, "error")) message("CDO error: ", res$message) else message(paste(res, collapse = "\n"))
  }
} else {
  message("ℹ️  CDO not found on PATH; will use Terra fallback.")
}

# Option B: Terra fallback
if (!merge_done) {
  message("🧱 Merging with Terra (this can take a while for many files) ...")
  r_stack <- terra::rast(nc_files)
  tm <- tryCatch(terra::time(r_stack), error = function(e) NULL)
  if (!is.null(tm)) {
    ord <- order(tm)
    r_stack <- r_stack[[ord]]
  }
  
  terra::writeCDF(
    r_stack,
    merged_filename,
    overwrite   = TRUE,
    varname     = variable_name,
    compression = 5
  )
  
  if (file.exists(merged_filename) && isTRUE(file.size(merged_filename) > 1024)) {
    message(sprintf("✅ Merge successful via Terra: %s", merged_filename))
  } else {
    stop("❌ Terra merge failed: merged file missing or too small.")
  }
}

message("🎉 All done.")

  