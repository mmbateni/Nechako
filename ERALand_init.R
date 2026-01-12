# ---- 1. Setup & Libraries ----
library(ecmwfr)    # Pure R interface to CDS (No Python setup needed)
library(terra)     # Memory-safe raster handling
library(lubridate) # Robust date handling (Handles leap years correctly)

# ---- 1. Credentials ----
# Prefer environment variables:
user_uid <- "mehdi.bateni@iusspavia.it"
user_key <- "f0b6de88-c627-4f98-98b4-22453fb77cac"


if (is.na(user_uid) || is.na(user_key) || !nzchar(user_uid) || !nzchar(user_key)) {
  stop("CDS credentials not found. Set CDS_USER and CDS_KEY env vars, or hardcode temporarily.")
}

# Register the key
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
keep_archives <- FALSE   # Set TRUE to keep .zip files after extraction

# ---- 3. Helper utilities ----

# Canonical monthly base name (no extension)
monthly_basename <- function(year, month) {
  sprintf("era5_%s_%d_%02d", variable_name, year, month)
}
monthly_nc_path <- function(year, month) file.path(output_directory, paste0(monthly_basename(year, month), ".nc"))
monthly_zip_path <- function(year, month) file.path(output_directory, paste0(monthly_basename(year, month), ".zip"))
monthly_any_pattern <- function(year, month) sprintf("^%s(\\..*)?$", monthly_basename(year, month))

# Detect ZIP by magic bytes (works even if extension is missing or misleading)
is_zip_file <- function(fp) {
  if (!file.exists(fp)) return(FALSE)
  con <- file(fp, "rb")
  on.exit(close(con), add = TRUE)
  header <- tryCatch(readBin(con, "raw", n = 4), error = function(e) raw(0))
  length(header) >= 4 && rawToChar(header[1:2]) == "PK"
}

# Extract a ZIP into .nc and rename to our canonical name
ensure_nc_from_archive <- function(zip_fp, desired_nc_fp) {
  if (!file.exists(zip_fp)) return(FALSE)
  # List the archive contents
  listing <- tryCatch(utils::unzip(zip_fp, list = TRUE), error = function(e) NULL)
  if (is.null(listing) || nrow(listing) == 0) return(FALSE)
  
  # Find a .nc inside the archive (common CDS behavior)
  inside_nc <- listing$Name[grepl("\\.nc$", listing$Name, ignore.case = TRUE)]
  if (length(inside_nc) == 0) {
    # If no explicit .nc, extract all and look in output dir
    utils::unzip(zip_fp, exdir = output_directory, overwrite = TRUE)
    candidates <- list.files(output_directory, pattern = "\\.nc$", full.names = TRUE)
    if (length(candidates) == 0) return(FALSE)
    # Use the most recent extracted .nc
    info <- file.info(candidates)
    inside_fp <- candidates[which.max(info$mtime)]
  } else {
    # Extract just the .nc we found
    utils::unzip(zip_fp, files = inside_nc, exdir = output_directory, overwrite = TRUE)
    inside_fp <- file.path(output_directory, inside_nc[1])
  }
  
  # Rename to the canonical path (if different)
  if (!file.exists(inside_fp)) return(FALSE)
  if (normalizePath(inside_fp, winslash = "/", mustWork = FALSE) !=
      normalizePath(desired_nc_fp, winslash = "/", mustWork = FALSE)) {
    ok <- tryCatch(file.rename(inside_fp, desired_nc_fp), error = function(e) FALSE)
    if (!ok) {
      # Fall back to copying if rename fails across partitions
      ok <- tryCatch(file.copy(inside_fp, desired_nc_fp, overwrite = TRUE), error = function(e) FALSE)
      if (!ok) return(FALSE)
      unlink(inside_fp)
    }
  }
  
  # Optionally remove the archive to avoid confusion
  if (!keep_archives) try(unlink(zip_fp), silent = TRUE)
  
  file.exists(desired_nc_fp) && file.size(desired_nc_fp) >= presence_min_bytes
}

# Ensure local .nc exists (extract from ZIP or “no-extension ZIP” if needed)
ensure_local_nc <- function(year, month) {
  nc_fp  <- monthly_nc_path(year, month)
  if (file.exists(nc_fp) && file.size(nc_fp) >= presence_min_bytes) return(TRUE)
  
  # Look for artifacts: .zip, .nc, or file with no extension but matching base
  pattern <- monthly_any_pattern(year, month)
  artifacts <- list.files(output_directory, pattern = pattern, full.names = TRUE)
  if (length(artifacts) == 0) return(FALSE)
  
  # If any artifact is a ZIP (even if the name ends with .nc or has no extension), extract it
  for (fp in artifacts) {
    if (is_zip_file(fp)) {
      if (ensure_nc_from_archive(fp, nc_fp)) return(TRUE)
    }
  }
  
  # If we already have .nc but too small, treat as not present
  file.exists(nc_fp) && file.size(nc_fp) >= presence_min_bytes
}

# Check if a month is available locally (ensures .nc presence before we decide)
has_month_locally <- function(year, month) {
  ensure_local_nc(year, month)
}

# Build the CDS request payload for a single year-month
build_request <- function(year, month, days_list, target_basename_zip) {
  list(
    product_type = "reanalysis",
    format       = "netcdf",                  # CDS typically returns a zip containing the .nc
    variable     = variable_name,
    year         = as.character(year),
    month        = sprintf("%02d", month),
    day          = days_list,
    time         = times,
    area         = area_bbox,
    dataset_short_name = dataset,
    target       = target_basename_zip        # We set the target explicitly to a ZIP filename
  )
}

# Months with small final .nc (< 1 MB) for reporting
small_months_for_year <- function(year, threshold_bytes = report_small_bytes) {
  small <- integer(0)
  for (m in months_of_interest) {
    nc_fp <- monthly_nc_path(year, m)
    # Try to ensure the .nc exists before checking size
    ensure_local_nc(year, m)
    if (file.exists(nc_fp)) {
      sz <- file.size(nc_fp)
      if (!is.na(sz) && sz < threshold_bytes) small <- c(small, m)
    }
  }
  small
}

# ---- 4. Download a month ONLY IF NEEDED (check & unzip BEFORE building request) ----
download_if_needed <- function(year, month) {
  base_name   <- monthly_basename(year, month)
  file_name_nc  <- paste0(base_name, ".nc")
  file_name_zip <- paste0(base_name, ".zip")
  full_path_nc  <- monthly_nc_path(year, month)
  full_path_zip <- monthly_zip_path(year, month)
  
  # ---- EXISTENCE/EXTRACTION CHECK BEFORE BUILDING THE CDS REQUEST ----
  if (has_month_locally(year, month)) {
    message(sprintf("✅ Local available: %s (%.2f MB) — no request sent.",
                    file_name_nc, file.size(full_path_nc) / (1024^2)))
    return(invisible(TRUE))
  }
  
  # Build days list only after confirming local absence
  current_date <- lubridate::make_date(year, month, 1)
  days_in_current_month <- lubridate::days_in_month(current_date)
  days_list <- sprintf("%02d", 1:days_in_current_month)
  
  # ---- Build the CDS request AFTER local absence confirmed ----
  req <- build_request(year, month, days_list, target_basename_zip = file_name_zip)
  
  # Retry logic
  max_attempts <- 3L
  for (attempt in 1:max_attempts) {
    message(sprintf("⬇️  Downloading %d-%02d (attempt %d/%d) ...", year, month, attempt, max_attempts))
    ok <- FALSE
    try({
      ecmwfr::wf_request(
        request = req,
        user    = user_uid,
        path    = output_directory,
        verbose = TRUE
      )
      ok <- TRUE
    }, silent = TRUE)
    
    # After download, normalize: if ZIP exists, extract to .nc; if not, try detecting ZIP magic on target
    # Handle cases where CDS saved a ZIP without the .zip extension:
    extracted <- FALSE
    if (file.exists(full_path_zip) && is_zip_file(full_path_zip)) {
      extracted <- ensure_nc_from_archive(full_path_zip, full_path_nc)
    } else if (file.exists(full_path_nc) && is_zip_file(full_path_nc)) {
      # The "nc" file is actually a ZIP; extract and rename
      extracted <- ensure_nc_from_archive(full_path_nc, full_path_nc)
    } else {
      # Look for any artifact matching the monthly base, and extract if ZIP
      pattern <- monthly_any_pattern(year, month)
      artifacts <- list.files(output_directory, pattern = pattern, full.names = TRUE)
      for (fp in artifacts) {
        if (is_zip_file(fp)) {
          if (ensure_nc_from_archive(fp, full_path_nc)) { extracted <- TRUE; break }
        }
      }
    }
    
    # Final verification: do we have a usable .nc of minimal size?
    if ((ok || extracted) && file.exists(full_path_nc) && file.size(full_path_nc) >= presence_min_bytes) {
      message(sprintf("💾 Ready: %s (%.2f MB)", file_name_nc, file.size(full_path_nc) / (1024^2)))
      return(invisible(TRUE))
    }
    
    if (attempt < max_attempts) {
      wait <- 5 * attempt
      message(sprintf("⚠️  Verification failed (no usable .nc). Retrying in %ds ...", wait))
      Sys.sleep(wait)
    }
  }
  
  warning(sprintf("❌ FAILED to prepare %d-%02d after %d attempts.", year, month, max_attempts))
  return(invisible(FALSE))
}

# ---- 5. Year loop: status + download-if-needed ----
for (year in start_year:end_year) {
  
  # Determine status based on finalized .nc presence
  present <- sapply(months_of_interest, function(m) has_month_locally(year, m))
  missing <- months_of_interest[!present]
  small   <- small_months_for_year(year, threshold_bytes = report_small_bytes)
  
  if (length(missing) == 0) {
    if (length(small) == 0) {
      message(sprintf("⏭️  Skipping %d — all 12 months present.", year))
    } else {
      message(sprintf("⏭️  Skipping %d — all 12 months present. ⚠️ Small-size months (<1MB): %s",
                      year, paste(sprintf("%02d", sort(small)), collapse = ", ")))
    }
    next
  } else {
    present_count <- 12 - length(missing)
    msg <- sprintf("📅 %d — %d/%d months present; missing: %s",
                   year, present_count, 12, paste(sprintf("%02d", sort(missing)), collapse = ", "))
    if (length(small) > 0) {
      msg <- paste0(msg, sprintf(" | ⚠️ Small-size months (<1MB): %s",
                                 paste(sprintf("%02d", sort(small)), collapse = ", ")))
    }
    message(msg)
  }
  
  # Iterate all months; for each, check/extract first, then only request if needed
  for (m in months_of_interest) {
    download_if_needed(year, m)
  }
}

message("✅ Downloads finished. Proceeding to merge ...")

# ---- 6. Smart Merge (CDO if available; Terra fallback) ----

# Collect only .nc files (ignore any leftover archives)
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
    message("⚠️  CDO merge did not produce a valid file. Will fall back to Terra.")
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

  