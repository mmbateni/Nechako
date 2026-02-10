# ---------- CONFIG ----------
library(ncdf4)

#setwd("D:/Nechako_Drought/monthly_data_direct/")
#nc_path <- "total_precipitation_monthly.nc"

setwd("D:/Nechako_Drought/spi_results/")
nc_path <- "spi_03month.nc"

# Define the range of time slices to export
# Option 1: Use index range (1 = first month, 900 = last month)
start_k <- 31L    # Starting time index (1 = January 1950)-841
end_k <- 35L     # Ending time index (12 = December 1950) - adjust as needed -900

# Option 2: Uncomment below to export ALL time slices (900 files!)
# end_k <- 900L

# Convert precipitation from meters -> millimeters (typical for ERA5 "tp")
convert_to_mm <- FALSE  # set TRUE if you want mm

# ---------- OPEN ----------
nc <- nc_open(nc_path)

# ---------- HELPERS ----------
nz <- function(x) if (is.null(x)) "" else as.character(x)

print_dims <- function(nc) {
  cat("===== DIMENSIONS =====\n")
  for (nm in names(nc$dim)) {
    d <- nc$dim[[nm]]
    cat(
      "\nName       :", nm,
      "\nLength     :", d$len,
      "\nUnits      :", nz(d$units),
      "\nUnlimited  :", isTRUE(d$unlim),
      "\n----------------------\n"
    )
  }
}

print_vars <- function(nc) {
  cat("\n===== VARIABLES =====\n")
  for (nm in names(nc$var)) {
    v <- nc$var[[nm]]
    dim_names <- paste(sapply(v$dim, `[[`, "name"), collapse = ", ")
    dim_lens  <- paste(sapply(v$dim, `[[`, "len"),  collapse = " x ")
    cat(
      "\nName        :", nm,
      "\nLong name   :", nz(v$longname),
      "\nUnits       :", nz(v$units),
      "\nDimensions  :", dim_names,
      "\nShape       :", dim_lens,
      "\n----------------------\n"
    )
  }
}

print_global_atts <- function(nc) {
  cat("\n===== GLOBAL ATTRIBUTES =====\n")
  if (length(nc$gatts) == 0L) {
    cat("(no global attributes)\n")
    return(invisible())
  }
  for (nm in names(nc$gatts)) {
    val <- nc$gatts[[nm]]$value
    cat(nm, ":", paste(as.character(val), collapse = " "), "\n")
  }
}

# ---------- PRINT SPECIFICATIONS ----------
print_dims(nc)
print_vars(nc)
print_global_atts(nc)

# ---------- DETECT TIME DIMENSION NAME ----------
# Find the time dimension (could be "valid_time", "time", etc.)
time_dim_name <- NULL
for (nm in names(nc$dim)) {
  d <- nc$dim[[nm]]
  # Check if units contain "since" (typical for time dimensions)
  if (!is.null(d$units) && grepl("since", d$units, ignore.case = TRUE)) {
    time_dim_name <- nm
    break
  }
}
# Fallback: if no "since" found, check common names
if (is.null(time_dim_name)) {
  common_names <- c("valid_time", "time", "Time", "TIME", "date")
  for (nm in common_names) {
    if (nm %in% names(nc$dim)) {
      time_dim_name <- nm
      break
    }
  }
}
# Last resort: use the longest dimension (usually time)
if (is.null(time_dim_name)) {
  dim_lens <- sapply(nc$dim, function(d) d$len)
  time_dim_name <- names(nc$dim)[which.max(dim_lens)]
}

cat(sprintf("\n===== DETECTED TIME DIMENSION: '%s' =====\n", time_dim_name))

# ---------- READ COORDS & TIME ----------
lon  <- ncvar_get(nc, "longitude")
lat  <- ncvar_get(nc, "latitude")
time <- ncvar_get(nc, time_dim_name)  # Use detected dimension name

# Decode time units for summary & filename
time_units <- nz(nc$dim[[time_dim_name]]$units)  # Use detected dimension name
origin <- as.POSIXct("1970-01-01 00:00:00", tz = "UTC")
mult <- 1
if (nz(time_units) != "" && grepl("since", time_units)) {
  parts  <- strsplit(time_units, " since ", fixed = TRUE)[[1]]
  unit_s <- tolower(trimws(parts[1]))
  origin <- as.POSIXct(trimws(parts[2]), tz = "UTC")
  mult <- switch(unit_s,
                 "seconds" = 1, "second" = 1,
                 "hours"   = 3600, "hour"  = 3600,
                 "days"    = 86400, "day"  = 86400,
                 1)
}
time_posix <- origin + time * mult

# Small summary so you can pick k
cat("\n===== TIME SUMMARY =====\n")
cat("Time length:", length(time), "\n")
if (length(time) > 0) {
  cat("Start:", format(min(time_posix), "%Y-%m-%d %H:%M:%S"), "|",
      "End:",   format(max(time_posix), "%Y-%m-%d %H:%M:%S"), "\n")
}
cat("------------------------\n")

# ---------- DETECT VARIABLE NAME ----------
# Find the main data variable (skip coordinate variables)
coord_vars <- c("longitude", "latitude", "lon", "lat", time_dim_name, "crs")
data_var_name <- NULL
for (nm in names(nc$var)) {
  if (!(nm %in% coord_vars)) {
    data_var_name <- nm
    break
  }
}

if (is.null(data_var_name)) {
  nc_close(nc)
  stop("Could not find a data variable in the file. Aborting.")
}

cat(sprintf("===== DETECTED DATA VARIABLE: '%s' =====\n", data_var_name))

# ---------- READ VARIABLE ----------
# Array is expected as [lon, lat, time]
tp <- ncvar_get(nc, data_var_name)

# Handle missing/fill values
fill <- nc$var[[data_var_name]]$missval
if (!is.null(fill) && !is.na(fill)) {
  tp[tp == fill] <- NA_real_
}

# Safety check for dimensions
if (length(dim(tp)) < 3) {
  nc_close(nc)
  stop(sprintf("Expected '%s' to have 3 dimensions [lon, lat, time]. Aborting.", data_var_name))
}
nt <- dim(tp)[3]

# Validate the time range
if (start_k < 1L || start_k > nt) {
  nc_close(nc)
  stop(sprintf("start_k=%d is out of range. Valid range is 1..%d", start_k, nt))
}
if (end_k < 1L || end_k > nt) {
  nc_close(nc)
  stop(sprintf("end_k=%d is out of range. Valid range is 1..%d", end_k, nt))
}
if (start_k > end_k) {
  nc_close(nc)
  stop(sprintf("start_k (%d) cannot be greater than end_k (%d)", start_k, end_k))
}

# Set row/column labels templates (same for all time slices)
lat_labels <- sprintf("%.3f", lat)
lon_labels <- sprintf("%.3f", lon)

# Calculate number of files to export
num_files <- end_k - start_k + 1
cat(sprintf("\n===== EXPORT SUMMARY =====\n"))
cat(sprintf("Exporting %d time slices from index %d to %d\n", num_files, start_k, end_k))
cat(sprintf("Date range: %s to %s\n", 
            format(time_posix[start_k], "%Y-%m-%d"),
            format(time_posix[end_k], "%Y-%m-%d")))
cat(sprintf("Output directory: %s\n", getwd()))
cat("------------------------\n")

# ---------- LOOP THROUGH TIME SLICES AND EXPORT ----------
cat("\nStarting export...\n")

for (k in start_k:end_k) {
  # Progress report
  progress <- (k - start_k + 1) / num_files * 100
  cat(sprintf("[%3.0f%%] Processing time slice %d of %d (%s)...", 
              progress, k, nt, format(time_posix[k], "%Y-%m")))
  
  # Slice the 3D array and transpose so rows = latitude, columns = longitude
  mat <- t(tp[, , k])
  
  if (isTRUE(convert_to_mm)) {
    mat <- mat * 1000  # meters -> millimeters
  }
  
  # Set row/column labels for readability in CSV
  rownames(mat) <- lat_labels
  colnames(mat) <- lon_labels
  
  # Name file using YYYYMM from the selected time
  stamp <- format(time_posix[k], "%Y%m")
  out_file <- sprintf("%s_matrix_%s.csv", data_var_name, stamp)  # Use detected variable name
  
  # Write the matrix CSV
  write.table(mat, file = out_file, sep = ",", row.names = TRUE,
              col.names = NA, quote = FALSE)
  
  cat(" Done\n")
}

# ---------- CLOSE ----------
nc_close(nc)

cat("\n===== EXPORT COMPLETE =====\n")
cat(sprintf("Successfully exported %d files from index %d to %d\n", 
            num_files, start_k, end_k))
cat(sprintf("Files were saved with format: %s_matrix_YYYYMM.csv\n", data_var_name))
cat(sprintf("Example file: %s_matrix_%s.csv\n", data_var_name, format(time_posix[start_k], "%Y%m")))
cat("====================================\n")