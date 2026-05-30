# =============================================================================
# Batch RVT → CSV Converter
# Converts all Raven .rvt ObservationData files in a folder to individual CSVs
# and one combined CSV with all stations.
#
# Handles:
#   - Multiple comment lines (# prefix)
#   - :ObservationData header (variable type, basin ID, units)
#   - Timestep header (start datetime, timestep in days, n records)
#   - Missing values (-1.2345 → NA)
#   - Any timestep (daily, hourly, sub-daily)
# =============================================================================


# ── CONFIG: set these two paths ───────────────────────────────────────────────

rvt_folder  <- "./Spatial/rvn"   # folder containing .rvt files
output_folder <- "csv_folder"                # where CSVs are saved (same folder)
                                            # change to a different path if needed


# ── FUNCTION: parse one .rvt file ────────────────────────────────────────────

parse_rvt <- function(filepath) {

  lines <- readLines(filepath, warn = FALSE)
  lines <- trimws(lines)

  # --- find :ObservationData header line ---
  obs_idx <- grep("^:ObservationData", lines, ignore.case = TRUE)
  if (length(obs_idx) == 0) {
    warning("No :ObservationData found in: ", basename(filepath), " — skipping.")
    return(NULL)
  }
  obs_idx <- obs_idx[1]

  # parse: :ObservationData  TYPE  BASIN_ID  UNITS
  obs_parts  <- strsplit(trimws(lines[obs_idx]), "\\s+")[[1]]
  data_type  <- if (length(obs_parts) >= 2) obs_parts[2] else "UNKNOWN"
  basin_id   <- if (length(obs_parts) >= 3) obs_parts[3] else NA
  units      <- if (length(obs_parts) >= 4) obs_parts[4] else "unknown"

  # --- find timestep header (next non-blank, non-comment line after :ObservationData) ---
  ts_idx <- obs_idx + 1
  while (ts_idx <= length(lines) &&
         (lines[ts_idx] == "" || startsWith(lines[ts_idx], "#"))) {
    ts_idx <- ts_idx + 1
  }
  if (ts_idx > length(lines)) {
    warning("No timestep header found in: ", basename(filepath), " — skipping.")
    return(NULL)
  }

  # parse: START_DATE START_TIME TIMESTEP_DAYS N_RECORDS
  ts_parts   <- strsplit(lines[ts_idx], "\\s+")[[1]]
  ts_parts   <- ts_parts[ts_parts != ""]
  start_dt   <- as.POSIXct(paste(ts_parts[1], ts_parts[2]),
                            format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  timestep_days <- as.numeric(ts_parts[3])   # 1=daily, 0.04167≈1h, etc.
  n_records     <- as.integer(ts_parts[4])

  # convert timestep to seconds for seq()
  timestep_sec  <- timestep_days * 86400

  # --- generate datetime sequence ---
  datetimes <- seq(from = start_dt,
                   by   = timestep_sec,
                   length.out = n_records)

  # --- read data values (lines after ts_idx, skip :EndObservationData) ---
  data_start <- ts_idx + 1
  data_lines <- lines[data_start:length(lines)]
  data_lines <- data_lines[!grepl("^:|^#|^$", data_lines)]   # drop keywords/comments/blanks

  values_raw <- suppressWarnings(as.numeric(data_lines))
  values_raw <- values_raw[!is.na(values_raw)]                # drop non-numeric

  # guard against length mismatch
  n_use <- min(n_records, length(values_raw))
  if (n_use < n_records) {
    warning(basename(filepath), ": expected ", n_records,
            " records but found ", n_use, " — truncating.")
  }

  # replace Raven missing value with NA
  values <- values_raw[seq_len(n_use)]
  values[values == -1.2345] <- NA

  # --- assemble data frame ---
  # Label datetime column by resolution
  if (timestep_days == 1) {
    time_col  <- as.Date(datetimes[seq_len(n_use)])
    time_name <- "Date"
  } else {
    time_col  <- datetimes[seq_len(n_use)]
    time_name <- "Datetime_UTC"
  }

  station_id <- tools::file_path_sans_ext(basename(filepath))

  df <- data.frame(
    Time      = time_col,
    Value     = values,
    Station   = station_id,
    BasinID   = basin_id,
    DataType  = data_type,
    Units     = units,
    Timestep  = paste0(timestep_days * 24, "h"),
    stringsAsFactors = FALSE
  )
  names(df)[1] <- time_name

  return(df)
}


# ── BATCH CONVERSION ─────────────────────────────────────────────────────────

rvt_files <- list.files(rvt_folder, pattern = "\\.rvt$",
                         full.names = TRUE, ignore.case = TRUE)

if (length(rvt_files) == 0) {
  stop("No .rvt files found in: ", rvt_folder)
}

cat(sprintf("Found %d .rvt file(s) in: %s\n\n", length(rvt_files), rvt_folder))

all_data   <- list()   # accumulate for combined CSV
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)   # recursive = TRUE handles nested paths
  cat("Created output folder:", output_folder, "\n\n")
}
n_ok <- 0
n_skip     <- 0

for (f in rvt_files) {
  station <- tools::file_path_sans_ext(basename(f))
  cat(sprintf("  Processing: %-20s  ", basename(f)))

  df <- tryCatch(parse_rvt(f), error = function(e) {
    cat(sprintf("ERROR: %s\n", e$message))
    NULL
  })

  if (is.null(df)) {
    cat("SKIPPED\n")
    n_skip <- n_skip + 1
    next
  }

  # individual CSV per station
  out_csv <- file.path(output_folder,
                       paste0(station, "_flow.csv"))
  write.csv(df, out_csv, row.names = FALSE)

  n_valid  <- sum(!is.na(df$Value))
  n_miss   <- sum( is.na(df$Value))
  cat(sprintf("→ %d rows  |  %d valid  |  %d missing  →  %s\n",
              nrow(df), n_valid, n_miss, basename(out_csv)))

  all_data[[station]] <- df
  n_ok <- n_ok + 1
}


# ── COMBINED CSV (all stations, long format) ──────────────────────────────────

if (n_ok > 0) {
  combined      <- do.call(rbind, all_data)
  rownames(combined) <- NULL
  combined_path <- file.path(output_folder, "ALL_stations_flow.csv")
  write.csv(combined, combined_path, row.names = FALSE)

  cat(sprintf(
    "\n=== Done ===\n  Converted : %d file(s)\n  Skipped   : %d file(s)\n  Combined  : %s  (%d rows total)\n",
    n_ok, n_skip, basename(combined_path), nrow(combined)
  ))
} else {
  cat("\nNo files were successfully converted.\n")
}
