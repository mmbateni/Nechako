####################################################################################
# HYDAT Database — Data Download & GRDC Coverage Validation
# Mirrors data_downloads_hy_validate_embedSt.R but uses tidyhydat (HYDAT) 
# instead of the WSC GeoMet API to fetch deep historical data.
####################################################################################

library(tidyhydat)
library(tidyverse)

# ── 0  Paths ────────────────────────────────────────────────────────────────────
setwd("D:/Nechako_Drought/Nechako/Hydrology/data_retrievalWaterSurveyofCanada")
output_dir <- "HYDAT_downloads_ver2"
report_dir <- "HYDAT_validation_ver2"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(report_dir)) dir.create(report_dir, recursive = TRUE)

message("Output  : ", output_dir)
message("Reports : ", report_dir)
message(strrep("-", 65))

# ── 1  Target stations & GRDC reference spans ──────────────────────────────────
# Exact same 7 stations and GRDC mappings as the GeoMet script
stations_df <- tribble(
  ~STATION_NUMBER, ~STATION_NAME,                      ~GRDC_NO, ~GRDC_RIVER,
  "08JA015",       "LAVENTIE CREEK NEAR THE MOUTH",    4207050,  "LAVENTIE CREEK",
  "08JB002",       "STELLAKO RIVER AT GLENANNAN",      4207100,  "STELLAKO RIVER",
  "08JC001",       "NECHAKO RIVER AT VANDERHOOF",      4207160,  "NECHAKO RIVER",
  "08JC002",       "NECHAKO RIVER AT ISLE PIERRE",     4207180,  "NECHAKO RIVER",
  "08KC001",       "SALMON RIVER NEAR PRINCE GEORGE",  4207220,  "SALMON RIVER",
  "08JE001",       "STUART RIVER NEAR FORT ST. JAMES", 4207150,  "STUART RIVER",
  "08JE004",       "TSILCOH RIVER NEAR THE MOUTH",     4207130,  "TSILCOH RIVER"
)

grdc_spans <- tribble(
  ~GRDC_NO, ~grdc_d_start,         ~grdc_d_end,           ~grdc_d_yrs,
  4207050,  as.Date("1976-01-01"), as.Date("2025-12-31"), 47,
  4207100,  as.Date("1929-01-01"), as.Date("2025-12-31"), 96,
  4207160,  as.Date("1915-01-01"), as.Date("2025-12-31"), 109,
  4207180,  as.Date("1950-01-01"), as.Date("2025-12-31"), 74,
  4207220,  as.Date("1953-01-01"), as.Date("2025-12-31"), 70,
  4207150,  as.Date("1929-01-01"), as.Date("2025-12-31"), 96,
  4207130,  as.Date("1975-01-01"), as.Date("2025-12-31"), 50
)

stations_df <- stations_df %>% dplyr::left_join(grdc_spans, by = "GRDC_NO")

# ── 2  One-time HYDAT database download ───────────────────────────────────────
if (!file.exists(hy_default_db())) {
  message("Downloading HYDAT database (one-time, ~1.5 GB) ...")
  download_hydat()
  message("HYDAT download complete.")
} else {
  message("HYDAT database already present.")
}

# ── 3  Fetch and Save Data via tidyhydat ──────────────────────────────────────
station_ids <- stations_df$STATION_NUMBER

# Helper to save data per station
save_hydat_data <- function(data, param_name, out_dir) {
  if (is.null(data) || nrow(data) == 0) return(invisible(NULL))
  
  col_name <- if (param_name == "DISCHARGE") "Discharge_cms" else "Level_m"
  
  data %>%
    dplyr::group_by(STATION_NUMBER) %>%
    dplyr::group_walk(function(df, grp) {
      stn <- grp$STATION_NUMBER
      df <- df %>% 
        dplyr::arrange(Date) %>% 
        dplyr::select(Date, Value) %>%  # Explicitly select needed columns
        dplyr::rename(!!col_name := Value)
      
      out_file <- file.path(out_dir, sprintf("%s_HYDAT_%s.csv", stn, param_name))
      readr::write_csv(df, out_file)
      message(sprintf("  -> %s: %d rows saved (%s to %s)", 
                      stn, nrow(df), min(df$Date), max(df$Date)))
    })
}
message("\nFetching daily DISCHARGE from HYDAT...")
flow_hist <- tryCatch(
  hy_daily_flows(station_number = station_ids),
  error = function(e) { message("ERROR fetching flow: ", e$message); NULL }
)
save_hydat_data(flow_hist, "DISCHARGE", output_dir)

message("\nFetching daily LEVEL from HYDAT...")
level_hist <- tryCatch(
  hy_daily_levels(station_number = station_ids),
  error = function(e) { message("ERROR fetching level: ", e$message); NULL }
)
save_hydat_data(level_hist, "LEVEL", output_dir)

# ── 4  Validation against GRDC time spans ─────────────────────────────────────
message(strrep("=", 65))
message("VALIDATION: HYDAT discharge span vs GRDC reference span")
message(strrep("=", 65))

# Calculate spans from downloaded flow data
flow_summary <- if (!is.null(flow_hist) && nrow(flow_hist) > 0) {
  flow_hist %>%
    dplyr::group_by(STATION_NUMBER) %>%
    dplyr::summarise(
      wsc_start = min(Date, na.rm = TRUE),
      wsc_end   = max(Date, na.rm = TRUE),
      wsc_rows  = dplyr::n(),
      .groups = "drop"
    )
} else {
  tibble(STATION_NUMBER = character(), wsc_start = as.Date(character()), 
         wsc_end = as.Date(character()), wsc_rows = integer())
}

validation <- stations_df %>%
  dplyr::left_join(flow_summary, by = "STATION_NUMBER") %>%
  dplyr::mutate(
    wsc_yrs       = as.integer(format(wsc_end, "%Y")) - as.integer(format(wsc_start, "%Y")) + 1,
    wsc_start_yr  = as.integer(format(wsc_start, "%Y")),
    wsc_end_yr    = as.integer(format(wsc_end, "%Y")),
    grdc_start_yr = as.integer(format(grdc_d_start, "%Y")),
    grdc_end_yr   = as.integer(format(grdc_d_end, "%Y")),
    covers_start  = !is.na(wsc_start) & wsc_start_yr <= grdc_start_yr,
    covers_end    = !is.na(wsc_end)   & wsc_end_yr   >= grdc_end_yr,
    passes        = covers_start & covers_end,
    start_gap_days = as.integer(wsc_start - grdc_d_start),
    end_gap_days   = as.integer(wsc_end - grdc_d_end)
  )

# Console summary
for (i in seq_len(nrow(validation))) {
  v <- validation[i, ]
  status <- if (isTRUE(v$passes)) "✅ PASS" else "❌ FAIL"
  
  message(sprintf("\n%s  %s — %s", status, v$STATION_NUMBER, v$STATION_NAME))
  message(sprintf("  GRDC ref : %s → %s  (%d yrs)", v$grdc_d_start, v$grdc_d_end, v$grdc_d_yrs))
  
  if (!is.na(v$wsc_start)) {
    message(sprintf("  HYDAT data: %s → %s  (%d yrs, %d rows)", 
                    v$wsc_start, v$wsc_end, v$wsc_yrs, v$wsc_rows))
    message(sprintf("  Start gap: %+d days  |  End gap: %+d days", v$start_gap_days, v$end_gap_days))
    
    if (!v$covers_start) message(sprintf("  ⚠  HYDAT starts %d days AFTER GRDC d_start", v$start_gap_days))
    if (!v$covers_end)   message(sprintf("  ⚠  HYDAT ends   %d days BEFORE GRDC d_end", -v$end_gap_days))
  } else {
    message("  ⚠  No HYDAT discharge data retrieved — cannot validate")
  }
}

# ── 5  Save validation report ─────────────────────────────────────────────────
report <- validation %>%
  dplyr::select(
    STATION_NUMBER, STATION_NAME, GRDC_NO, GRDC_RIVER,
    grdc_d_start, grdc_d_end, grdc_d_yrs,
    wsc_start, wsc_end, wsc_yrs, wsc_rows,
    covers_start, covers_end, passes,
    start_gap_days, end_gap_days
  ) %>%
  dplyr::mutate(dplyr::across(c(grdc_d_start, grdc_d_end, wsc_start, wsc_end), as.character))

report_file <- file.path(report_dir, "GRDC_HYDAT_coverage_validation.csv")
readr::write_csv(report, report_file)

message(strrep("-", 65))
message("Validation report saved: ", report_file)
n_pass <- sum(validation$passes, na.rm = TRUE)
n_fail <- nrow(validation) - n_pass
message(sprintf("Summary: %d PASS  |  %d FAIL  (out of %d stations)", n_pass, n_fail, nrow(validation)))
message(strrep("=", 65))
message("✅  HYDAT Download & validation complete.")