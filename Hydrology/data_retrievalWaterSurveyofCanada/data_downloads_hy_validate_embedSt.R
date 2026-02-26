####################################################################################
# WSC GeoMet API — Data Download & GRDC Coverage Validation
#
# Stations (WSC ↔ GRDC matched):
#   08JA015  Laventie Creek Near the Mouth       ↔  GRDC 4207050
#   08JB002  Stellako River at Glenannan          ↔  GRDC 4207100
#   08JC001  Nechako River at Vanderhoof          ↔  GRDC 4207160
#   08JC002  Nechako River at Isle Pierre         ↔  GRDC 4207180
#   08JC005  Chilako River Near Prince George     ↔  GRDC 4207220
#   08JE001  Stuart River Near Fort St. James     ↔  GRDC 4207150
#   08JE004  Tsilcoh River Near the Mouth         ↔  GRDC 4207130
#
# Validation rule: WSC daily discharge record must cover AT LEAST the full
# time span recorded in the GRDC-Monthly.nc file (d_start / d_end from
# GRDC_Stations catalogue). A station PASSES if:
#   WSC start date  <=  GRDC d_start  AND  WSC end date  >=  GRDC d_end
####################################################################################

library(tidyverse)
library(httr)
library(jsonlite)
library(purrr)

# ── 0  Paths ────────────────────────────────────────────────────────────────────
setwd("D:/Nechako_Drought/Nechako/Hydrology/data_retrievalWaterSurveyofCanada")
output_dir  <- "WSC_downloads_ver2"
report_dir  <- "WSC_validation_ver2"
if (!dir.exists(output_dir))  dir.create(output_dir,  recursive = TRUE)
if (!dir.exists(report_dir))  dir.create(report_dir,  recursive = TRUE)

message("Output  : ", output_dir)
message("Reports : ", report_dir)
message(strrep("-", 65))

# ── 1  Target stations ──────────────────────────────────────────────────────────
stations_df <- tribble(
  ~STATION_NUMBER, ~STATION_NAME,                           ~GRDC_NO, ~GRDC_RIVER,
  "08JA015", "LAVENTIE CREEK NEAR THE MOUTH",        4207050, "LAVENTIE CREEK",
  "08JB002", "STELLAKO RIVER AT GLENANNAN",           4207100, "STELLAKO RIVER",
  "08JC001", "NECHAKO RIVER AT VANDERHOOF",           4207160, "NECHAKO RIVER",
  "08JC002", "NECHAKO RIVER AT ISLE PIERRE",          4207180, "NECHAKO RIVER",
  "08JC005", "CHILAKO RIVER NEAR PRINCE GEORGE",      4207220, "SALMON RIVER",
  "08JE001", "STUART RIVER NEAR FORT ST. JAMES",      4207150, "STUART RIVER",
  "08JE004", "TSILCOH RIVER NEAR THE MOUTH",          4207130, "TSILCOH RIVER"
)

# ── 2  GRDC reference time spans (from GRDC_Stations.xlsx d_start / d_end)
#       These are the MINIMUM spans the WSC record must cover ──────────────────
grdc_spans <- tribble(
  ~GRDC_NO, ~grdc_d_start,           ~grdc_d_end,             ~grdc_d_yrs,
  4207050,  as.Date("1976-01-01"),   as.Date("2022-12-31"),   47,
  4207100,  as.Date("1929-01-01"),   as.Date("2024-12-31"),   96,
  4207160,  as.Date("1915-01-01"),   as.Date("2023-12-31"),   109,
  4207180,  as.Date("1950-01-01"),   as.Date("2023-12-31"),   74,
  4207220,  as.Date("1953-01-01"),   as.Date("2022-12-31"),   70,
  4207150,  as.Date("1929-01-01"),   as.Date("2024-12-31"),   96,
  4207130,  as.Date("1975-01-01"),   as.Date("2024-12-31"),   50
)

stations_df <- stations_df %>% left_join(grdc_spans, by = "GRDC_NO")

# ── 3  API helpers ──────────────────────────────────────────────────────────────
BASE_URL <- "https://api.weather.gc.ca/collections/hydrometric-daily-mean/items"

safe_get_json <- function(url, timeout_seconds = 120) {
  resp <- tryCatch(GET(url, timeout(timeout_seconds)), error = function(e) e)
  if (inherits(resp, "error")) {
    warning("HTTP error: ", resp$message); return(NULL)
  }
  if (http_error(resp)) {
    warning("HTTP ", status_code(resp), " — ", url); return(NULL)
  }
  raw_text <- content(resp, "text", encoding = "UTF-8")
  parsed   <- tryCatch(fromJSON(raw_text, simplifyVector = FALSE),
                       error = function(e) { warning("JSON parse error: ", e$message); NULL })
  list(parsed = parsed, raw_text = raw_text)
}

# Paginated fetch — GeoMet caps at 10 000 records per request
fetch_wsc_paginated <- function(station_no, param, page_limit = 10000) {
  param <- toupper(param)
  stopifnot(param %in% c("DISCHARGE", "LEVEL"))
  
  all_rows <- list()
  offset   <- 0
  repeat {
    url <- paste0(BASE_URL,
                  "?STATION_NUMBER=", URLencode(station_no, reserved = TRUE),
                  "&limit=",  page_limit,
                  "&offset=", offset)
    res <- safe_get_json(url)
    if (is.null(res) || is.null(res$parsed)) break
    
    features <- res$parsed$features
    if (is.null(features) || length(features) == 0) break
    
    page_rows <- map_dfr(features, function(f) {
      props <- lapply(f$properties, function(x) if (is.null(x)) NA else x)
      as_tibble(props)
    })
    names(page_rows) <- toupper(names(page_rows))
    
    # Normalise date column
    if (!"DATE" %in% names(page_rows) && "DATETIME" %in% names(page_rows))
      page_rows <- rename(page_rows, DATE = DATETIME)
    if (!"DATE" %in% names(page_rows)) break
    
    # Keep only rows for the requested parameter
    if (!param %in% names(page_rows)) {
      warning(sprintf("'%s' not in response for %s", param, station_no))
      break
    }
    
    page_rows <- page_rows %>%
      select(Date = DATE, Value = all_of(param)) %>%
      mutate(Date = as.Date(str_sub(as.character(Date), 1, 10))) %>%
      filter(!is.na(Value))
    
    all_rows <- c(all_rows, list(page_rows))
    
    # If fewer rows than the page limit, we've reached the end
    if (length(features) < page_limit) break
    offset <- offset + page_limit
    message(sprintf("    [%s / %s] fetched %d rows so far ...",
                    station_no, param, offset))
  }
  
  if (length(all_rows) == 0) return(NULL)
  bind_rows(all_rows) %>% arrange(Date) %>% distinct(Date, .keep_all = TRUE)
}

# ── 4  Download loop ─────────────────────────────────────────────────────────────
download_station <- function(STATION_NUMBER, STATION_NAME,
                             GRDC_NO, GRDC_RIVER,
                             grdc_d_start, grdc_d_end, grdc_d_yrs) {
  message(sprintf("\n[%s]  %s", STATION_NUMBER, STATION_NAME))
  message(sprintf("  GRDC %d (%s) — reference span: %s → %s (%d yrs)",
                  GRDC_NO, GRDC_RIVER, grdc_d_start, grdc_d_end, grdc_d_yrs))
  
  result <- list(station_no   = STATION_NUMBER,
                 station_name = STATION_NAME,
                 grdc_no      = GRDC_NO)
  
  for (param in c("DISCHARGE", "LEVEL")) {
    message(sprintf("  Fetching %s ...", param))
    df <- fetch_wsc_paginated(STATION_NUMBER, param)
    
    if (!is.null(df) && nrow(df) > 0) {
      col_name <- if (param == "DISCHARGE") "Discharge_cms" else "Level_m"
      df       <- rename(df, !!col_name := Value)
      
      out_file <- file.path(output_dir, sprintf("%s_WSC_%s.csv", STATION_NUMBER, param))
      write_csv(df, out_file)
      message(sprintf("  -> %d rows saved: %s", nrow(df), out_file))
      message(sprintf("     WSC span: %s  →  %s  (%d yrs)",
                      min(df$Date), max(df$Date),
                      as.integer(format(max(df$Date), "%Y")) -
                        as.integer(format(min(df$Date), "%Y")) + 1))
      
      result[[paste0(tolower(param), "_rows")]]  <- nrow(df)
      result[[paste0(tolower(param), "_start")]] <- min(df$Date)
      result[[paste0(tolower(param), "_end")]]   <- max(df$Date)
    } else {
      message(sprintf("  -> WARNING: no %s data found", param))
      result[[paste0(tolower(param), "_rows")]]  <- 0L
      result[[paste0(tolower(param), "_start")]] <- NA_Date_
      result[[paste0(tolower(param), "_end")]]   <- NA_Date_
    }
  }
  result
}

results_raw <- pmap(stations_df, download_station)

# ── 5  Validation against GRDC time spans ──────────────────────────────────────
message(strrep("=", 65))
message("VALIDATION: WSC discharge span vs GRDC reference span")
message(strrep("=", 65))

validation <- stations_df %>%
  mutate(
    wsc_start = map_vec(results_raw, ~ .x$discharge_start),
    wsc_end   = map_vec(results_raw, ~ .x$discharge_end),
    wsc_rows  = map_int(results_raw, ~ .x$discharge_rows %||% 0L)
  ) %>%
  mutate(
    wsc_yrs         = as.integer(format(wsc_end, "%Y")) -
      as.integer(format(wsc_start, "%Y")) + 1,
    covers_start    = !is.na(wsc_start) & wsc_start <= grdc_d_start,
    covers_end      = !is.na(wsc_end)   & wsc_end   >= grdc_d_end,
    passes          = covers_start & covers_end,
    start_gap_days  = as.integer(wsc_start - grdc_d_start),  # +ve = WSC starts later
    end_gap_days    = as.integer(wsc_end   - grdc_d_end)     # +ve = WSC ends later
  )

# Console summary
for (i in seq_len(nrow(validation))) {
  v <- validation[i, ]
  status <- if (isTRUE(v$passes)) "✅ PASS" else "❌ FAIL"
  message(sprintf(
    "\n%s  %s — %s",
    status, v$STATION_NUMBER, v$STATION_NAME))
  message(sprintf(
    "  GRDC ref : %s → %s  (%d yrs)",
    v$grdc_d_start, v$grdc_d_end, v$grdc_d_yrs))
  if (!is.na(v$wsc_start)) {
    message(sprintf(
      "  WSC data : %s → %s  (%d yrs,  %d rows)",
      v$wsc_start, v$wsc_end, v$wsc_yrs, v$wsc_rows))
    message(sprintf(
      "  Start gap: %+d days  |  End gap: %+d days  [negative = WSC extends further]",
      v$start_gap_days, v$end_gap_days))
    if (!v$covers_start)
      message(sprintf("  ⚠  WSC starts %d days AFTER GRDC d_start",  v$start_gap_days))
    if (!v$covers_end)
      message(sprintf("  ⚠  WSC ends   %d days BEFORE GRDC d_end",  -v$end_gap_days))
  } else {
    message("  ⚠  No WSC discharge data retrieved — cannot validate")
  }
}

# ── 6  Save validation report ───────────────────────────────────────────────────
report <- validation %>%
  select(
    STATION_NUMBER, STATION_NAME, GRDC_NO, GRDC_RIVER,
    grdc_d_start, grdc_d_end, grdc_d_yrs,
    wsc_start, wsc_end, wsc_yrs, wsc_rows,
    covers_start, covers_end, passes,
    start_gap_days, end_gap_days
  ) %>%
  mutate(across(c(grdc_d_start, grdc_d_end, wsc_start, wsc_end), as.character))

report_file <- file.path(report_dir, "GRDC_WSC_coverage_validation.csv")
write_csv(report, report_file)
message(strrep("-", 65))
message("Validation report saved: ", report_file)

n_pass <- sum(validation$passes, na.rm = TRUE)
n_fail <- nrow(validation) - n_pass
message(sprintf("Summary: %d PASS  |  %d FAIL  (out of %d stations)",
                n_pass, n_fail, nrow(validation)))
message(strrep("=", 65))
message("✅  Download & validation complete.")