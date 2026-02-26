# ── Historical WSC Data via tidyhydat (HYDAT database) ──────────────────────
# Saves separate CSV files per station for both daily flow and daily level.
# Coverage goes back to the 1900s — far beyond the GeoMet API (~2011 limit).
# ────────────────────────────────────────────────────────────────────────────

library(tidyhydat)
library(tidyverse)

# ── 0  Paths ──────────────────────────────────────────────────────────────────
setwd("D:/Nechako_Drought/Nechako/Hydrology/data_retrievalWaterSurveyofCanada")

flow_dir  <- "HYDAT_downloads/flow"
level_dir <- "HYDAT_downloads/level"
dir.create(flow_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(level_dir, recursive = TRUE, showWarnings = FALSE)

message("Flow  output : ", flow_dir)
message("Level output : ", level_dir)
message(strrep("-", 65))

# ── 1  One-time HYDAT database download (~1.5 GB) ─────────────────────────────
# Comment this out after the first run — the database persists locally.
if (!file.exists(hy_default_db())) {
  message("Downloading HYDAT database (one-time, ~1.5 GB) ...")
  download_hydat()
  message("HYDAT download complete.")
} else {
  message("HYDAT database already present: ", hy_default_db())
}

# ── 2  Target stations ────────────────────────────────────────────────────────
stations <- tribble(
  ~STATION_NUMBER, ~STATION_NAME,
  "08JA015", "LAVENTIE CREEK NEAR THE MOUTH",
  "08JB002", "STELLAKO RIVER AT GLENANNAN",
  "08JC001", "NECHAKO RIVER AT VANDERHOOF",
  "08JC002", "NECHAKO RIVER AT ISLE PIERRE",
  "08JC005", "CHILAKO RIVER NEAR PRINCE GEORGE",
  "08JE001", "STUART RIVER NEAR FORT ST. JAMES",
  "08JE004", "TSILCOH RIVER NEAR THE MOUTH"
)

station_ids <- stations$STATION_NUMBER
message("Stations to process: ", paste(station_ids, collapse = ", "))
message(strrep("-", 65))

# ── 3  Helper: save one variable (flow or level) per station ──────────────────
save_per_station <- function(data, var_label, out_dir, col_rename) {
  # data      : full data frame returned by hy_daily_flows() / hy_daily_levels()
  # var_label : "DISCHARGE" or "LEVEL" — used in filename
  # out_dir   : destination directory
  # col_rename: named vector c(old = new) for the value column
  
  if (is.null(data) || nrow(data) == 0) {
    message(sprintf("  No %s data returned — skipping.", var_label))
    return(invisible(NULL))
  }
  
  data %>%
    rename(any_of(col_rename)) %>%          # rename Value → Discharge_cms / Level_m
    group_by(STATION_NUMBER) %>%
    group_walk(function(df, grp) {
      stn  <- grp$STATION_NUMBER
      name <- stations$STATION_NAME[stations$STATION_NUMBER == stn]
      df   <- df %>% arrange(Date) %>% select(-any_of("STATION_NUMBER"))
      
      out_file <- file.path(out_dir, sprintf("%s_HYDAT_%s.csv", stn, var_label))
      write_csv(df, out_file)
      
      message(sprintf("  [%s] %s", stn, name))
      message(sprintf("    -> %d rows | %s → %s | saved: %s",
                      nrow(df), min(df$Date), max(df$Date), out_file))
    })
}

# ── 4  Download & save daily flow (discharge) ─────────────────────────────────
message("\nFetching daily FLOW (discharge) from HYDAT ...")
flow_hist <- tryCatch(
  hy_daily_flows(station_number = station_ids),
  error = function(e) { message("ERROR fetching flow: ", e$message); NULL }
)

save_per_station(
  data       = flow_hist,
  var_label  = "DISCHARGE",
  out_dir    = flow_dir,
  col_rename = c(Discharge_cms = "Value")
)

# ── 5  Download & save daily water level ──────────────────────────────────────
message("\nFetching daily LEVEL from HYDAT ...")
level_hist <- tryCatch(
  hy_daily_levels(station_number = station_ids),
  error = function(e) { message("ERROR fetching level: ", e$message); NULL }
)

save_per_station(
  data       = level_hist,
  var_label  = "LEVEL",
  out_dir    = level_dir,
  col_rename = c(Level_m = "Value")
)

# ── 6  Summary table ──────────────────────────────────────────────────────────
message(strrep("=", 65))
message("SUMMARY")
message(strrep("=", 65))

summarise_coverage <- function(data, var_label) {
  if (is.null(data) || nrow(data) == 0) return(tibble())
  data %>%
    group_by(STATION_NUMBER) %>%
    summarise(
      !!paste0(var_label, "_start") := min(Date, na.rm = TRUE),
      !!paste0(var_label, "_end")   := max(Date, na.rm = TRUE),
      !!paste0(var_label, "_rows")  := n(),
      .groups = "drop"
    )
}

summary_tbl <- stations %>%
  left_join(summarise_coverage(flow_hist,  "flow"),  by = "STATION_NUMBER") %>%
  left_join(summarise_coverage(level_hist, "level"), by = "STATION_NUMBER")

print(summary_tbl, n = Inf)

summary_file <- "HYDAT_downloads/HYDAT_coverage_summary.csv"
write_csv(summary_tbl, summary_file)
message("\nCoverage summary saved: ", summary_file)
message(strrep("=", 65))
message("✅  HYDAT download complete.")