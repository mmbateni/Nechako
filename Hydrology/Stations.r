####################################################################################
# Nechako Basin ??? Hydrometric Station Discovery & Kenney Dam Classification
#
# Based on: Stations.r
# Modification: adds data-availability check and classifies every station
# into one of three categories relative to Kenney Dam / Skins Lake Spillway:
#
#   1 ??? UPSTREAM  of Kenney Dam : within the Nechako Reservoir catchment
#       (Ootsa, Whitesail, Tahtsa, Tetachuck, Natalkuz, Eutsuk lakes and
#        their tributaries; dam/spillway gauges). NOTE: Francois Lake,
#        Fraser Lake, the Endako River, and the Nadina River/Lake are NOT
#        part of the reservoir ??? they are an independent natural lake chain
#        and have been moved to category 3 (see below).
#
#   2 ??? DOWNSTREAM of Kenney Dam : on the regulated Nechako River / Cheslatta
#       River system below the Skins Lake Spillway outlet ??? including
#       Cheslatta Lake (sits below the spillway, above Cheslatta Falls),
#       Cheslatta Falls, and onward through Fort Fraser, Vanderhoof, and
#       Isle Pierre ??? whose flow is the direct, dam-controlled release.
#       CORRECTED: Cheslatta Lake (08JA018) was previously misclassified as
#       upstream/reservoir; it is fed by the spillway release, not impounded
#       within the reservoir itself.
#
#   3 ??? NOT AFFECTED by dam : tributaries with independent natural catchments
#       that join the Nechako or Stuart River system downstream, including:
#       Chilako, Driftwood, Pinchi, Tsilcoh, small creeks near Vanderhoof/PG,
#       the Stuart River (independent Stuart Lake watershed), and the
#       Francois Lake / Fraser Lake / Endako River / Nadina River???Lake /
#       Stellako River / Nautley River system (a separate, natural lake
#       chain that joins the Nechako at Fort Fraser ??? not impounded by or
#       regulated by Kenney Dam).
#       North Beach Creek (08JB013) is left in category 1 provisionally;
#       its relationship to the dam vs. the Francois/Fraser system is
#       unconfirmed ??? revisit if more information becomes available.
#
# Kenney Dam coordinates : ~53.58 N, 124.97 W  (station 08JA022)
# Skins Lake Spillway    : ~53.77 N, 125.99 W  (station 08JA013 / 08JA023)
#
# NOTE: all dplyr verbs are explicitly namespaced (dplyr::select, etc.) because
# raster and other spatial packages loaded in this R session mask dplyr::select(),
# dplyr::filter(), etc. and produce "unused arguments" errors without namespacing.
####################################################################################

library(tidyhydat)
library(dplyr)
library(tidyr)
library(readr)
# Set working directory to the script's location
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
} else {
  args <- commandArgs(trailingOnly = FALSE)
  script_path <- sub("--file=", "", args[grep("--file=", args)])
  if (length(script_path) > 0) setwd(dirname(script_path))
}
# ?????? Download HYDAT database the first time (uncomment if needed) ??????????????????????????????????????????
# download_hydat()

# ?????? 1  All Nechako basin stations (WSC sub-basin prefix 08J) ???????????????????????????????????????????????????
nechako_gauges <- hy_stations(prov_terr_state_loc = "BC") %>%
  dplyr::filter(grepl("^08J", STATION_NUMBER))

message("Total 08J stations in HYDAT : ", nrow(nechako_gauges))

# ?????? 2  Data availability: keep stations that have discharge (Q) or level (H) ???
data_range <- hy_stn_data_range(station_number = nechako_gauges$STATION_NUMBER)

avail <- data_range %>%
  dplyr::filter(DATA_TYPE %in% c("Q", "H")) %>%
  dplyr::select(STATION_NUMBER, DATA_TYPE, Year_from, Year_to) %>%
  tidyr::pivot_wider(
    id_cols     = STATION_NUMBER,
    names_from  = DATA_TYPE,
    values_from = c(Year_from, Year_to),
    values_fn   = first          # guard against duplicate rows per station/type
  ) %>%
  dplyr::rename(
    flow_start  = Year_from_Q,
    flow_end    = Year_to_Q,
    level_start = Year_from_H,
    level_end   = Year_to_H
  )

stations_with_data <- nechako_gauges %>%
  dplyr::inner_join(avail, by = "STATION_NUMBER") %>%
  dplyr::filter(!is.na(flow_start) | !is.na(level_start))

message("Stations with discharge or level data : ", nrow(stations_with_data))
message(strrep("-", 72))

# ?????? 3  Station classification ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
#
# ?????? UPSTREAM (Nechako Reservoir watershed) ????????????????????????????????????????????????????????????????????????????????????????????????????????????
# Includes: reservoir lake gauges, all tributaries flowing INTO the reservoir
# (Ootsa / Whitesail / Tahtsa / Natalkuz / Eutsuk arms), and the dam / spillway
# infrastructure gauges. Francois Lake, Fraser Lake, Endako River, and Nadina
# River/Lake are EXCLUDED ??? see category 3 below.
upstream_ids <- c(
  # Reservoir inflows ??? western arms (Ootsa / Whitesail / Tahtsa group)
  "08JA002",  # OOTSA RIVER AT OOTSA LAKE
  "08JA003",  # WHITESAIL RIVER NEAR OOTSA LAKE
  "08JA004",  # TETACHUCK RIVER NEAR OOTSA LAKE
  "08JA005",  # TAHTSA RIVER NEAR OOTSA LAKE
  "08JA006",  # TAHTSA RIVER AT OUTLET OF TAHTSA LAKE
  "08JA021",  # COLES CREEK ABOVE TROITSA CREEK          (~127.2 W, Tahtsa arm)
  "08JA025",  # KASALKA CREEK AT THE MOUTH               (~127.1 W, Tahtsa arm)
  "08JA026",  # MOUNT BAPTISTE CREEK BELOW 1040 M        (~127.4 W, Tahtsa arm)
  "08JA027",  # NECHAKO RESERVOIR AT TAHTSA REACH        (reservoir gauge)
  "08JA029",  # WHITESAIL MIDDLE CREEK NEAR TAHTSA REACH (~127.0 W)
  "08JA030",  # TAHTSA LAKE NEAR KEMANO                  (Kemano diversion arm)
  "08JA015",  # LAVENTIE CREEK NEAR THE MOUTH            (~127.5 W, reservoir)
  "08JA016",  # MACIVOR CREEK NEAR THE MOUTH             (~126.4 W, reservoir)
  "08JA028",  # EUTSUK RIVER AT OUTLET OF EUTSUK LAKE    (~126.1 W, reservoir)
  "08JA020",  # CHELASLIE RIVER NEAR THE MOUTH           (~125.9 W, reservoir)
  # Nechako River within the reservoir (Natalkuz / Ootsa area)
  "08JA007",  # NECHAKO RIVER AT OUTLET OF NATALKUZ LAKE (~125.1 W, in reservoir)
  "08JA011",  # NECHAKO RIVER NEAR OOTSA LAKE            (within reservoir)
  # Cheslatta / Entiako area (west of dam, part of reservoir system)
  # NOTE: Cheslatta Lake (08JA018) moved OUT of this list ??? see downstream_ids
  # below. It lies below the Skins Lake Spillway outlet (on the Cheslatta
  # River, above Cheslatta Falls), so it receives the dam-controlled spillway
  # release rather than being impounded within the reservoir itself.
  "08JA009",  # CHESLATTA RIVER NEAR OOTSA LAKE
  "08JA014",  # VAN TINE CREEK NEAR THE MOUTH            (~125.4 W)
  "08JA024",  # ENTIAKO RIVER NEAR THE MOUTH             (~125.4 W, reservoir)
  # Dam / spillway infrastructure gauges
  "08JA013",  # SKINS LAKE SPILLWAY, NECHAKO RESERVOIR   (spillway flow ??? ACTIVE)
  "08JA022",  # NECHAKO RESERVOIR AT KENNEY DAM          (dam gauge)
  "08JA023",  # NECHAKO RESERVOIR AT SKINS LAKE SPILLWAY (reservoir level ??? ACTIVE)
  # NOTE: North Beach Creek's status relative to the dam is unconfirmed ???
  # kept here provisionally; revisit if you find it actually drains to the
  # independent Francois/Fraser system below instead of the reservoir proper.
  "08JB013"   # NORTH BEACH CREEK ABOVE ALLIN CREEK      (drains to Burns/Francois ??? UNCONFIRMED)
)

# ?????? DOWNSTREAM (regulated Nechako River ??? direct dam-controlled release) ??????????????????
# Includes: the Cheslatta River/Lake reach immediately below the Skins Lake
# Spillway outlet, then the Nechako mainstem from Cheslatta Falls east through
# Fort Fraser, Vanderhoof, and Isle Pierre ??? the entire reach that directly
# carries Kenney Dam / Skins Lake Spillway releases.
downstream_ids <- c(
  # Cheslatta Lake ??? immediately below the Skins Lake Spillway outlet,
  # upstream of Cheslatta Falls. CORRECTED: previously misclassified as
  # upstream/reservoir; it is fed by the spillway release itself.
  "08JA018",  # CHESLATTA LAKE AT WEST END            (~125.6 W, below spillway)
  # Nechako mainstem ??? near-dam and upper regulated reach
  "08JA010",  # NECHAKO RIVER BELOW BIG BEND CREEK   (~124.86 W, east of dam)
  "08JA017",  # NECHAKO RIVER BELOW CHESLATTA FALLS   (ACTIVE ??? regulated release point)
  "08JA019",  # CHEDAKUZ CREEK AT THE MOUTH           (joins Nechako just below dam)
  "08JA001",  # NECHAKO RIVER AT FORT FRASER
  # Nechako mainstem ??? Vanderhoof to Isle Pierre
  "08JC001",  # NECHAKO RIVER AT VANDERHOOF           (ACTIVE)
  "08JC002"   # NECHAKO RIVER AT ISLE PIERRE          (ACTIVE)
)

# Everything else = NOT AFFECTED (independent catchments joining Nechako/Stuart
# downstream; their own flows are entirely natural and unregulated by the dam).
#
# NOTE: Stellako River (08JB002), Nautley River (08JB003), and Stuart River
# (08JE001) were moved OUT of downstream_ids. Stellako/Nautley convey the
# natural, independent Francois Lake -> Fraser Lake outflow into the Nechako
# at Fort Fraser ??? that lake chain is not part of the Nechako Reservoir behind
# Kenney Dam, so its discharge is not dam-regulated. Stuart River drains the
# separate Stuart Lake watershed and joins the Nechako downstream, but its
# flow is likewise entirely natural. All three now fall through to category 3.

# ?????? Apply classification ??????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
stations_classified <- stations_with_data %>%
  dplyr::mutate(
    dam_category = dplyr::case_when(
      STATION_NUMBER %in% upstream_ids   ~ "1 - Upstream of Kenney Dam",
      STATION_NUMBER %in% downstream_ids ~ "2 - Downstream of Kenney Dam",
      TRUE                                ~ "3 - Not Affected by Dam"
    )
  ) %>%
  dplyr::arrange(dam_category, STATION_NUMBER)

# ?????? 4  Console summary ????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
cat("\n", strrep("=", 72), "\n", sep = "")
cat("NECHAKO BASIN ??? Hydrometric Stations with Data vs Kenney Dam Position\n")
cat(strrep("=", 72), "\n\n")

category_labels <- c(
  "1 - Upstream of Kenney Dam",
  "2 - Downstream of Kenney Dam",
  "3 - Not Affected by Dam"
)

for (lbl in category_labels) {
  sub    <- dplyr::filter(stations_classified, dam_category == lbl)
  header <- gsub("^[0-9]+ - ", "", lbl)
  cat(sprintf("?????? %s  (%d stations) ??????\n", header, nrow(sub)))
  
  if (nrow(sub) == 0) { cat("  (none)\n\n"); next }
  
  for (i in seq_len(nrow(sub))) {
    r      <- sub[i, ]
    status <- if (!is.na(r$HYD_STATUS) && r$HYD_STATUS == "ACTIVE") "???" else "???"
    q_str  <- if (!is.na(r$flow_start))  sprintf("Q:%d-%d", r$flow_start,  r$flow_end)  else "Q:---"
    h_str  <- if (!is.na(r$level_start)) sprintf("H:%d-%d", r$level_start, r$level_end) else "H:---"
    cat(sprintf("  %s %s  %-48s  %s  %s\n",
                status, r$STATION_NUMBER, r$STATION_NAME, q_str, h_str))
  }
  cat("  (??? ACTIVE  ??? DISCONTINUED)\n\n")
}

# ?????? 5  Summary counts ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
n_up   <- sum(stations_classified$dam_category == "1 - Upstream of Kenney Dam")
n_down <- sum(stations_classified$dam_category == "2 - Downstream of Kenney Dam")
n_none <- sum(stations_classified$dam_category == "3 - Not Affected by Dam")

message(strrep("-", 72))
message(sprintf(
  "Total with data: %d  |  Upstream: %d  |  Downstream: %d  |  Not affected: %d",
  nrow(stations_classified), n_up, n_down, n_none))
message(strrep("=", 72))

# ?????? 6  Save classified station list ??????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
out_df <- stations_classified %>%
  dplyr::select(
    dam_category,
    STATION_NUMBER, STATION_NAME, HYD_STATUS,
    LATITUDE, LONGITUDE, DRAINAGE_AREA_GROSS,
    flow_start, flow_end,
    level_start, level_end
  )

write_csv(out_df, "nechako_stations_classified.csv")
message("Saved: nechako_stations_classified.csv")

# ?????? 7  Drought-suitable subset (2022???2024 analysis) ??????????????????????????????????????????????????????????????????????????????
#
# Criteria for a station to be usable in the 2022-2024 drought investigation:
#   (a) at least 10 years of record (combined span of whichever of Q/H is
#       longer ??? using the earliest start and latest end across the two
#       data types kept on the station), AND
#   (b) the record fully COVERS the 2022-2024 period, i.e. starts at or
#       before 2022 and ends at or after 2024 (not just overlapping it).
#
# If your HYDAT copy hasn't been updated through 2024 yet, most stations will
# fail criterion (b) on the "end >= 2024" side ??? rerun download_hydat() to
# refresh, or switch to the "overlap" version (commented below) if partial
# coverage is acceptable for your analysis.

drought_start <- 2022
drought_end   <- 2024
min_years     <- 10

stations_drought <- stations_classified %>%
  dplyr::mutate(
    rec_start  = pmin(flow_start, level_start, na.rm = TRUE),
    rec_end    = pmax(flow_end,   level_end,   na.rm = TRUE),
    n_years    = rec_end - rec_start + 1,
    has_10yr   = n_years >= min_years,
    # Full coverage of 2022-2024:
    covers_2022_2024 = rec_start <= drought_start & rec_end >= drought_end
    # Overlap-only alternative (uncomment to use instead of full coverage):
    # covers_2022_2024 = rec_start <= drought_end & rec_end >= drought_start
  ) %>%
  dplyr::filter(has_10yr, covers_2022_2024) %>%
  dplyr::arrange(dam_category, STATION_NUMBER)

message(strrep("-", 72))
message("Stations with >= ", min_years, " years of data AND covering ",
        drought_start, "-", drought_end, " : ", nrow(stations_drought))

n_up_d   <- sum(stations_drought$dam_category == "1 - Upstream of Kenney Dam")
n_down_d <- sum(stations_drought$dam_category == "2 - Downstream of Kenney Dam")
n_none_d <- sum(stations_drought$dam_category == "3 - Not Affected by Dam")

message(sprintf(
  "  Upstream: %d  |  Downstream: %d  |  Not affected: %d",
  n_up_d, n_down_d, n_none_d))
message(strrep("-", 72))

for (lbl in category_labels) {
  sub    <- dplyr::filter(stations_drought, dam_category == lbl)
  header <- gsub("^[0-9]+ - ", "", lbl)
  cat(sprintf("?????? %s  (%d stations, drought-eligible) ??????\n", header, nrow(sub)))
  
  if (nrow(sub) == 0) { cat("  (none)\n\n"); next }
  
  for (i in seq_len(nrow(sub))) {
    r      <- sub[i, ]
    status <- if (!is.na(r$HYD_STATUS) && r$HYD_STATUS == "ACTIVE") "???" else "???"
    cat(sprintf("  %s %s  %-48s  %d-%d  (%d yrs)\n",
                status, r$STATION_NUMBER, r$STATION_NAME,
                r$rec_start, r$rec_end, r$n_years))
  }
  cat("\n")
}

out_drought <- stations_drought %>%
  dplyr::select(
    dam_category,
    STATION_NUMBER, STATION_NAME, HYD_STATUS,
    LATITUDE, LONGITUDE, DRAINAGE_AREA_GROSS,
    flow_start, flow_end,
    level_start, level_end,
    rec_start, rec_end, n_years
  )

write_csv(out_drought, "nechako_stations_drought_2022_2024.csv")
message("Saved: nechako_stations_drought_2022_2024.csv")

# ?????? 8  Missing-data audit (days/months) for drought-candidate stations ?????????????????????
#
# Rather than checking just the 2022-2024 window, this audits each candidate's
# FULL native period of record ??? its own flow_start..flow_end for discharge (Q)
# and its own level_start..level_end for level (H), since the two series on a
# station don't always span the same years. For each, it pulls the ACTUAL
# daily values from HYDAT and reports how many days are missing and how many
# distinct calendar months contain at least one gap. A day counts as "missing"
# if there is no row for that date in HYDAT, or the row's Value is NA.
# Output is one compact line per data type per station (no exhaustive
# month-by-month listing ??? see the saved CSV if you need the full date list).

flow_period  <- stations_drought %>%
  dplyr::filter(!is.na(flow_start)) %>%
  dplyr::select(STATION_NUMBER, p_start = flow_start, p_end = flow_end)

level_period <- stations_drought %>%
  dplyr::filter(!is.na(level_start)) %>%
  dplyr::select(STATION_NUMBER, p_start = level_start, p_end = level_end)

# Fetch once across the superset of years needed, then check each station
# against its OWN native start/end within that fetched data.
fetch_daily <- function(period_df, fetch_fn) {
  if (nrow(period_df) == 0) {
    return(tibble::tibble(STATION_NUMBER = character(), Date = as.Date(character()), Value = numeric()))
  }
  fetch_fn(
    station_number = period_df$STATION_NUMBER,
    start_date = as.Date(sprintf("%d-01-01", min(period_df$p_start))),
    end_date   = as.Date(sprintf("%d-12-31", max(period_df$p_end)))
  )
}

daily_q <- fetch_daily(flow_period,  tidyhydat::hy_daily_flows)
daily_h <- fetch_daily(level_period, tidyhydat::hy_daily_levels)

# One compact row per station: expected/missing days over ITS OWN native
# period, plus how many distinct months that missingness touches.
missing_summary <- function(daily_df, period_df) {
  dplyr::bind_rows(lapply(seq_len(nrow(period_df)), function(i) {
    sid <- period_df$STATION_NUMBER[i]
    d0  <- as.Date(sprintf("%d-01-01", period_df$p_start[i]))
    d1  <- as.Date(sprintf("%d-12-31", period_df$p_end[i]))
    full_dates <- seq(d0, d1, by = "day")
    
    present <- daily_df %>%
      dplyr::filter(STATION_NUMBER == sid, !is.na(Value)) %>%
      dplyr::pull(Date)
    missing <- full_dates[!full_dates %in% present]
    
    tibble::tibble(
      STATION_NUMBER = sid,
      period         = sprintf("%d-%d", period_df$p_start[i], period_df$p_end[i]),
      n_expected     = length(full_dates),
      n_missing      = length(missing),
      pct_missing    = round(100 * length(missing) / length(full_dates), 1),
      n_months_gappy = length(unique(format(missing, "%Y-%m")))
    )
  }))
}

q_missing <- missing_summary(daily_q, flow_period)
h_missing <- missing_summary(daily_h, level_period)

cat("\n", strrep("=", 72), "\n", sep = "")
cat("MISSING-DATA AUDIT ??? full period of record, per station\n")
cat(strrep("=", 72), "\n\n")

for (lbl in category_labels) {
  sub <- dplyr::filter(stations_drought, dam_category == lbl)
  if (nrow(sub) == 0) next
  header <- gsub("^[0-9]+ - ", "", lbl)
  cat(sprintf("?????? %s ??????\n", header))
  
  for (i in seq_len(nrow(sub))) {
    r <- sub[i, ]
    
    q_str <- if (r$STATION_NUMBER %in% flow_period$STATION_NUMBER) {
      qm <- dplyr::filter(q_missing, STATION_NUMBER == r$STATION_NUMBER)
      sprintf("Q %s: %d/%dmo missing (%.1f%%)", qm$period, qm$n_missing, qm$n_months_gappy, qm$pct_missing)
    } else {
      "Q: ???"
    }
    
    h_str <- if (r$STATION_NUMBER %in% level_period$STATION_NUMBER) {
      hm <- dplyr::filter(h_missing, STATION_NUMBER == r$STATION_NUMBER)
      sprintf("H %s: %d/%dmo missing (%.1f%%)", hm$period, hm$n_missing, hm$n_months_gappy, hm$pct_missing)
    } else {
      "H: ???"
    }
    
    cat(sprintf("  %s  %-48s  %s   %s\n", r$STATION_NUMBER, r$STATION_NAME, q_str, h_str))
  }
  cat("\n")
}
cat(strrep("=", 72), "\n", sep = "")
# Legend: "N/Mmo missing (P%)" = N days missing, spread across M distinct
# calendar months, P% of that series' full native period.

# ?????? 9  Save missing-data audit ????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
out_missing <- stations_drought %>%
  dplyr::select(dam_category, STATION_NUMBER, STATION_NAME) %>%
  dplyr::left_join(
    q_missing %>% dplyr::rename_with(~ paste0("Q_", .x), -STATION_NUMBER),
    by = "STATION_NUMBER"
  ) %>%
  dplyr::left_join(
    h_missing %>% dplyr::rename_with(~ paste0("H_", .x), -STATION_NUMBER),
    by = "STATION_NUMBER"
  )

write_csv(out_missing, "nechako_stations_missing_data_full_period.csv")
message("Saved: nechako_stations_missing_data_full_period.csv")