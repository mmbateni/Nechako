# Climate Normals scraper (multi-period) with name-based stnID resolver and custom output path
# Requires: tidyverse, httr, rvest, xml2
library(tidyverse)
library(httr)
library(rvest)
library(xml2)

# --- CONFIG: change these to your exact paths/values ---
stations_file <- "Stations.csv"
output_dir    <- "D:/baten/Nechako_Drought_Project/Meteorology/Normaldata_retrievalECCCanada/normals_downloads"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

periods <- c("1991_2020", "1981_2010", "1971_2000", "1961_1990", "1951_1980", "1941_1970", "1931_1960")

# polite pause and retry settings
pause_seconds <- 0.8
max_retries <- 3
retry_backoff <- function(attempt) 1 + attempt * 1.5

# user agent (include contact if you like)
ua <- user_agent("R (rvest) - climate normals fetch - contact: your_email@example.com")

# --- helper functions ---

# Build normals page URL for a given period token and stnID
period_to_page <- function(period_token) paste0("results_", period_token, "_e.html")
build_normals_url <- function(stnID, period_token, prov = "BC") {
  page <- period_to_page(period_token)
  base <- paste0("https://climate.weather.gc.ca/climate_normals/", page)
  paste0(base,
         "?stnID=", URLencode(as.character(stnID)),
         "&prov=", URLencode(prov),
         "&stationName=&searchType=stnName&txtCentralLatMin=&txtCentralLatSec=&txtCentralLongMin=&txtCentralLongSec=&searchMethod=contains")
}

# Safe GET with retries for transient errors
safe_get <- function(url, max_retries = 3, ua = ua) {
  attempt <- 1
  repeat {
    resp <- tryCatch(GET(url, ua, timeout(60)), error = function(e) e)
    if (inherits(resp, "error")) {
      if (attempt >= max_retries) return(list(ok = FALSE, error = conditionMessage(resp), resp = NULL))
      Sys.sleep(retry_backoff(attempt))
      attempt <- attempt + 1
      next
    }
    # HTTP-level handling
    sc <- status_code(resp)
    if (sc >= 500 && attempt < max_retries) {
      Sys.sleep(retry_backoff(attempt))
      attempt <- attempt + 1
      next
    }
    return(list(ok = TRUE, resp = resp, status = sc))
  }
}

# Parse normals page: return temp and precip tables (raw)
parse_normals_page <- function(html_text) {
  doc <- read_html(html_text)
  find_table_by_heading <- function(keyword) {
    headings <- doc %>% html_nodes(xpath = sprintf("//*[contains(translate(text(), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), '%s')]", tolower(keyword)))
    if (length(headings) == 0) return(NULL)
    for (h in headings) {
      tbl <- h %>% html_node(xpath = "following::table[1]")
      if (!is.null(tbl)) return(tbl)
    }
    NULL
  }
  temp_tbl_node <- find_table_by_heading("monthly normals") %||% find_table_by_heading("temperature")
  precip_tbl_node <- find_table_by_heading("precipitation") %||% find_table_by_heading("precipitation normals")
  temp_df <- if (!is.null(temp_tbl_node)) tryCatch(temp_tbl_node %>% html_table(fill = TRUE) %>% as_tibble(), error = function(e) NULL) else NULL
  precip_df <- if (!is.null(precip_tbl_node)) tryCatch(precip_tbl_node %>% html_table(fill = TRUE) %>% as_tibble(), error = function(e) NULL) else NULL
  list(temp = temp_df, precip = precip_df, doc = doc)
}

# Tidy table into long form if it has a MONTH column
tidy_normals_table <- function(df) {
  if (is.null(df)) return(NULL)
  df <- df %>% select(where(~ !all(is.na(.)))) %>% filter(if_any(everything(), ~ !is.na(.)))
  names(df) <- names(df) %>% str_trim() %>% str_replace_all("\\s+", "_") %>% toupper()
  if ("MONTH" %in% names(df)) {
    long <- df %>% pivot_longer(-MONTH, names_to = "STAT", values_to = "VALUE") %>% mutate(MONTH = str_trim(MONTH))
    return(long)
  }
  df
}

# Search the Climate Normals site by station name to find candidate stnID(s)
# This function queries the site search page and extracts links with stnID parameters.
resolve_stn_by_name <- function(station_name, prov = "BC") {
  # Build a search URL that returns a list of stations (the site uses a search page; we try a few common endpoints)
  # Primary attempt: use the normals search results page for the latest period (1981_2010) and search by stationName param
  search_base <- "https://climate.weather.gc.ca/climate_normals/search_e.html"
  # The site accepts 'searchType=stnName' and 'txtStationName' in some forms; try both patterns
  q1 <- paste0("https://climate.weather.gc.ca/climate_normals/results_1981_2010_e.html?stationName=", URLencode(station_name), "&searchType=stnName&prov=", prov)
  q2 <- paste0(search_base, "?searchType=stnName&txtStationName=", URLencode(station_name), "&prov=", prov)
  candidates <- character(0)
  for (q in c(q1, q2)) {
    res <- safe_get(q, max_retries = 2)
    if (!res$ok) next
    page <- read_html(content(res$resp, "text", encoding = "UTF-8"))
    # find links that include "stnID="
    links <- page %>% html_nodes("a") %>% html_attr("href") %>% na.omit()
    stn_links <- links[grepl("stnID=", links)]
    if (length(stn_links) > 0) {
      # extract numeric stnID values
      ids <- unique(gsub(".*stnID=([0-9]+).*", "\\1", stn_links))
      candidates <- c(candidates, ids)
    }
    if (length(candidates) > 0) break
  }
  unique(candidates)
}

# --- Load stations ---
stations <- read_csv(stations_file, show_col_types = FALSE) %>%
  filter(!is.na(STATION_ID)) %>%
  distinct(STATION_ID, STATION_NAME, .keep_all = TRUE) %>%
  mutate(STATION_ID = as.character(STATION_ID))

message("Stations to process: ", nrow(stations))
message("Periods: ", paste(periods, collapse = ", "))

# --- Main loop ---
results <- tibble(STATION_ID = character(), STATION_NAME = character(), PERIOD = character(),
                  used_stnID = character(), temp_rows = integer(), precip_rows = integer(), note = character())

for (i in seq_len(nrow(stations))) {
  st <- stations[i, ]
  st_id <- st$STATION_ID
  st_name <- st$STATION_NAME
  
  # candidate stnIDs: start with the provided STATION_ID
  candidate_ids <- c(st_id)
  
  # We'll cache resolved IDs per station to avoid repeated name searches
  resolved_cache <- character(0)
  
  for (period in periods) {
    message(sprintf("[%d/%d] %s (%s) period %s", i, nrow(stations), st_name, st_id, period))
    
    success <- FALSE
    last_note <- NA_character_
    used_id <- NA_character_
    
    # iterate candidate stnIDs (first the original, then any resolved ones)
    for (cid in unique(candidate_ids)) {
      used_id <- cid
      url <- build_normals_url(cid, period)
      res <- safe_get(url, max_retries = max_retries)
      if (!res$ok) {
        # request-level failure
        last_note <- paste0("request_error:", res$error %||% "unknown")
        next
      }
      # HTTP status handling
      sc <- res$status
      if (sc == 404) {
        last_note <- "http_404"
        next
      }
      if (sc >= 400) {
        last_note <- paste0("http_", sc)
        next
      }
      # parse page
      html_text <- content(res$resp, "text", encoding = "UTF-8")
      parsed <- parse_normals_page(html_text)
      temp_df <- tidy_normals_table(parsed$temp)
      precip_df <- tidy_normals_table(parsed$precip)
      
      saved_temp_rows <- if (!is.null(temp_df)) nrow(temp_df) else 0L
      saved_precip_rows <- if (!is.null(precip_df)) nrow(precip_df) else 0L
      
      if ((saved_temp_rows + saved_precip_rows) > 0) {
        # save files
        if (!is.null(temp_df) && nrow(temp_df) > 0) {
          temp_out_file <- file.path(output_dir, paste0(st_id, "_NORMALS_TEMP_", period, "_stnID_", cid, ".csv"))
          write_csv(temp_df, temp_out_file)
        }
        if (!is.null(precip_df) && nrow(precip_df) > 0) {
          precip_out_file <- file.path(output_dir, paste0(st_id, "_NORMALS_PRECIP_", period, "_stnID_", cid, ".csv"))
          write_csv(precip_df, precip_out_file)
        }
        results <- bind_rows(results, tibble(STATION_ID = st_id, STATION_NAME = st_name, PERIOD = period,
                                             used_stnID = cid, temp_rows = saved_temp_rows, precip_rows = saved_precip_rows, note = "ok"))
        success <- TRUE
        break
      } else {
        last_note <- "no_tables_found"
        # continue to next candidate id
      }
    } # end candidate loop
    
    if (!success) {
      # attempt to resolve by station name (only once per station)
      if (length(resolved_cache) == 0) {
        message("  -> Attempting name-based stnID resolution for station: ", st_name)
        ids_found <- resolve_stn_by_name(st_name, prov = "BC")
        if (length(ids_found) > 0) {
          message("  -> Resolved candidate stnID(s): ", paste(ids_found, collapse = ", "))
          candidate_ids <- c(candidate_ids, ids_found)
          resolved_cache <- ids_found
          # retry this period with new candidates
          # decrement loop index to retry same period with new candidates
          # (we implement by re-running the candidate loop above via a simple repeat)
          # For simplicity, re-run candidate loop once more here:
          for (cid in unique(candidate_ids)) {
            used_id <- cid
            url <- build_normals_url(cid, period)
            res <- safe_get(url, max_retries = max_retries)
            if (!res$ok) {
              last_note <- paste0("request_error:", res$error %||% "unknown")
              next
            }
            sc <- res$status
            if (sc == 404) {
              last_note <- "http_404"
              next
            }
            if (sc >= 400) {
              last_note <- paste0("http_", sc)
              next
            }
            html_text <- content(res$resp, "text", encoding = "UTF-8")
            parsed <- parse_normals_page(html_text)
            temp_df <- tidy_normals_table(parsed$temp)
            precip_df <- tidy_normals_table(parsed$precip)
            saved_temp_rows <- if (!is.null(temp_df)) nrow(temp_df) else 0L
            saved_precip_rows <- if (!is.null(precip_df)) nrow(precip_df) else 0L
            if ((saved_temp_rows + saved_precip_rows) > 0) {
              if (!is.null(temp_df) && nrow(temp_df) > 0) {
                temp_out_file <- file.path(output_dir, paste0(st_id, "_NORMALS_TEMP_", period, "_stnID_", cid, ".csv"))
                write_csv(temp_df, temp_out_file)
              }
              if (!is.null(precip_df) && nrow(precip_df) > 0) {
                precip_out_file <- file.path(output_dir, paste0(st_id, "_NORMALS_PRECIP_", period, "_stnID_", cid, ".csv"))
                write_csv(precip_df, precip_out_file)
              }
              results <- bind_rows(results, tibble(STATION_ID = st_id, STATION_NAME = st_name, PERIOD = period,
                                                   used_stnID = cid, temp_rows = saved_temp_rows, precip_rows = saved_precip_rows, note = "ok_resolved"))
              success <- TRUE
              break
            } else {
              last_note <- "no_tables_found_after_resolve"
            }
          } # end retry with resolved ids
        } else {
          last_note <- "resolve_none_found"
        }
      } # end name-resolve block
    }
    
    if (!success) {
      results <- bind_rows(results, tibble(STATION_ID = st_id, STATION_NAME = st_name, PERIOD = period,
                                           used_stnID = used_id %||% NA_character_, temp_rows = 0L, precip_rows = 0L, note = last_note %||% "failed"))
    }
    
    Sys.sleep(pause_seconds)
  } # end period loop
} # end station loop

# --- Save summary ---
summary_file <- file.path(output_dir, "normals_fetch_summary.csv")
write_csv(results, summary_file)
message("Done. Summary written to: ", summary_file)
message("Saved normals CSVs (if any) to: ", normalizePath(output_dir))
