# ============================================================================
# SURFACE-WATER SUPPLY INDEX (SWSI) — NECHAKO BASIN
# Complete Garen (1993) Implementation with GloLakes v2 Reservoir Storage
#
# Reference: Garen, D.C. (1993). Revised Surface-Water Supply Index for 
#            Western United States. J. Water Resources Planning & Management,
#            119(4), 437-454. https://doi.org/10.1061/(ASCE)0733-9496(1993)119:4(437)
#
# Data Sources:
#   (1) GloLakes v2 — Lake storage 1984-present (Hou et al. 2024, ESSD)
#   (2) ERA5-Land   — SWE & Precipitation (monthly_data_direct/)
#   (3) WSC         — Naturalized flows & discharge data
#   (4) BC RFC      — Snow pillow SWE data (Hydrology/SWE_data/)
#
# Time Span: 1984-2025 (maximized based on GloLakes start date)
# ============================================================================

rm(list = ls())

# ============================================================================
# SECTION 1: PACKAGES
# ============================================================================
packages_needed <- c("tidyverse", "lubridate", "zoo", "ncdf4", "lmomco")
for (pkg in packages_needed) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.rstudio.com/")
    library(pkg, character.only = TRUE)
  }
}

# ============================================================================
# SECTION 2: FILE PATHS & SETTINGS (CORRECTED)
# ============================================================================
setwd("D:/Nechako_Drought/Nechako")

# ── Input Directories (CORRECTED PATHS) ─────────────────────────────────────
INPUT_GLOLAKES_DIR  <- "GloLakes/nechako_extracted"
INPUT_ERA5_DIR      <- "monthly_data_direct"  # ✓ CORRECTED
INPUT_NATURAL_DIR   <- "Hydrology/data_retrievalWaterSurveyofCanada/naturalized_flows"  # ✓ CORRECTED
INPUT_WSC_DIR       <- "Hydrology/data_retrievalWaterSurveyofCanada/data_downloads_geomet_api"
INPUT_SNOWPILLOW_DIR <- "Hydrology/SWE_data"  # ✓ CORRECTED

# ── Output Directory ─────────────────────────────────────────────────────────
MAIN_OUTPUT_DIR <- "swsi_results_complete"
for (sub in c("", "processed_data", "era5_extracted", "forecast_models", 
              "swsi_timeseries", "drought_events", "reports", "diagnostics")) {
  dir.create(file.path(MAIN_OUTPUT_DIR, sub), showWarnings = FALSE, recursive = TRUE)
}

# ── Analysis Parameters ──────────────────────────────────────────────────────
WATER_YEAR_START_MONTH <- 10L
ACCUMULATION_MONTHS    <- c(10, 11, 12, 1, 2, 3)   # Oct-Mar
MELT_MONTHS            <- c(4, 5, 6, 7, 8, 9)       # Apr-Sep
MIN_YEARS_FOR_DIST     <- 10
MIN_YEARS_REGRESSION   <- 15
MIN_COMPLETENESS_PCT   <- 70
MAX_FILL_DAYS          <- 14
MIN_DAYS_PER_MONTH     <- 15

# ── ERA5-Land Basin Bounding Box ─────────────────────────────────────────────
BASIN_LAT_MIN  <- 52.5
BASIN_LAT_MAX  <- 55.5
BASIN_LON_MIN  <- -127.0
BASIN_LON_MAX  <- -120.5
BASIN_AREA_KM2 <- 47000  # Nechako above Isle Pierre

# ── Drought Triggers (Garen p.444) ───────────────────────────────────────────
TRIGGER_WATCH     <- -1.0
TRIGGER_EMERGENCY <- -2.0

# ── Reservoir Data Status ────────────────────────────────────────────────────
RESERVOIR_DATA_AVAILABLE <- TRUE  # ✓ GloLakes extraction provides storage

# ── Primary Station for Streamflow Forecast ──────────────────────────────────
PRIMARY_STATION <- "08JC002"  # Nechako R above Nautley R (naturalized available)

# ============================================================================
# SECTION 3: GLOLAKES RESERVOIR STORAGE LOADING
# ============================================================================
load_reservoir_storage <- function() {
  cat("\n[GloLakes] Loading reservoir storage from best-estimate extraction...\n")
  
  f <- file.path(INPUT_GLOLAKES_DIR, "BEST_ESTIMATE_basin_total_Mm3.csv")
  
  if (!file.exists(f)) {
    stop(paste("GloLakes best-estimate file not found:", f,
               "\nRun GloLakes_Nechako_Extraction_v2.R first."))
  }
  
  raw <- read.csv(f, stringsAsFactors = FALSE)
  
  reservoir <- raw %>%
    mutate(
      date       = as.Date(date),
      storage    = as.numeric(total_storage_Mm3),
      year       = year(date),
      month      = month(date),
      water_year = if_else(month >= WATER_YEAR_START_MONTH, year, year - 1L),
      n_lakes    = as.numeric(n_lakes_contributing),
      products   = products_used
    ) %>%
    filter(!is.na(storage), storage >= 0, n_lakes >= 1) %>%
    select(date, storage, year, month, water_year, n_lakes, products)
  
  cat(sprintf("  GloLakes storage loaded: %d months (%s to %s)\n",
              nrow(reservoir), min(reservoir$date), max(reservoir$date)))
  cat(sprintf("  Storage range: %.0f – %.0f Mm³ | mean %.0f Mm³\n",
              min(reservoir$storage), max(reservoir$storage), 
              mean(reservoir$storage, na.rm = TRUE)))
  
  # Save diagnostic plot
  if (nrow(reservoir) > 12) {
    p <- ggplot(reservoir, aes(x = date, y = storage)) +
      geom_line(colour = "#1f78b4", linewidth = 0.8) +
      geom_point(colour = "#1f78b4", size = 0.8, alpha = 0.6) +
      labs(title = "Nechako Basin — GloLakes Reservoir Storage",
           subtitle = sprintf("Best-estimate from %d lakes | Q1+Q2 (absolute) only", 
                              mean(reservoir$n_lakes, na.rm = TRUE)),
           x = NULL, y = "Storage (Mm³)",
           caption = "Source: GloLakes v2 (Hou et al. 2024)") +
      theme_bw() + theme(legend.position = "bottom")
    ggsave(file.path(MAIN_OUTPUT_DIR, "diagnostics", "plot_reservoir_storage.png"),
           p, width = 12, height = 5, dpi = 150)
  }
  
  return(reservoir)
}

# ============================================================================
# SECTION 4: ERA5-LAND CLIMATE DATA LOADING
# ============================================================================
ERA5_EPOCH <- as.POSIXct("1900-01-01 00:00:00", tz = "UTC")

load_era5_basin_mean <- function(nc_path, var_name,
                                 lat_min = BASIN_LAT_MIN, lat_max = BASIN_LAT_MAX,
                                 lon_min = BASIN_LON_MIN, lon_max = BASIN_LON_MAX) {
  cat(sprintf("  ERA5: reading '%s' from %s\n", var_name, basename(nc_path)))
  
  nc <- ncdf4::nc_open(nc_path)
  on.exit(ncdf4::nc_close(nc))
  
  lat      <- ncdf4::ncvar_get(nc, "latitude")
  lon      <- ncdf4::ncvar_get(nc, "longitude")
  time_raw <- ncdf4::ncvar_get(nc, "valid_time")
  
  dates <- as.Date(ERA5_EPOCH + as.difftime(time_raw, units = "hours"))
  
  lat_idx <- which(lat >= lat_min & lat <= lat_max)
  lon_idx <- which(lon >= lon_min & lon <= lon_max)
  
  if (length(lat_idx) == 0 | length(lon_idx) == 0)
    stop(sprintf("No ERA5 cells found for '%s' within bounding box", var_name))
  
  start_vec <- c(min(lon_idx), min(lat_idx), 1)
  count_vec <- c(length(lon_idx), length(lat_idx), length(dates))
  data_3d   <- ncdf4::ncvar_get(nc, var_name, start = start_vec, count = count_vec)
  
  data_3d[data_3d >= 9e36] <- NA
  
  lat_w <- cos(lat[lat_idx] * pi / 180)
  
  basin_mean <- vapply(seq_len(length(dates)), function(t) {
    slice <- data_3d[, , t]
    col_means <- apply(slice, 1, function(col_vals) {
      valid <- !is.na(col_vals)
      if (sum(valid) == 0) return(NA_real_)
      weighted.mean(col_vals[valid], lat_w[valid])
    })
    mean(col_means, na.rm = TRUE)
  }, numeric(1))
  
  out <- data.frame(
    date       = dates,
    year       = year(dates),
    month      = month(dates),
    water_year = ifelse(month(dates) >= WATER_YEAR_START_MONTH, 
                        year(dates), year(dates) - 1L),
    stringsAsFactors = FALSE
  )
  out[[var_name]] <- basin_mean
  
  cat(sprintf("    Non-NA values: %d / %d\n", sum(!is.na(basin_mean)), length(basin_mean)))
  return(out)
}

load_era5_all <- function(era5_dir = INPUT_ERA5_DIR) {
  cat("\n[ERA5-Land] Loading basin-averaged monthly climate data...\n")
  
  # ✓ CONFIRMED FILES IN monthly_data_direct/
  f_snowc <- file.path(era5_dir, "snow_cover_monthly.nc")
  f_sd    <- file.path(era5_dir, "snow_depth_water_equivalent_monthly.nc")
  f_tp    <- file.path(era5_dir, "total_precipitation_monthly.nc")
  
  for (f in c(f_snowc, f_sd, f_tp))
    if (!file.exists(f)) stop(paste("ERA5 file missing:", f))
  
  df_snowc <- load_era5_basin_mean(f_snowc, "snowc")
  df_sd    <- load_era5_basin_mean(f_sd, "sd")
  df_tp    <- load_era5_basin_mean(f_tp, "tp")
  
  era5 <- df_snowc %>%
    inner_join(df_sd[, c("date","sd")], by = "date") %>%
    inner_join(df_tp[, c("date","tp")], by = "date") %>%
    mutate(
      swe_mm          = sd * 1000,
      swe_Mm3         = sd * BASIN_AREA_KM2 * 1e6 / 1e6,
      precip_mm       = tp * 1000,
      precip_Mm3      = tp * BASIN_AREA_KM2 * 1e6 / 1e6,
      snow_cover_pct  = snowc * 100
    ) %>%
    arrange(date) %>%
    group_by(water_year) %>%
    mutate(precip_ytd_Mm3 = cumsum(precip_Mm3),
           precip_ytd_mm  = cumsum(precip_mm)) %>%
    ungroup()
  
  cat(sprintf("  ERA5 merged: %d months (%s to %s)\n",
              nrow(era5), min(era5$date), max(era5$date)))
  return(era5)
}

# ============================================================================
# SECTION 5: SNOW PILLOW BIAS CORRECTION
# ============================================================================
load_all_snow_pillows <- function(pillow_dir = INPUT_SNOWPILLOW_DIR) {
  if (!dir.exists(pillow_dir)) {
    cat("  Snow pillow directory not found — skipping bias correction\n")
    return(NULL)
  }
  
  cat("\n[Snow Pillows] Loading BC RFC snow pillow SWE data...\n")
  
  all_files  <- list.files(pillow_dir, full.names = TRUE, pattern = "\\.csv$")
  if (length(all_files) == 0) {
    cat("  No CSV files found — skipping bias correction\n")
    return(NULL)
  }
  
  all_loaded <- list()
  for (f in all_files) {
    raw <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(raw) || nrow(raw) == 0) next
    
    colnames(raw) <- trimws(colnames(raw))
    col_start <- grep("start|Start|Date", colnames(raw), value = TRUE)[1]
    col_swe   <- grep("SW|swe|SWE|Average", colnames(raw), value = TRUE, ignore.case = TRUE)[1]
    
    if (is.na(col_start) || is.na(col_swe)) next
    
    df <- data.frame(
      date_start = tryCatch(lubridate::parse_date_time(raw[[col_start]], 
                                                       orders = c("d/m/Y H:M","m/d/Y H:M","Y-m-d H:M","Y-m-d"), tz = "UTC"),
                            error = function(e) NA),
      swe_mm     = suppressWarnings(as.numeric(raw[[col_swe]])),
      stringsAsFactors = FALSE
    ) %>%
      filter(!is.na(date_start), !is.na(swe_mm), swe_mm > 0) %>%
      mutate(year = year(date_start), month = month(date_start))
    
    if (nrow(df) > 0) {
      monthly <- df %>%
        group_by(year, month) %>%
        summarise(swe_mm_mean = mean(swe_mm, na.rm = TRUE), n_readings = n(), .groups = "drop") %>%
        filter(n_readings >= 5) %>%
        mutate(date = as.Date(paste0(year,"-",sprintf("%02d",month),"-01")))
      
      if (nrow(monthly) > 0) all_loaded[[basename(f)]] <- monthly
    }
  }
  
  if (length(all_loaded) == 0) return(NULL)
  
  combined <- bind_rows(all_loaded)
  basin_mean <- combined %>%
    group_by(year, month, date) %>%
    summarise(swe_pillow_mm = mean(swe_mm_mean, na.rm = TRUE),
              n_stations_avail = n(), .groups = "drop")
  
  cat(sprintf("  Basin-mean pillow SWE: %d months\n", nrow(basin_mean)))
  return(list(by_station = combined, basin_mean = basin_mean))
}

bias_correct_era5_swe <- function(era5_data, pillow_data) {
  if (is.null(pillow_data)) {
    cat("  No snow pillow data — ERA5 SWE used without bias correction\n")
    era5_data$swe_corrected_mm  <- era5_data$swe_mm
    era5_data$swe_corrected_Mm3 <- era5_data$swe_Mm3
    era5_data$swe_bias_factor   <- 1.0
    return(era5_data)
  }
  
  cat("\n[SWE Bias Correction] Estimating ERA5 SWE bias from snow pillows...\n")
  
  paired <- era5_data %>%
    select(year, month, swe_mm, swe_Mm3) %>%
    inner_join(pillow_data$basin_mean %>% select(year, month, swe_pillow_mm),
               by = c("year","month")) %>%
    filter(swe_mm > 1.0, swe_pillow_mm > 1.0)
  
  if (nrow(paired) < MIN_YEARS_FOR_DIST) {
    cat("  Insufficient overlap for bias correction — using raw ERA5\n")
    era5_data$swe_corrected_mm  <- era5_data$swe_mm
    era5_data$swe_corrected_Mm3 <- era5_data$swe_Mm3
    era5_data$swe_bias_factor   <- 1.0
    return(era5_data)
  }
  
  bf <- paired %>%
    group_by(month) %>%
    summarise(bias_factor = median(swe_pillow_mm / swe_mm, na.rm = TRUE),
              n_pairs = n(), .groups = "drop") %>%
    mutate(bias_factor_used = if_else(n_pairs >= MIN_YEARS_FOR_DIST, bias_factor, 1.0))
  
  era5_corrected <- era5_data %>%
    left_join(bf[, c("month","bias_factor_used")], by = "month") %>%
    mutate(
      bias_factor_used  = replace_na(bias_factor_used, 1.0),
      swe_corrected_mm  = swe_mm * bias_factor_used,
      swe_corrected_Mm3 = swe_Mm3 * bias_factor_used,
      swe_bias_factor   = bias_factor_used
    )
  
  write.csv(bf, file.path(MAIN_OUTPUT_DIR, "era5_extracted", "swe_bias_factors.csv"),
            row.names = FALSE)
  return(era5_corrected)
}

# ============================================================================
# SECTION 6: WSC DISCHARGE DATA LOADING
# ============================================================================
load_discharge_data <- function(station_id, data_path = INPUT_WSC_DIR) {
  patterns <- c(
    file.path(data_path, paste0(station_id, "_WSC_DISCHARGE.csv")),
    file.path(data_path, paste0(station_id, "_discharge.csv")),
    file.path(data_path, paste0(station_id, ".csv"))
  )
  
  raw <- NULL
  for (p in patterns) {
    if (file.exists(p)) { raw <- read.csv(p, stringsAsFactors=FALSE); break }
  }
  
  if (is.null(raw)) {
    warning(sprintf("No discharge file for %s", station_id))
    return(NULL)
  }
  
  colnames(raw)  <- tolower(colnames(raw))
  date_col       <- grep("date", colnames(raw), value=TRUE)[1]
  disch_col      <- grep("discharge|flow|q_|value", colnames(raw), value=TRUE)[1]
  
  if (is.na(date_col) || is.na(disch_col)) {
    warning(sprintf("Cannot identify date/discharge columns for %s", station_id))
    return(NULL)
  }
  
  raw$date       <- as.Date(raw[[date_col]])
  raw$discharge  <- suppressWarnings(as.numeric(raw[[disch_col]]))
  raw            <- raw[!is.na(raw$date) & !is.na(raw$discharge), ]
  raw$year       <- year(raw$date)
  raw$month      <- month(raw$date)
  raw$station_id <- station_id
  raw$flow_type  <- "observed"
  
  cat(sprintf("    Observed: %d records (%d-%d)\n",
              nrow(raw), min(raw$year), max(raw$year)))
  return(raw)
}

load_naturalized_flow <- function(station_id, nat_path = INPUT_NATURAL_DIR) {
  fname <- file.path(nat_path, paste0("dat_natural_", station_id, ".csv"))
  if (!file.exists(fname)) {
    warning(sprintf("Naturalized flow file not found: %s", fname))
    return(NULL)
  }
  
  raw <- read.csv(fname, stringsAsFactors=FALSE)
  colnames(raw) <- tolower(trimws(colnames(raw)))
  
  data <- raw %>%
    mutate(date     = as.Date(date),
           discharge = suppressWarnings(as.numeric(value))) %>%
    filter(!is.na(date), tolower(symbol) == "s", !is.na(discharge), discharge >= 0) %>%
    mutate(year       = year(date),
           month      = month(date),
           station_id = station_id,
           flow_type  <- "naturalized") %>%
    select(date, discharge, year, month, station_id, flow_type)
  
  cat(sprintf("    Naturalized: %d records (%d-%d)\n",
              nrow(data), min(data$year), max(data$year)))
  return(data)
}

fill_missing_data <- function(data, max_fill = MAX_FILL_DAYS) {
  if (is.null(data) || nrow(data) == 0) return(NULL)
  
  full_d     <- seq(min(data$date), max(data$date), by="day")
  cdata      <- merge(data.frame(date=full_d), data, by="date", all.x=TRUE)
  cdata      <- cdata[order(cdata$date), ]
  q          <- cdata$discharge
  n          <- length(q)
  runs       <- rle(is.na(q))
  ends       <- cumsum(runs$lengths)
  starts     <- ends - runs$lengths + 1
  
  filled <- 0L
  for (k in seq_along(runs$lengths)) {
    if (!runs$values[k]) next
    rs <- starts[k]; re <- ends[k]; rl <- runs$lengths[k]
    if (rl > max_fill) next
    li <- rs - 1L; ri <- re + 1L
    if (li < 1L || ri > n || is.na(q[li]) || is.na(q[ri])) next
    ctx <- max(1L,li-15L):min(n,ri+15L)
    obs <- !is.na(q[ctx])
    if (sum(obs) < 2) next
    sf       <- splinefun(ctx[obs], q[ctx[obs]], method="monoH.FC")
    q[rs:re] <- pmax(0, sf(rs:re))
    filled   <- filled + rl
  }
  
  cdata$discharge <- q
  cdata$year  <- year(cdata$date)
  cdata$month <- month(cdata$date)
  if ("station_id" %in% colnames(data)) cdata$station_id <- data$station_id[1]
  if ("flow_type"  %in% colnames(data)) cdata$flow_type  <- data$flow_type[1]
  
  cat(sprintf("      Filled: %d values (PCHIP)\n", filled))
  return(cdata)
}

aggregate_to_monthly_volume <- function(daily_data) {
  if (is.null(daily_data) || nrow(daily_data) == 0) return(NULL)
  
  daily_data %>%
    group_by(year, month) %>%
    summarise(n_valid        = sum(!is.na(discharge)),
              discharge_mean = mean(discharge, na.rm=TRUE),
              .groups = "drop") %>%
    filter(n_valid >= MIN_DAYS_PER_MONTH) %>%
    mutate(
      days_in_mo = days_in_month(as.Date(paste0(year,"-",sprintf("%02d",month),"-01"))),
      volume_Mm3 = discharge_mean * days_in_mo * 86400 / 1e6,
      date       = as.Date(paste0(year,"-",sprintf("%02d",month),"-15")),
      water_year = if_else(month >= WATER_YEAR_START_MONTH, year, year - 1L)
    )
}

# ============================================================================
# SECTION 7: STREAMFLOW FORECAST MODELS
# ============================================================================
build_forecast_models <- function(era5_data, monthly_vol_primary) {
  cat("\n[Forecast] Building seasonal streamflow regression models...\n")
  
  apr_sep_total <- monthly_vol_primary %>%
    filter(month %in% MELT_MONTHS) %>%
    group_by(water_year) %>%
    summarise(apr_sep_vol_Mm3 = sum(volume_Mm3, na.rm=TRUE),
              n_melt_complete = sum(!is.na(volume_Mm3)),
              .groups = "drop") %>%
    filter(n_melt_complete == 6)
  
  cum_melt <- monthly_vol_primary %>%
    filter(month %in% MELT_MONTHS) %>%
    arrange(water_year, month) %>%
    group_by(water_year) %>%
    mutate(cum_obs_to_month = cumsum(replace_na(volume_Mm3, 0))) %>%
    ungroup() %>%
    select(water_year, month, cum_obs_to_month)
  
  prev_q <- monthly_vol_primary %>%
    arrange(year, month) %>%
    mutate(prev_volume_Mm3 = lag(volume_Mm3, 1)) %>%
    select(year, month, prev_volume_Mm3)
  
  master <- era5_data %>%
    select(year, month, water_year, swe_corrected_Mm3, precip_ytd_Mm3) %>%
    left_join(prev_q,       by = c("year","month")) %>%
    left_join(apr_sep_total, by = "water_year")
  
  models     <- vector("list", 12)
  names(models) <- month.abb
  forecasts_list <- list()
  
  for (m in 1:12) {
    cat(sprintf("  %s: ", month.abb[m]))
    
    train <- master %>%
      filter(month == m, !is.na(swe_corrected_Mm3), !is.na(precip_ytd_Mm3),
             !is.na(prev_volume_Mm3), !is.na(apr_sep_vol_Mm3))
    
    if (m %in% MELT_MONTHS) {
      train <- train %>%
        left_join(cum_melt %>% filter(month == m), by = c("water_year","month")) %>%
        mutate(target_Mm3 = pmax(0, apr_sep_vol_Mm3 - replace_na(cum_obs_to_month, 0)))
    } else {
      train <- train %>% mutate(target_Mm3 = apr_sep_vol_Mm3)
    }
    
    if (nrow(train) < MIN_YEARS_REGRESSION) {
      cat(sprintf("skipped — only %d years\n", nrow(train)))
      models[[m]] <- NULL
      next
    }
    
    fit <- lm(target_Mm3 ~ swe_corrected_Mm3 + precip_ytd_Mm3 + prev_volume_Mm3, data = train)
    
    loo_pred <- numeric(nrow(train))
    for (j in seq_len(nrow(train))) {
      fit_j      <- lm(target_Mm3 ~ swe_corrected_Mm3 + precip_ytd_Mm3 + prev_volume_Mm3,
                       data = train[-j, ])
      loo_pred[j] <- pmax(0, predict(fit_j, newdata = train[j, ]))
    }
    
    ss_res  <- sum((train$target_Mm3 - loo_pred)^2)
    ss_tot  <- sum((train$target_Mm3 - mean(train$target_Mm3))^2)
    r2_loo  <- 1 - ss_res / ss_tot
    rmse_loo <- sqrt(mean((train$target_Mm3 - loo_pred)^2))
    
    cat(sprintf("n=%d | LOO R²=%.2f | RMSE=%.0f Mm³\n", nrow(train), r2_loo, rmse_loo))
    
    models[[m]] <- fit
    
    all_preds <- master %>%
      filter(month == m, !is.na(swe_corrected_Mm3), !is.na(precip_ytd_Mm3),
             !is.na(prev_volume_Mm3)) %>%
      mutate(
        forecast_Mm3  = pmax(0, predict(fit, newdata = .)),
        loocv_r2      = round(r2_loo, 3),
        loocv_rmse    = round(rmse_loo, 1),
        forecast_type = if_else(m %in% MELT_MONTHS, "remaining_Apr-Sep", "full_Apr-Sep")
      ) %>%
      select(year, month, water_year, forecast_Mm3, loocv_r2, loocv_rmse, forecast_type)
    
    forecasts_list[[m]] <- all_preds
  }
  
  forecast_df <- bind_rows(forecasts_list)
  write.csv(forecast_df %>%
              group_by(month) %>%
              summarise(loocv_r2=first(loocv_r2), loocv_rmse=first(loocv_rmse),
                        n_years=n(), .groups="drop"),
            file.path(MAIN_OUTPUT_DIR, "forecast_models", "regression_diagnostics.csv"),
            row.names = FALSE)
  
  saveRDS(models, file.path(MAIN_OUTPUT_DIR, "forecast_models", "regression_models.rds"))
  return(list(models = models, forecasts = forecast_df))
}

# ============================================================================
# SECTION 8: NONEXCEEDANCE PROBABILITY & COMPONENT INDEX
# ============================================================================
compute_weibull_prob <- function(hist_vals, x_new) {
  v <- sort(hist_vals[!is.na(hist_vals)])
  n <- length(v)
  if (n == 0 || is.na(x_new)) return(NA_real_)
  P <- seq_len(n) / (n + 1) * 100
  pmax(0, pmin(100, approx(v, P, xout=x_new, rule=2, method="linear")$y))
}

fit_parametric_cdf <- function(hist_vals) {
  v <- hist_vals[!is.na(hist_vals) & hist_vals > 0]
  if (length(v) < MIN_YEARS_FOR_DIST) return(NULL)
  
  lmom <- tryCatch(lmomco::lmoms(v, nmom=4), error=function(e) NULL)
  if (is.null(lmom)) return(NULL)
  
  fits <- lapply(c("gamma", "lnorm", "norm"), function(d) {
    tryCatch({
      if (d == "gamma") {
        par <- lmomco::pargam(lmom); cdf <- function(x) lmomco::cdfgam(x, par)
      } else if (d == "lnorm") {
        par <- lmomco::parlnorm(lmom); cdf <- function(x) lmomco::cdflnorm(x, par)
      } else {
        par <- lmomco::parnor(lmom); cdf <- function(x) lmomco::cdfnor(x, par)
      }
      list(dist=d, cdf=cdf)
    }, error=function(e) list(dist=d, cdf=NULL))
  })
  
  best <- fits[[which.min(sapply(fits, function(x) if(is.null(x$cdf)) Inf else 0))]]
  if (is.null(best$cdf)) return(NULL)
  list(best_dist=best$dist, cdf_fn=best$cdf)
}

compute_component_index <- function(df, value_col, label) {
  df[[paste0("P_",   label)]] <- NA_real_
  df[[paste0("idx_", label)]] <- NA_real_
  df[[paste0("dist_",label)]] <- NA_character_
  
  for (m in 1:12) {
    row_m    <- which(df$month == m)
    hist_all <- df[[value_col]][row_m]
    if (sum(!is.na(hist_all)) < MIN_YEARS_FOR_DIST) next
    
    param_fit <- fit_parametric_cdf(hist_all)
    
    for (r in row_m) {
      v_r <- df[[value_col]][r]
      if (is.na(v_r)) next
      
      P_emp <- compute_weibull_prob(hist_all, v_r)
      
      if (!is.null(param_fit)) {
        P_use  <- tryCatch(pmax(0, pmin(100, param_fit$cdf_fn(v_r)*100)),
                           error=function(e) P_emp)
        dist_used <- param_fit$best_dist
      } else {
        P_use     <- P_emp
        dist_used <- "empirical"
      }
      
      df[[paste0("P_",   label)]][r] <- P_use
      df[[paste0("idx_", label)]][r] <- (P_use - 50) / 12
      df[[paste0("dist_",label)]][r] <- dist_used
    }
  }
  return(df)
}

# ============================================================================
# SECTION 9: FULL SWSI COMPUTATION (Garen 1993)
# ============================================================================
compute_swsi <- function(era5_data, monthly_vol, forecast_result, 
                         reservoir_data, station_id) {
  cat(sprintf("\n[SWSI] Computing for station %s...\n", station_id))
  
  combined <- monthly_vol %>%
    select(year, month, water_year, date, volume_Mm3) %>%
    left_join(forecast_result$forecasts %>% 
                select(year, month, forecast_Mm3, forecast_type),
              by = c("year","month")) %>%
    left_join(era5_data %>% select(year, month, swe_corrected_Mm3, 
                                   swe_corrected_mm, precip_ytd_Mm3, 
                                   precip_ytd_mm, snow_cover_pct),
              by = c("year","month"))
  
  cat("  Component 1: observed streamflow...\n")
  combined <- compute_component_index(combined, "volume_Mm3", "streamflow")
  
  cat("  Component 2: forecast streamflow...\n")
  combined <- compute_component_index(combined, "forecast_Mm3", "forecast")
  
  cat("  Component 3: SWE...\n")
  combined <- compute_component_index(combined, "swe_corrected_Mm3", "swe")
  
  cat("  Component 4: precipitation YTD...\n")
  combined <- compute_component_index(combined, "precip_ytd_Mm3", "precip")
  
  if (RESERVOIR_DATA_AVAILABLE && !is.null(reservoir_data)) {
    cat("  Component 5: reservoir storage...\n")
    combined <- combined %>%
      left_join(reservoir_data %>% select(year, month, storage),
                by = c("year","month"))
    combined <- compute_component_index(combined, "storage", "reservoir")
    
    combined <- combined %>%
      mutate(total_water_Mm3 = replace_na(storage, 0) + replace_na(forecast_Mm3, 0))
    combined <- compute_component_index(combined, "total_water_Mm3", "total")
    
    combined$swsi      <- combined$idx_total
    combined$swsi_type <- "FULL"
  } else {
    combined$swsi      <- combined$idx_forecast
    combined$swsi_type <- "PARTIAL_no_reservoir"
    cat("  [NOTE] Partial SWSI — reservoir storage not available\n")
  }
  
  combined$drought_category <- cut(
    combined$swsi,
    breaks = c(-Inf, -4.0, -3.0, -2.0, -1.0, 0.0, 2.0, Inf),
    labels = c("Extreme drought","Severe drought","Moderate drought",
               "Drought watch","Below normal","Near normal","Abundant supply"),
    right = TRUE
  )
  combined$station_id <- station_id
  
  cat(sprintf("  SWSI computed for %d months\n", sum(!is.na(combined$swsi))))
  return(combined)
}

# ============================================================================
# SECTION 10: DROUGHT EVENT IDENTIFICATION
# ============================================================================
identify_swsi_droughts <- function(swsi_df, station_id=NULL) {
  if (is.null(swsi_df)) return(NULL)
  
  data <- swsi_df[order(swsi_df$date), ]
  data$in_drought <- !is.na(data$swsi) & data$swsi <= TRIGGER_WATCH
  data$event_id <- with(rle(data$in_drought), rep(seq_along(lengths), lengths))
  
  event_info <- data %>%
    group_by(event_id) %>%
    filter(in_drought == TRUE) %>%
    summarise(
      start_date       = min(date),
      end_date         = max(date),
      duration_months  = n(),
      start_year       = first(year),
      water_year_start = first(water_year),
      swsi_min         = min(swsi, na.rm=TRUE),
      swsi_mean        = mean(swsi, na.rm=TRUE),
      idx_swe_min      = min(idx_swe, na.rm=TRUE),
      idx_precip_min   = min(idx_precip, na.rm=TRUE),
      idx_forecast_min = min(idx_forecast, na.rm=TRUE),
      .groups = "drop"
    ) %>%
    filter(duration_months >= 1) %>%
    mutate(
      event_id_new = row_number(),
      peak_category = case_when(
        swsi_min <= -4.0 ~ "Extreme drought",
        swsi_min <= -3.0 ~ "Severe drought",
        swsi_min <= -2.0 ~ "Moderate drought",
        TRUE             ~ "Drought watch"
      ),
      primary_driver = case_when(
        idx_swe_min <= idx_forecast_min & idx_swe_min <= idx_precip_min ~ "SWE deficit",
        idx_precip_min <= idx_forecast_min ~ "Precipitation deficit",
        TRUE ~ "Streamflow forecast deficit"
      )
    )
  
  if (!is.null(station_id)) event_info$station_id <- station_id
  
  n_emergency <- sum(event_info$swsi_min <= TRIGGER_EMERGENCY, na.rm=TRUE)
  cat(sprintf("  Events: %d | Emergency: %d\n", nrow(event_info), n_emergency))
  
  return(list(monthly_data = data, drought_events = event_info,
              n_events = nrow(event_info), n_emergency = n_emergency))
}

# ============================================================================
# SECTION 11: MAIN PROCESSING
# ============================================================================
cat("\n", strrep("=",70), "\n",
    "SWSI COMPUTATION — Nechako Basin (Garen 1993)\n",
    "Data: GloLakes v2 (1984+) + ERA5-Land + WSC\n",
    strrep("=",70), "\n", sep="")

# ── A: Load Reservoir Storage (GloLakes) ─────────────────────────────────────
reservoir_data <- load_reservoir_storage()

# ── B: Load ERA5-Land Climate ────────────────────────────────────────────────
era5_raw  <- load_era5_all(INPUT_ERA5_DIR)
pillow_data <- load_all_snow_pillows(INPUT_SNOWPILLOW_DIR)
era5_data   <- bias_correct_era5_swe(era5_raw, pillow_data)

write.csv(era5_data, file.path(MAIN_OUTPUT_DIR, "era5_extracted", 
                               "era5_basin_monthly_corrected.csv"),
          row.names = FALSE)

# ── C: Load Streamflow Data ──────────────────────────────────────────────────
stations <- data.frame(
  StationID   = c("08JC001", "08JC002", "08JB002", "08JA015"),
  Name        = c("Nechako R at Vanderhoof", "Nechako R above Nautley R",
                  "Nechako R at Isle Pierre", "Stuart Lake outlet"),
  Has_NatFlow = c(TRUE, TRUE, FALSE, FALSE),
  Use_in_SWSI = c(TRUE, TRUE, TRUE, TRUE),
  stringsAsFactors = FALSE
)

all_monthly_vol <- list()
forecast_result <- NULL

for (i in seq_len(nrow(stations))) {
  if (!stations$Use_in_SWSI[i]) next
  
  sid <- stations$StationID[i]
  cat(sprintf("\n%s — %s\n", sid, stations$Name[i]))
  
  if (stations$Has_NatFlow[i]) {
    daily_data <- load_naturalized_flow(sid, INPUT_NATURAL_DIR)
  } else {
    daily_data <- load_discharge_data(sid, INPUT_WSC_DIR)
  }
  
  if (is.null(daily_data)) next
  
  data_filled <- fill_missing_data(daily_data)
  monthly_vol <- aggregate_to_monthly_volume(data_filled)
  all_monthly_vol[[sid]] <- monthly_vol
  
  if (sid == PRIMARY_STATION) {
    forecast_result <- build_forecast_models(era5_data, monthly_vol)
  }
}

if (is.null(forecast_result))
  stop(sprintf("Primary station %s not processed", PRIMARY_STATION))

# ── D: Compute SWSI ──────────────────────────────────────────────────────────
all_swsi     <- list()
all_droughts <- list()

for (sid in names(all_monthly_vol)) {
  swsi_df <- compute_swsi(era5_data = era5_data,
                          monthly_vol = all_monthly_vol[[sid]],
                          forecast_result = forecast_result,
                          reservoir_data = reservoir_data,
                          station_id = sid)
  all_swsi[[sid]] <- swsi_df
  
  droughts <- identify_swsi_droughts(swsi_df, station_id=sid)
  all_droughts[[sid]] <- droughts
  
  write.csv(swsi_df, file.path(MAIN_OUTPUT_DIR, "swsi_timeseries", 
                               paste0(sid,"_swsi.csv")), row.names = FALSE)
  
  if (!is.null(droughts$drought_events) && nrow(droughts$drought_events) > 0) {
    write.csv(droughts$drought_events, 
              file.path(MAIN_OUTPUT_DIR, "drought_events", paste0(sid,"_droughts.csv")),
              row.names = FALSE)
  }
}

# ── E: Save Master Results ───────────────────────────────────────────────────
all_results <- list(
  stations        = stations,
  era5_data       = era5_data,
  reservoir_data  = reservoir_data,
  monthly_volumes = all_monthly_vol,
  forecast_result = forecast_result,
  swsi            = all_swsi,
  droughts        = all_droughts,
  metadata = list(
    index_type       = "SWSI — Surface-Water Supply Index",
    reference        = "Garen (1993) J.WRPM 119(4):437-454",
    basin            = "Nechako River Basin, BC, Canada",
    swsi_formula     = "SWSI = (P - 50) / 12",
    swsi_type        = "FULL (reservoir + forecast streamflow)",
    time_span        = sprintf("%s to %s", 
                               min(reservoir_data$date), 
                               max(reservoir_data$date)),
    streamflow_source = "WSC daily discharge; naturalized for 08JC001/08JC002",
    reservoir_source  = "GloLakes v2 (Hou et al. 2024)",
    era5_variables    = c("snowc", "sd (SWE)", "tp (precip)"),
    processed_date    = Sys.Date(),
    script_version    = "3.0_complete_GloLakes_corrected_paths"
  )
)

saveRDS(all_results, file.path(MAIN_OUTPUT_DIR, "all_results_swsi.rds"))

# ── F: Summary Report ────────────────────────────────────────────────────────
cat(sprintf("\n%s\n", strrep("=",70)))
cat("SWSI ANALYSIS SUMMARY — Nechako Basin\n")
cat(sprintf("%s\n", strrep("=",70)))
cat(sprintf("SWSI type:    FULL (reservoir + forecast streamflow)\n"))
cat(sprintf("Time span:    %s to %s (%.1f years)\n",
            min(reservoir_data$date), max(reservoir_data$date),
            as.numeric(difftime(max(reservoir_data$date), 
                                min(reservoir_data$date), units="days"))/365.25))
cat(sprintf("Reservoir:    GloLakes v2 (mean %.0f Mm³, %d lakes)\n",
            mean(reservoir_data$storage, na.rm=TRUE),
            mean(reservoir_data$n_lakes, na.rm=TRUE)))
cat(sprintf("ERA5 period:  %s to %s\n", min(era5_data$date), max(era5_data$date)))
cat(sprintf("Stations:     %d processed\n", length(all_swsi)))
cat(sprintf("\nOutputs → %s/\n", MAIN_OUTPUT_DIR))
cat("  swsi_timeseries/    Monthly SWSI + component indexes per station\n")
cat("  drought_events/     Drought event tables\n")
cat("  forecast_models/    Regression diagnostics\n")
cat("  all_results_swsi.rds  Master combined output\n")

cat(sprintf("\n%s\n", strrep("=",70)))
cat("Garen (1993) SWSI Categories:\n")
cat("  >= +2.0       Abundant supply\n")
cat("  -2.0 to +2.0  Near normal\n")
cat("  -3.0 to -2.0  Moderate drought\n")
cat("  -4.0 to -3.0  Severe drought\n")
cat("  < -4.0        Extreme drought\n")
cat(sprintf("Triggers:  Watch ≤ %.1f  |  Emergency ≤ %.1f\n",
            TRIGGER_WATCH, TRIGGER_EMERGENCY))
cat(sprintf("%s\n", strrep("=",70)))