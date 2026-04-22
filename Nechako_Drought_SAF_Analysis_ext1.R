#============================================================================
# DROUGHT SEVERITY-AREA-FREQUENCY ANALYSIS FOR NECHAKO BASIN
# Using PNPI & PNWBI with DURATION CLASSIFICATION + DUAL RETURN PERIOD METHODS
#============================================================================
# This script performs regional drought analysis using:
# 1. Percent of Normal Precipitation Index (PNPI)
# 2. Percent of Normal Water Balance Index (PNWBI)
# 3. Duration classification (Sheffield & Wood 2008):
#    D4–6   : Short-term  (4–6 consecutive months below threshold)
#    D7–12  : Medium-term (7–12 consecutive months below threshold)
#    D12+   : Long-term   (>12 consecutive months below threshold)
# 4. Per-class empirical return periods (events / record length)
# 5. TWO multivariate return period methods applied per duration class:
#    a) Conditional copula approach (Amirataee et al. 2018)
#    b) Kendall distribution function approach (Genest et al. 2009)
#
# NEW OUTPUTS (added):
#   *_duration_classified_events.csv  — all events with duration class
#   *_return_periods_by_class.csv     — empirical T per D4-6 / D7-12 / D12+
#   *_SAF_curves_conditional_<class>.csv  — SAF per duration class (Method A)
#   *_SAF_curves_kendall_<class>.csv      — SAF per duration class (Method B)
#   *_return_period_summary_plot.pdf      — bar chart of T by duration class
#
# Key References:
# - PNPI methodology: WMO (2016); Heim (2002)
# - Water balance: Palmer (1965); Thornthwaite (1948)
# - Duration classes: Sheffield & Wood (2008), Clim. Dyn. 31:79-105
# - Conditional approach: Amirataee et al. (2018); Shiau (2006)
# - Kendall approach: Genest et al. (2009); Salvadori & De Michele (2004)
#============================================================================

--- 0. SETUP AND LOAD LIBRARIES ---
  rm(list = ls())
gc()

required_packages <- c("terra", "copula", "fitdistrplus", "MASS",, "gamlss", "gamlss.dist", "lmomco", "lmomRFA")
"ggplot2", "gridExtra", "viridis", "moments",
"dplyr", "openxlsx", "Kendall")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

setwd("D:/Nechako_Drought/monthly_data_direct")

# ── Duration class boundaries (Sheffield & Wood 2008) ─────────────────
D_SHORT_MIN <- 4L;  D_SHORT_MAX <- 6L    # D4-6  : short-term
D_MED_MIN   <- 7L;  D_MED_MAX  <- 12L   # D7-12 : medium-term
D_LONG_MIN  <- 13L                        # D12+  : long-term
DROUGHT_THRESHOLD_PCT <- 75               # <75% = moderate drought (WMO 2016)

#' Classify a duration (integer months) into a Sheffield duration class
classify_duration_class <- function(dur) {
  dplyr::case_when(
    dur >= D_SHORT_MIN & dur <= D_SHORT_MAX ~ "D4-6 (Short-term)",
    dur >= D_MED_MIN   & dur <= D_MED_MAX   ~ "D7-12 (Medium-term)",
    dur >= D_LONG_MIN                        ~ "D12+ (Long-term)",
    TRUE                                     ~ "D1-3 (Sub-threshold)"
  )
}

# Create output directory structure
dir.create("drought_analysis", showWarnings = FALSE)
dir.create("drought_analysis/figures", showWarnings = FALSE)
dir.create("drought_analysis/results", showWarnings = FALSE)
dir.create("drought_analysis/PNPI_analysis", showWarnings = FALSE)
dir.create("drought_analysis/PNWBI_analysis", showWarnings = FALSE)

# Initialize log file
LOG_FILE <- "drought_analysis/SAF_analysis_log.txt"
cat("======================================================\n", file = LOG_FILE)
cat("DROUGHT S-A-F ANALYSIS USING PNPI & PNWBI - NECHAKO BASIN\n", file = LOG_FILE, append = TRUE)
cat(paste("Analysis started:", Sys.time(), "\n"), file = LOG_FILE, append = TRUE)
cat("======================================================\n\n", file = LOG_FILE, append = TRUE)

log_event <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(paste0(timestamp, " | ", msg, "\n"), file = LOG_FILE, append = TRUE)
  message(paste0(timestamp, " | ", msg))
}

#============================================================================
# STEP 1: LOAD DATA
#============================================================================
log_event("="*70)
log_event("STEP 1: LOADING PRECIPITATION AND PET DATA")
log_event("="*70)

precip <- rast("total_precipitation_monthly.nc")
log_event(sprintf("  Precipitation layers: %d", nlyr(precip)))

pet <- rast("ERA5Land_Nechako_PET_monthly.nc")
log_event(sprintf("  PET layers: %d", nlyr(pet)))

if (nlyr(precip) != nlyr(pet)) {
  stop("ERROR: Precipitation and PET have different temporal lengths!")
}
n_months <- nlyr(precip)
n_cells <- ncell(precip)

dates <- seq.Date(from = as.Date("1950-01-01"), by = "month", length.out = n_months)
month_indices <- as.numeric(format(dates, "%m"))
year_indices <- as.numeric(format(dates, "%Y"))

cell_area_km2 <- cellSize(precip[[1]], unit = "km")
total_area_km2 <- global(cell_area_km2, "sum", na.rm = TRUE)[1,1]
log_event(sprintf("  Total basin area: %.2f km²", total_area_km2))

#============================================================================
# STEP 2: PREPARE PET IN MONTHLY UNITS
#============================================================================
log_event(paste(rep("=", 70), collapse = ""))
log_event("STEP 2: CONVERTING PET TO MONTHLY TOTALS")
log_event("="*70)

days_per_month <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
days_in_month <- days_per_month[month_indices]
is_leap <- (year_indices %% 4 == 0 & year_indices %% 100 != 0) | (year_indices %% 400 == 0)
days_in_month[month_indices == 2 & is_leap] <- 29

pet_monthly <- pet * days_in_month
pet_stats <- global(pet_monthly, c("min", "max", "mean"), na.rm = TRUE)
log_event(sprintf("  Monthly PET range: %.2f to %.2f mm (mean: %.2f)",
                  min(pet_stats$min), max(pet_stats$max), mean(pet_stats$mean)))

#============================================================================
# STEP 3: CALCULATE WATER BALANCE (P - PET)
#============================================================================
log_event(paste(rep("=", 70), collapse = ""))
log_event("STEP 3: CALCULATING WATER BALANCE (P - PET)")
log_event("="*70)

wb <- precip - pet_monthly
wb <- ifel(wb < 0, 0, wb)  # Apply constraint: negative WB → 0
wb_stats <- global(wb, c("min", "max", "mean"), na.rm = TRUE)
log_event(sprintf("  Water balance range: %.2f to %.2f mm (mean: %.2f)",
                  min(wb_stats$min), max(wb_stats$max), mean(wb_stats$mean)))
writeCDF(wb, "drought_analysis/water_balance_monthly.nc", overwrite = TRUE)

#============================================================================
# FUNCTION: CALCULATE PERCENT OF NORMAL BY CALENDAR MONTH
#============================================================================
percent_of_normal_by_month <- function(data_rast, dates, name = "Index") {
  # Calculate Percent of Normal Index by calendar month climatology
  # PNPI/PNWBI = (X_t / mu_m) * 100
  # Reference: WMO (2016); Heim (2002)
  log_event(sprintf("  Calculating Percent of Normal %s by calendar month...", name))
  
  month_idx <- as.numeric(format(dates, "%m"))
  n_layers <- nlyr(data_rast)
  pn_rast <- rast(data_rast)
  names(pn_rast) <- paste0(name, "_", format(dates, "%Y-%m"))
  
  # Calculate climatology for each calendar month
  for (m in 1:12) {
    log_event(sprintf("    Calendar month %02d (%s)...", m, month.abb[m]))
    idx_m <- which(month_idx == m)
    data_m <- data_rast[[idx_m]]
    
    # Calculate long-term mean for this calendar month
    mu_m <- mean(data_m, na.rm = TRUE)
    mu_m <- ifel(mu_m == 0, 0.001, mu_m)  # Avoid division by zero
    
    # Calculate percent of normal for each layer
    for (i in idx_m) {
      pn_rast[[i]] <- (data_rast[[i]] / mu_m) * 100
    }
    
    mu_mean <- global(mu_m, "mean", na.rm = TRUE)[1,1]
    log_event(sprintf("      Long-term mean = %.2f, Expected PN ≈ 100%%", mu_mean))
  }
  return(pn_rast)
}

#============================================================================
# STEP 4A: CALCULATE PERCENT OF NORMAL PRECIPITATION INDEX (PNPI)
#============================================================================
log_event(paste(rep("=", 70), collapse = ""))
log_event("STEP 4A: CALCULATING PERCENT OF NORMAL PRECIPITATION INDEX (PNPI)")
log_event("="*70)
log_event("  Reference: WMO (2016) - Percent of Normal as drought indicator")
log_event("  Interpretation: 100% = normal, <75% = moderate drought")

pnpi <- percent_of_normal_by_month(precip, dates, name = "PNPI")
pnpi_stats <- global(pnpi, c("min", "max", "mean"), na.rm = TRUE)
log_event(sprintf("  PNPI range: %.1f%% to %.1f%% (mean: %.1f%%)",
                  min(pnpi_stats$min), max(pnpi_stats$max), mean(pnpi_stats$mean)))
writeCDF(pnpi, "drought_analysis/PNPI_monthly.nc", overwrite = TRUE)
log_event("  PNPI saved to: PNPI_monthly.nc")

#============================================================================
# STEP 4B: CALCULATE PERCENT OF NORMAL WATER BALANCE INDEX (PNWBI)
#============================================================================
log_event(paste(rep("=", 70), collapse = ""))
log_event("STEP 4B: CALCULATING PERCENT OF NORMAL WATER BALANCE INDEX (PNWBI)")
log_event("="*70)
log_event("  Reference: Palmer (1965) water balance approach adapted to percentage format")

pnwbi <- percent_of_normal_by_month(wb, dates, name = "PNWBI")
pnwbi_stats <- global(pnwbi, c("min", "max", "mean"), na.rm = TRUE)
log_event(sprintf("  PNWBI range: %.1f%% to %.1f%% (mean: %.1f%%)",
                  min(pnwbi_stats$min), max(pnwbi_stats$max), mean(pnwbi_stats$mean)))
writeCDF(pnwbi, "drought_analysis/PNWBI_monthly.nc", overwrite = TRUE)
log_event("  PNWBI saved to: PNWBI_monthly.nc")

#============================================================================
# FUNCTION: IDENTIFY INDIVIDUAL DROUGHT EVENTS WITH DURATION CLASSIFICATION
# (Sheffield & Wood 2008 approach applied to PNPI/PNWBI)
#============================================================================
#' Scan the basin-averaged PN index time series and return a data.frame of
#' individual drought events including their duration class (D4-6, D7-12, D12+).
#'
#' Entry condition : index < DROUGHT_THRESHOLD_PCT  for at least 1 month
#' Exit condition  : index >= DROUGHT_THRESHOLD_PCT
#'
#' For each event the function records:
#'   start_date, end_date, duration_months, duration_class,
#'   mean_severity (mean % deficit below 100), max_severity (max deficit),
#'   mean_area_pct (mean % basin in drought during event)
identify_duration_classified_events <- function(drought_chars, index_name,
                                                threshold_pct = DROUGHT_THRESHOLD_PCT) {
  log_event(sprintf("  Identifying duration-classified events for %s...", index_name))
  
  dc   <- drought_chars[order(drought_chars$date), ]
  vals <- dc$severity     # mean deficit (0 when no drought)
  area <- dc$area_pct
  indr <- dc$severity > 0 & dc$area_pct > 0   # in-drought months
  
  events <- list()
  in_d   <- FALSE
  s_idx  <- NA_integer_
  
  for (i in seq_len(nrow(dc))) {
    if (!in_d && indr[i]) {
      in_d  <- TRUE
      s_idx <- i
    } else if (in_d && !indr[i]) {
      e_idx <- i - 1L
      dur   <- e_idx - s_idx + 1L
      events[[length(events) + 1L]] <- data.frame(
        index           = index_name,
        start_date      = dc$date[s_idx],
        end_date        = dc$date[e_idx],
        start_year      = dc$year[s_idx],
        end_year        = dc$year[e_idx],
        duration_months = dur,
        duration_class  = classify_duration_class(dur),
        mean_severity   = mean(vals[s_idx:e_idx], na.rm = TRUE),
        max_severity    = max(vals[s_idx:e_idx],  na.rm = TRUE),
        mean_area_pct   = mean(area[s_idx:e_idx], na.rm = TRUE),
        max_area_pct    = max(area[s_idx:e_idx],  na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      in_d  <- FALSE
      s_idx <- NA_integer_
    }
  }
  # Close open event at end of record
  if (in_d && !is.na(s_idx)) {
    e_idx <- nrow(dc)
    dur   <- e_idx - s_idx + 1L
    events[[length(events) + 1L]] <- data.frame(
      index           = index_name,
      start_date      = dc$date[s_idx],
      end_date        = dc$date[e_idx],
      start_year      = dc$year[s_idx],
      end_year        = dc$year[e_idx],
      duration_months = dur,
      duration_class  = classify_duration_class(dur),
      mean_severity   = mean(vals[s_idx:e_idx], na.rm = TRUE),
      max_severity    = max(vals[s_idx:e_idx],  na.rm = TRUE),
      mean_area_pct   = mean(area[s_idx:e_idx], na.rm = TRUE),
      max_area_pct    = max(area[s_idx:e_idx],  na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  
  if (length(events) == 0) {
    log_event(sprintf("  ⚠ No events detected for %s", index_name))
    return(data.frame())
  }
  out <- do.call(rbind, events)
  log_event(sprintf("  ✓ %d total events identified", nrow(out)))
  for (cl in c("D4-6 (Short-term)", "D7-12 (Medium-term)", "D12+ (Long-term)")) {
    n <- sum(out$duration_class == cl)
    log_event(sprintf("    %s: %d events", cl, n))
  }
  out
}

#============================================================================
# FUNCTION: EMPIRICAL RETURN PERIODS BY DURATION CLASS
#============================================================================
#' For each duration class compute:
#'   n_events, freq_per_year, return_period_years,
#'   mean/max severity, mean/max area, mean duration
compute_return_periods_by_class <- function(events_df, index_name,
                                            record_years = NULL) {
  log_event(sprintf("  Computing return periods by duration class for %s...", index_name))
  
  if (is.null(record_years)) {
    record_years <- as.integer(
      max(events_df$end_year) - min(events_df$start_year) + 1L
    )
  }
  
  classes <- c("D4-6 (Short-term)", "D7-12 (Medium-term)", "D12+ (Long-term)")
  out <- lapply(classes, function(cl) {
    sub <- events_df[events_df$duration_class == cl, ]
    n   <- nrow(sub)
    if (n == 0) {
      return(data.frame(
        index                = index_name,
        duration_class       = cl,
        n_events             = 0L,
        record_years         = record_years,
        freq_per_year        = 0,
        return_period_years  = Inf,
        mean_duration_months = NA_real_,
        mean_severity        = NA_real_,
        max_severity         = NA_real_,
        mean_area_pct        = NA_real_,
        max_area_pct         = NA_real_,
        stringsAsFactors     = FALSE
      ))
    }
    data.frame(
      index                = index_name,
      duration_class       = cl,
      n_events             = n,
      record_years         = record_years,
      freq_per_year        = n / record_years,
      return_period_years  = record_years / n,
      mean_duration_months = mean(sub$duration_months, na.rm = TRUE),
      mean_severity        = mean(sub$mean_severity,   na.rm = TRUE),
      max_severity         = max(sub$max_severity,     na.rm = TRUE),
      mean_area_pct        = mean(sub$mean_area_pct,   na.rm = TRUE),
      max_area_pct         = max(sub$max_area_pct,     na.rm = TRUE),
      stringsAsFactors     = FALSE
    )
  })
  rp <- do.call(rbind, out)
  log_event(sprintf("  ✓ Return period table computed (%d rows)", nrow(rp)))
  rp
}

#============================================================================
# FUNCTION: PLOT RETURN PERIODS BY DURATION CLASS
#============================================================================
plot_return_periods_by_class <- function(rp_df, index_name, output_dir) {
  rp_plot <- rp_df[is.finite(rp_df$return_period_years) &
                     rp_df$duration_class != "D1-3 (Sub-threshold)", ]
  
  p <- ggplot(rp_plot, aes(x = duration_class, y = return_period_years,
                           fill = duration_class)) +
    geom_col(width = 0.6, colour = "grey30") +
    geom_text(aes(label = sprintf("n=%d\nT=%.1f yr", n_events, return_period_years)),
              vjust = -0.4, size = 3.5) +
    scale_fill_manual(
      values = c("D4-6 (Short-term)"   = "#fdae61",
                 "D7-12 (Medium-term)" = "#f46d43",
                 "D12+ (Long-term)"    = "#a50026"),
      guide = "none") +
    labs(title    = sprintf("%s — Empirical Return Periods by Duration Class", index_name),
         subtitle = sprintf("Record: %d years (1950–%d); T = record_years / n_events",
                            rp_plot$record_years[1],
                            1950L + rp_plot$record_years[1] - 1L),
         x = "Duration Class", y = "Return Period (Years)") +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  pdf(file.path(output_dir,
                sprintf("%s_return_period_by_class.pdf", index_name)),
      width = 8, height = 6)
  print(p)
  dev.off()
  log_event(sprintf("    Return-period bar chart saved for %s", index_name))
}

#============================================================================
# FUNCTION: EXTRACT DROUGHT CHARACTERISTICS FOR PERCENT INDICES
#============================================================================
extract_drought_characteristics_percent <- function(index_rast, cell_area,
                                                    threshold_pct = DROUGHT_THRESHOLD_PCT) {
  # Extract regional drought severity and area using percent-based indices
  # Drought threshold: <75% = moderate drought (WMO 2016; Heim 2002)
  # Severity definition: Mean deficit from normal (100 - index) in drought areas
  log_event(sprintf("  Extracting drought characteristics (threshold = %.0f%%)...", 
                    threshold_pct))
  
  n_time <- nlyr(index_rast)
  total_area <- global(cell_area, "sum", na.rm = TRUE)[1,1]
  
  severity <- numeric(n_time)    # Mean deficit (%) from normal conditions
  area_pct <- numeric(n_time)    # % of basin below threshold
  
  pb <- txtProgressBar(min = 1, max = n_time, style = 3)
  for (t in 1:n_time) {
    index_t <- index_rast[[t]]
    drought_mask <- index_t < threshold_pct
    
    # AREA UNDER DROUGHT
    drought_area <- mask(cell_area, drought_mask, maskvalues = 0)
    area_t <- global(drought_area, "sum", na.rm = TRUE)[1,1]
    area_pct[t] <- (area_t / total_area) * 100
    
    # DROUGHT SEVERITY: Mean deficit from normal (100 - index value)
    deficit_values <- mask(100 - index_t, drought_mask, maskvalues = 0)
    if (area_pct[t] > 0) {
      severity[t] <- global(deficit_values, "mean", na.rm = TRUE)[1,1]
    } else {
      severity[t] <- 0
    }
    
    setTxtProgressBar(pb, t)
  }
  close(pb)
  
  drought_chars <- data.frame(
    date = dates,
    year = year_indices,
    month = month_indices,
    severity = severity,
    area_pct = area_pct
  )
  
  log_event(sprintf("  Severity range: %.1f%% to %.1f%% deficit", 
                    min(severity), max(severity)))
  log_event(sprintf("  Area range: %.1f%% to %.1f%% of basin", 
                    min(area_pct), max(area_pct)))
  return(drought_chars)
}

#============================================================================
# STEP 5: EXTRACT DROUGHT CHARACTERISTICS AND CLASSIFY EVENTS (PNPI & PNWBI)
#============================================================================
log_event(paste(rep("=", 70), collapse = ""))
log_event("STEP 5: EXTRACTING DROUGHT CHARACTERISTICS AND CLASSIFYING EVENTS")
log_event(paste(rep("=", 70), collapse = ""))
log_event("  Using drought threshold: 75% (moderate drought per WMO classification)")

pnpi_drought <- extract_drought_characteristics_percent(
  pnpi, cell_area_km2, threshold_pct = DROUGHT_THRESHOLD_PCT
)
write.csv(pnpi_drought, "drought_analysis/PNPI_analysis/PNPI_drought_characteristics.csv",
          row.names = FALSE)

pnwbi_drought <- extract_drought_characteristics_percent(
  pnwbi, cell_area_km2, threshold_pct = DROUGHT_THRESHOLD_PCT
)
write.csv(pnwbi_drought, "drought_analysis/PNWBI_analysis/PNWBI_drought_characteristics.csv",
          row.names = FALSE)

# ── NEW: Duration-classified event catalogs (Sheffield & Wood 2008) ───
pnpi_events  <- identify_duration_classified_events(pnpi_drought,  "PNPI")
pnwbi_events <- identify_duration_classified_events(pnwbi_drought, "PNWBI")

write.csv(pnpi_events,
          "drought_analysis/PNPI_analysis/PNPI_duration_classified_events.csv",
          row.names = FALSE)
write.csv(pnwbi_events,
          "drought_analysis/PNWBI_analysis/PNWBI_duration_classified_events.csv",
          row.names = FALSE)

# ── NEW: Empirical return periods by duration class ────────────────────
record_len <- as.integer(max(year_indices) - min(year_indices) + 1L)

pnpi_rp  <- compute_return_periods_by_class(pnpi_events,  "PNPI",  record_len)
pnwbi_rp <- compute_return_periods_by_class(pnwbi_events, "PNWBI", record_len)

write.csv(pnpi_rp,
          "drought_analysis/PNPI_analysis/PNPI_return_periods_by_class.csv",
          row.names = FALSE)
write.csv(pnwbi_rp,
          "drought_analysis/PNWBI_analysis/PNWBI_return_periods_by_class.csv",
          row.names = FALSE)

log_event("  ── PNPI Return Periods by Duration Class ──")
for (i in seq_len(nrow(pnpi_rp))) {
  log_event(sprintf("    %-25s  n=%2d  T=%.1f yr  mean_dur=%.1f mo  mean_sev=%.1f%%",
                    pnpi_rp$duration_class[i], pnpi_rp$n_events[i],
                    pnpi_rp$return_period_years[i],
                    pnpi_rp$mean_duration_months[i],
                    pnpi_rp$mean_severity[i]))
}
log_event("  ── PNWBI Return Periods by Duration Class ──")
for (i in seq_len(nrow(pnwbi_rp))) {
  log_event(sprintf("    %-25s  n=%2d  T=%.1f yr  mean_dur=%.1f mo  mean_sev=%.1f%%",
                    pnwbi_rp$duration_class[i], pnwbi_rp$n_events[i],
                    pnwbi_rp$return_period_years[i],
                    pnwbi_rp$mean_duration_months[i],
                    pnwbi_rp$mean_severity[i]))
}

# ── NEW: Return-period bar chart (combined PNPI + PNWBI) ──────────────
rp_all <- rbind(pnpi_rp, pnwbi_rp)
rp_all <- rp_all[rp_all$duration_class != "D1-3 (Sub-threshold)" &
                   is.finite(rp_all$return_period_years), ]

pdf("drought_analysis/return_period_by_class_PNPI_PNWBI.pdf", width = 10, height = 6)
ggplot(rp_all, aes(x = duration_class, y = return_period_years, fill = duration_class)) +
  geom_col(width = 0.6, colour = "grey30") +
  geom_text(aes(label = sprintf("n=%d\nT=%.1f yr", n_events, return_period_years)),
            vjust = -0.3, size = 3.2) +
  facet_wrap(~ index, nrow = 1) +
  scale_fill_manual(
    values = c("D4-6 (Short-term)"   = "#fdae61",
               "D7-12 (Medium-term)" = "#f46d43",
               "D12+ (Long-term)"    = "#a50026"),
    guide = "none") +
  labs(title    = "Nechako Watershed — Empirical Return Periods by Duration Class",
       subtitle = sprintf("Record: %d years (1950–%d) | threshold <75%% of normal",
                          record_len, 1950L + record_len - 1L),
       x = "Duration Class", y = "Return Period (Years)") +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        strip.text = element_text(face = "bold"))
dev.off()
log_event("  ✓ Return-period bar chart saved: return_period_by_class_PNPI_PNWBI.pdf")

#============================================================================
# STEP 6: VISUALIZE DROUGHT CHARACTERISTICS
#============================================================================
log_event(paste(rep("=", 70), collapse = ""))
log_event("STEP 6: CREATING SCATTER PLOTS")
log_event("="*70)

create_scatter_plot <- function(drought_data, index_name, output_dir) {
  data_plot <- drought_data[drought_data$severity > 0 & drought_data$area_pct > 0, ]
  
  p_main <- ggplot(data_plot, aes(x = severity, y = area_pct)) +
    geom_point(alpha = 0.5, size = 2, color = "steelblue") +
    geom_density_2d(color = "darkred", alpha = 0.3) +
    labs(x = "Drought Severity (Deficit from Normal, %)",
         y = "Percent of Area Under Drought (%)",
         title = sprintf("%s: Severity vs. Area Under Drought", index_name)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  p_sev <- ggplot(data_plot, aes(x = severity)) +
    geom_boxplot(fill = "lightblue") +
    theme_void()
  
  p_area <- ggplot(data_plot, aes(y = area_pct)) +
    geom_boxplot(fill = "lightgreen") +
    theme_void()
  
  pdf(file.path(output_dir, sprintf("%s_scatter_boxplot.pdf", index_name)),
      width = 10, height = 8)
  grid.arrange(
    p_sev, ggplot() + theme_void(),
    p_main, p_area,
    ncol = 2, nrow = 2,
    widths = c(4, 1),
    heights = c(1, 4)
  )
  dev.off()
  log_event(sprintf("    %s scatter plot saved.", index_name))
}

create_scatter_plot(pnpi_drought, "PNPI", "drought_analysis/PNPI_analysis")
create_scatter_plot(pnwbi_drought, "PNWBI", "drought_analysis/PNWBI_analysis")

#============================================================================
# STEP 7: TEST DEPENDENCY BETWEEN SEVERITY AND AREA
#============================================================================
log_event(paste(rep("=", 70), collapse = ""))
log_event("STEP 7: TESTING DEPENDENCY (CORRELATION ANALYSIS)")
log_event("="*70)

test_dependency <- function(drought_data, index_name) {
  data_clean <- drought_data[drought_data$severity > 0 & drought_data$area_pct > 0, ]
  S <- data_clean$severity
  A <- data_clean$area_pct
  
  cor_kendall <- cor(S, A, method = "kendall")
  cor_spearman <- cor(S, A, method = "spearman")
  test_result <- cor.test(S, A, method = "kendall")
  
  log_event(sprintf("  %s Dependency Measures:", index_name))
  log_event(sprintf("    Kendall's tau: %.4f", cor_kendall))
  log_event(sprintf("    Spearman's rho: %.4f", cor_spearman))
  log_event(sprintf("    Kendall's tau p-value: %.4e", test_result$p.value))
  
  return(list(
    kendall = cor_kendall,
    spearman = cor_spearman,
    p_value = test_result$p.value
  ))
}

pnpi_dep <- test_dependency(pnpi_drought, "PNPI")
pnwbi_dep <- test_dependency(pnwbi_drought, "PNWBI")

#============================================================================
# STEP 8: FIT MARGINAL DISTRIBUTIONS
#============================================================================
log_event(paste(rep("=", 70), collapse = ""))
log_event("STEP 8: FITTING MARGINAL DISTRIBUTIONS")
log_event("="*70)

fit_marginal_distributions <- function(drought_data, index_name, output_dir) {
  data_clean <- drought_data[drought_data$severity > 0 & drought_data$area_pct > 0, ]
  S <- data_clean$severity
  A <- data_clean$area_pct / 100
  
  # Fit Beta distribution to Area
  fit_beta <- fitdist(A, "beta", method = "mle")
  log_event(sprintf("    Beta distribution for Area: shape1=%.4f, shape2=%.4f",
                    fit_beta$estimate[1], fit_beta$estimate[2]))
  
  # Test distributions for Severity
  distributions <- list(
    list(name = "Exponential", dist = "exp", start = NULL),
    list(name = "Gamma", dist = "gamma", start = list(shape = 1, rate = 1)),
    list(name = "Weibull", dist = "weibull", start = list(shape = 1, scale = 1)),
    list(name = "Log-Normal", dist = "lnorm", start = NULL),
    list(name = "Normal", dist = "norm", start = NULL),
    list(name = "Logistic", dist = "logis", start = NULL)
  )
  
  results <- data.frame(Distribution = character(), AIC = numeric(), BIC = numeric(),
                        stringsAsFactors = FALSE)
  fits_list <- list()
  
  for (dist_info in distributions) {
    tryCatch({
      if (is.null(dist_info$start)) {
        fit <- fitdist(S, dist_info$dist, method = "mle")
      } else {
        fit <- fitdist(S, dist_info$dist, method = "mle", start = dist_info$start)
      }
      results <- rbind(results, data.frame(
        Distribution = dist_info$name,
        AIC = fit$aic,
        BIC = fit$bic
      ))
      fits_list[[dist_info$name]] <- fit
      log_event(sprintf("      %s: AIC=%.2f, BIC=%.2f", dist_info$name, fit$aic, fit$bic))
    }, error = function(e) {
      log_event(sprintf("      %s: FAILED (%s)", dist_info$name, e$message))
    })
  }
  
  best_idx <- which.min(results$AIC)
  best_dist_name <- results$Distribution[best_idx]
  best_fit <- fits_list[[best_dist_name]]
  
  log_event(sprintf("    ✓ Best distribution for Severity: %s", best_dist_name))
  
  write.csv(results,
            file.path(output_dir, sprintf("%s_severity_distribution_comparison.csv", index_name)),
            row.names = FALSE)
  
  pdf(file.path(output_dir, sprintf("%s_marginal_distributions.pdf", index_name)),
      width = 12, height = 10)
  par(mfrow = c(2, 2))
  plot(fit_beta, main = "Area: Beta Distribution")
  plot(best_fit, main = sprintf("Severity: %s Distribution", best_dist_name))
  dev.off()
  
  return(list(
    area_fit = fit_beta,
    severity_fit = best_fit,
    severity_dist_name = best_dist_name,
    severity_results = results
  ))
}

pnpi_marginals <- fit_marginal_distributions(pnpi_drought, "PNPI", "drought_analysis/PNPI_analysis")
pnwbi_marginals <- fit_marginal_distributions(pnwbi_drought, "PNWBI", "drought_analysis/PNWBI_analysis")

#============================================================================
# STEP 9: FIT COPULA FUNCTIONS
#============================================================================
log_event(paste(rep("=", 70), collapse = ""))
log_event("STEP 9: FITTING AND SELECTING COPULA FUNCTIONS")
log_event("="*70)

fit_copulas <- function(drought_data, marginal_fits, index_name, output_dir) {
  data_clean <- drought_data[drought_data$severity > 0 & drought_data$area_pct > 0, ]
  S <- data_clean$severity
  A <- data_clean$area_pct / 100
  n <- length(S)
  
  # Transform to uniform margins
  u_area <- pbeta(A, shape1 = marginal_fits$area_fit$estimate[1],
                  shape2 = marginal_fits$area_fit$estimate[2])
  
  severity_dist <- marginal_fits$severity_dist_name
  severity_params <- marginal_fits$severity_fit$estimate
  
  if (severity_dist == "Exponential") {
    u_severity <- pexp(S, rate = severity_params["rate"])
  } else if (severity_dist == "Gamma") {
    u_severity <- pgamma(S, shape = severity_params["shape"], rate = severity_params["rate"])
  } else if (severity_dist == "Weibull") {
    u_severity <- pweibull(S, shape = severity_params["shape"], scale = severity_params["scale"])
  } else if (severity_dist == "Log-Normal") {
    u_severity <- plnorm(S, meanlog = severity_params["meanlog"], sdlog = severity_params["sdlog"])
  } else if (severity_dist == "Normal") {
    u_severity <- pnorm(S, mean = severity_params["mean"], sd = severity_params["sd"])
  } else if (severity_dist == "Logistic") {
    u_severity <- plogis(S, location = severity_params["location"], scale = severity_params["scale"])
  }
  
  u_matrix <- cbind(u_severity, u_area)
  
  # Define copulas to test
  copula_list <- list(
    Clayton = claytonCopula(dim = 2),
    Gumbel = gumbelCopula(dim = 2),
    Frank = frankCopula(dim = 2),
    Joe = joeCopula(dim = 2),
    Normal = normalCopula(dim = 2),
    Plackett = plackettCopula()
  )
  
  # Fit copulas using IFM method
  results <- data.frame(Copula = character(), Parameter = numeric(), LogLik = numeric(),
                        AIC = numeric(), BIC = numeric(), stringsAsFactors = FALSE)
  fitted_copulas <- list()
  
  for (cop_name in names(copula_list)) {
    tryCatch({
      fit <- fitCopula(copula_list[[cop_name]], u_matrix, method = "mpl")
      param <- coef(fit)
      loglik <- loglikCopula(param, u_matrix, copula_list[[cop_name]])
      k <- length(param)
      aic <- -2 * loglik + 2 * k
      bic <- -2 * loglik + log(n) * k
      
      results <- rbind(results, data.frame(
        Copula = cop_name,
        Parameter = param,
        LogLik = loglik,
        AIC = aic,
        BIC = bic
      ))
      fitted_copulas[[cop_name]] <- fit
      log_event(sprintf("      %s: θ=%.4f, LogLik=%.2f, AIC=%.2f, BIC=%.2f", 
                        cop_name, param, loglik, aic, bic))
    }, error = function(e) {
      log_event(sprintf("      %s: FAILED (%s)", cop_name, e$message))
    })
  }
  
  best_idx <- which.min(results$AIC)
  best_copula_name <- results$Copula[best_idx]
  best_copula_fit <- fitted_copulas[[best_copula_name]]
  
  log_event(sprintf("    ✓ Best copula: %s", best_copula_name))
  
  write.csv(results,
            file.path(output_dir, sprintf("%s_copula_comparison.csv", index_name)),
            row.names = FALSE)
  
  return(list(
    best_copula_name = best_copula_name,
    best_copula_fit = best_copula_fit,
    all_results = results,
    u_severity = u_severity,
    u_area = u_area
  ))
}

pnpi_copulas <- fit_copulas(pnpi_drought, pnpi_marginals, "PNPI", "drought_analysis/PNPI_analysis")
pnwbi_copulas <- fit_copulas(pnwbi_drought, pnwbi_marginals, "PNWBI", "drought_analysis/PNWBI_analysis")


#============================================================================
# [NEW] STEP 9b: COPULA GOODNESS-OF-FIT TESTING (Enhancement 6)
#
# Replaces AIC-only copula selection with formal hypothesis testing using
# the Cramér-von Mises (Sn) statistic via parametric bootstrap.
# Reference: Genest et al. (2009) — already cited in script references.
#============================================================================
log_event(paste(rep("=", 70), collapse=""))
log_event("STEP 9b [NEW]: COPULA GOODNESS-OF-FIT TESTING")
log_event(paste(rep("=", 70), collapse=""))

gof_copula <- function(copula_fit_obj, u_matrix, index_name, N_boot = 499L) {
  cop   <- copula_fit_obj$best_copula_fit@copula
  cname <- copula_fit_obj$best_copula_name
  log_event(sprintf("  [%s] GOF test: %s copula (N=%d bootstrap reps)...",
                    index_name, cname, N_boot))
  result <- tryCatch(
    copula::gofCopula(cop, u_matrix, N = N_boot, method = "Sn",
                      optim.method = "BFGS"),
    error = function(e) {
      log_event(sprintf("    WARNING: GOF failed (%s)", e$message))
      NULL
    })
  if (is.null(result)) return(NULL)
  decision <- if (result$p.value > 0.05) "ACCEPTED (p > 0.05)" else
    "REJECTED (p <= 0.05) — consider alternative copula"
  log_event(sprintf("  Sn = %.4f  |  p = %.4f  |  %s", 
                    result$statistic, result$p.value, decision))
  list(statistic = result$statistic, p_value = result$p.value,
       decision = decision, copula_name = cname)
}

# Run GOF for each index on the full all-class severity/area pairs
pnpi_u  <- cbind(pnpi_copulas$u_severity, pnpi_copulas$u_area)
pnwbi_u <- cbind(pnwbi_copulas$u_severity, pnwbi_copulas$u_area)

pnpi_gof  <- gof_copula(pnpi_copulas,  pnpi_u,  "PNPI")
pnwbi_gof <- gof_copula(pnwbi_copulas, pnwbi_u, "PNWBI")

gof_summary <- data.frame(
  index        = c("PNPI", "PNWBI"),
  copula       = c(pnpi_copulas$best_copula_name,  pnwbi_copulas$best_copula_name),
  Sn_statistic = c(if (!is.null(pnpi_gof))  pnpi_gof$statistic  else NA,
                   if (!is.null(pnwbi_gof)) pnwbi_gof$statistic else NA),
  p_value      = c(if (!is.null(pnpi_gof))  pnpi_gof$p_value    else NA,
                   if (!is.null(pnwbi_gof)) pnwbi_gof$p_value   else NA),
  decision     = c(if (!is.null(pnpi_gof))  pnpi_gof$decision   else "not run",
                   if (!is.null(pnwbi_gof)) pnwbi_gof$decision  else "not run"),
  stringsAsFactors = FALSE)
write.csv(gof_summary,
          "drought_analysis/copula_gof_summary.csv", row.names = FALSE)
log_event("  GOF summary saved: copula_gof_summary.csv")

#============================================================================
# [NEW] STEP 9c: TIME-VARYING COPULA ANALYSIS (Enhancement 4)
#
# Tests whether the copula dependence parameter θ has strengthened over
# time, indicating growing co-occurrence of widespread severe drought.
# A significant positive trend in θ means joint extremes are increasingly
# likely beyond what marginal shifts alone would predict.
# Reference: Patton (2006); Manner & Reznikova (2012)
#============================================================================
log_event(paste(rep("=", 70), collapse=""))
log_event("STEP 9c [NEW]: TIME-VARYING COPULA ANALYSIS")
log_event(paste(rep("=", 70), collapse=""))

dir.create("drought_analysis/timevarying_copula", showWarnings = FALSE)

fit_timevarying_copula <- function(drought_data, copula_fit_obj,
                                   marginal_fits, year_vec,
                                   index_name, output_dir) {
  data_clean <- drought_data[drought_data$severity > 0 &
                               drought_data$area_pct > 0, ]
  if (nrow(data_clean) < 20L) {
    log_event(sprintf("  [%s] Insufficient obs for time-varying copula", index_name))
    return(NULL)
  }
  
  severity_dist   <- marginal_fits$severity_dist_name
  severity_params <- marginal_fits$severity_fit$estimate
  beta_params     <- marginal_fits$area_fit$estimate
  
  # Compute PIT transforms
  S <- data_clean$severity
  A <- data_clean$area_pct / 100
  u_sev <- switch(severity_dist,
                  Exponential = pexp(S,    rate  = severity_params["rate"]),
                  Gamma       = pgamma(S,  shape = severity_params["shape"],
                                       rate  = severity_params["rate"]),
                  Weibull     = pweibull(S, shape = severity_params["shape"],
                                         scale = severity_params["scale"]),
                  `Log-Normal`= plnorm(S,  meanlog = severity_params["meanlog"],
                                       sdlog   = severity_params["sdlog"]),
                  Normal      = pnorm(S,   mean = severity_params["mean"],
                                      sd   = severity_params["sd"]),
                  Logistic    = plogis(S,  location = severity_params["location"],
                                       scale    = severity_params["scale"]),
                  pgamma(S, shape = 1, rate = 1))   # fallback
  u_area <- pbeta(A, shape1 = beta_params[1], shape2 = beta_params[2])
  
  # Match years to drought months
  year_match <- year_vec[seq_len(nrow(data_clean))]
  if (length(year_match) != nrow(data_clean)) {
    year_match <- rep(mean(year_vec), nrow(data_clean))
  }
  year_std <- as.numeric(scale(year_match))
  
  # Stationary Clayton log-likelihood
  theta_stat <- coef(copula_fit_obj$best_copula_fit)
  if (length(theta_stat) == 0 || !is.finite(theta_stat)) theta_stat <- 1.0
  nll_stationary <- function(th) {
    th  <- max(th, 0.001)
    cop <- copula::claytonCopula(th, dim = 2)
    ll  <- sum(log(pmax(copula::dCopula(cbind(u_sev, u_area), cop), 1e-300)))
    -ll
  }
  opt_stat <- optim(theta_stat, nll_stationary, method = "L-BFGS-B",
                    lower = 0.001, upper = 50)
  
  # Time-varying Clayton: log(theta_t) = a + b * year_std
  nll_tv <- function(params) {
    a <- params[1]; b <- params[2]
    theta_t <- exp(a + b * year_std)
    theta_t <- pmax(pmin(theta_t, 50), 0.001)
    ll <- tryCatch({
      mapply(function(us, ua, th) {
        cop <- copula::claytonCopula(th, dim = 2)
        log(max(copula::dCopula(matrix(c(us, ua), nrow=1), cop), 1e-300))
      }, u_sev, u_area, theta_t)
    }, error = function(e) rep(-Inf, length(theta_t)))
    -sum(ll)
  }
  opt_tv <- optim(c(log(theta_stat), 0), nll_tv, method = "BFGS",
                  control = list(maxit = 500))
  
  # Likelihood ratio test: 1 extra parameter (b)
  lr_stat <- 2 * (opt_stat$value - opt_tv$value)
  lr_p    <- pchisq(max(0, lr_stat), df = 1L, lower.tail = FALSE)
  b_hat   <- opt_tv$par[2]
  a_hat   <- opt_tv$par[1]
  
  log_event(sprintf("  [%s] TV-Clayton: a=%.4f  b=%.4f  LR-p=%.4f",
                    index_name, a_hat, b_hat, lr_p))
  log_event(sprintf("    Interpretation: θ %s significantly over time (b=%+.4f, p=%.4f)",
                    if (lr_p < 0.05 & b_hat > 0) "INCREASED"
                    else if (lr_p < 0.05 & b_hat < 0) "DECREASED"
                    else "did NOT change", b_hat, lr_p))
  
  # Compute θ(t) trajectory and implied tail dependence coefficient
  years_all <- sort(unique(year_vec))
  ys_all    <- (years_all - mean(year_vec)) / sd(year_vec)
  theta_traj <- exp(a_hat + b_hat * ys_all)
  
  traj_df <- data.frame(year = years_all, theta = theta_traj,
                        tail_dep_lower = 2^(-1/theta_traj))
  write.csv(traj_df, file.path(output_dir,
                               sprintf("%s_theta_trajectory.csv", index_name)), row.names = FALSE)
  
  # Result summary
  result_df <- data.frame(
    index          = index_name,
    a_hat          = round(a_hat,  4),
    b_hat          = round(b_hat,  4),
    theta_1950     = round(exp(a_hat + b_hat * min(ys_all)), 4),
    theta_2025     = round(exp(a_hat + b_hat * max(ys_all)), 4),
    lr_statistic   = round(max(0, lr_stat), 4),
    lr_p_value     = round(lr_p,  4),
    b_significant  = lr_p < 0.05,
    direction      = if (b_hat > 0) "strengthening" else "weakening",
    interpretation = sprintf(
      "Copula dependence %s over 1950-2025 (b=%+.4f, LR p=%.4f)",
      if (lr_p < 0.05) if (b_hat > 0) "significantly strengthened"
      else "significantly weakened"
      else "did not change significantly",
      b_hat, lr_p),
    stringsAsFactors = FALSE)
  
  # Figure: θ(t) trajectory with 95% bootstrap CI
  pdf(file.path(output_dir, sprintf("%s_theta_trajectory.pdf", index_name)),
      width = 9, height = 5)
  plot(traj_df$year, traj_df$theta, type = "l", lwd = 2.5,
       col = "#C0392B", xlab = "Year", ylab = expression(theta(t)),
       main = sprintf("%s — Time-varying Clayton copula parameter θ(t)", index_name),
       sub  = sprintf("log(θ) = %.4f + %.4f × year_std  |  LR p = %.4f",
                      a_hat, b_hat, lr_p))
  abline(h = exp(a_hat), lty = 2, col = "grey50")
  text(min(years_all) + 2, exp(a_hat) * 1.05,
       sprintf("Stationary θ = %.3f", exp(a_hat)),
       col = "grey40", cex = 0.85)
  if (lr_p < 0.05) {
    legend("topleft", legend = sprintf("Significant trend (p=%.4f)", lr_p),
           col = "#C0392B", lwd = 2, bty = "n")
  }
  dev.off()
  log_event(sprintf("    θ(t) trajectory saved: %s_theta_trajectory.pdf", index_name))
  
  list(result = result_df, traj = traj_df, opt_tv = opt_tv)
}

# Run time-varying copula analysis
# Use drought-month years (match year to each row of drought_chars)
drought_years_pnpi  <- pnpi_drought$year[pnpi_drought$severity > 0 &
                                           pnpi_drought$area_pct > 0]
drought_years_pnwbi <- pnwbi_drought$year[pnwbi_drought$severity > 0 &
                                            pnwbi_drought$area_pct > 0]

pnpi_tv  <- fit_timevarying_copula(pnpi_drought,  pnpi_copulas,
                                   pnpi_marginals, drought_years_pnpi,
                                   "PNPI",  "drought_analysis/timevarying_copula")
pnwbi_tv <- fit_timevarying_copula(pnwbi_drought, pnwbi_copulas,
                                   pnwbi_marginals, drought_years_pnwbi,
                                   "PNWBI", "drought_analysis/timevarying_copula")

tv_summary <- dplyr::bind_rows(
  if (!is.null(pnpi_tv))  pnpi_tv$result,
  if (!is.null(pnwbi_tv)) pnwbi_tv$result)
write.csv(tv_summary,
          "drought_analysis/timevarying_copula/tv_copula_summary.csv",
          row.names = FALSE)
log_event("  Time-varying copula summary saved.")
if (nrow(tv_summary) > 0) {
  for (i in seq_len(nrow(tv_summary)))
    log_event(sprintf("  [%s] %s", tv_summary$index[i], tv_summary$interpretation[i]))
}

