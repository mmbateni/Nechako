####################################################################################
# w5_event_ranking.R  ·  DROUGHT EVENT RANKING & 2022-2025 CONTEXTUALIZATION
####################################################################################

source("DROUGHT_ANALYSIS_utils.R")
utils_load_packages(c("terra", "ggplot2", "dplyr", "data.table", "lubridate",
                      "openxlsx", "patchwork", "scales", "sf", "viridis",
                      "gridExtra"))   # <-- added gridExtra for the summary table

# Define the %||% operator (used later for NULL/NA fallback)
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x[1])) y else x

if (!dir.exists(WD_PATH)) stop("Working directory not found: ", WD_PATH)
setwd(WD_PATH)

####################################################################################
# CONFIGURATION
####################################################################################
RANKING_DIR <- file.path(WD_PATH, "temporal_drought", "event_ranking")
dir.create(RANKING_DIR, showWarnings = FALSE, recursive = TRUE)

# Analysis period
HISTORICAL_START <- 1950
HISTORICAL_END   <- 2025
RECENT_DROUGHT_START <- as.Date("2022-01-01")
RECENT_DROUGHT_END   <- as.Date("2025-12-31")

# Duration classes (Sheffield & Wood 2008)
DURATION_CLASSES <- list(
  short   = c(4, 6, "D4-6 (Short-term)"),
  medium  = c(7, 12, "D7-12 (Medium-term)"),
  long    = c(13, 999, "D12+ (Long-term)")
)

# Indices and scales to analyze (match w1 output)
INDICES_TO_ANALYZE <- list(
  list(index = "spi", scales = c(3, 6, 12)),
  list(index = "spei", scales = c(3, 6, 12)),
  list(index = "swei", scales = c(SWEI_SCALE))
)

# Significance level for return period
ALPHA <- 0.05

####################################################################################
# HELPER FUNCTIONS
####################################################################################

#' Classify drought duration according to Sheffield & Wood (2008)
classify_duration <- function(duration_months) {
  sapply(duration_months, function(d) {
    if (d >= 4 & d <= 6) return("D4-6 (Short-term)")
    if (d >= 7 & d <= 12) return("D7-12 (Medium-term)")
    if (d >= 13) return("D12+ (Long-term)")
    return("D1-3 (Sub-threshold)")
  })
}

#' Calculate severity (intensity × duration)
calc_severity <- function(intensity, duration) {
  intensity * duration
}

#' Estimate return period from empirical frequency
calc_return_period <- function(event_rank, n_years_record) {
  # Weibull plotting position: T = (n+1)/m
  (n_years_record + 1) / event_rank
}

#' Load event catalog from w1 output
#' FIX: use same filename pattern as w1_trend_test.R (uppercase index, hyphen)
load_event_catalog <- function(index_type, scale) {
  label <- sprintf("%s-%02d", toupper(index_type), scale)
  file_path <- file.path(TREND_DIR, "basin_extent",
                         sprintf("%s_basin_event_catalog.csv", label))
  
  if (!file.exists(file_path)) {
    cat(sprintf("  ⚠ No event catalog for %s\n", label))
    return(NULL)
  }
  
  df <- tryCatch(
    read.csv(file_path, stringsAsFactors = FALSE),
    error = function(e) {
      cat(sprintf("  ❌ Error loading %s: %s\n", label, e$message))
      return(NULL)
    }
  )
  
  if (is.null(df) || nrow(df) == 0) {
    cat(sprintf("  ⚠ Empty catalog for %s\n", label))
    return(NULL)
  }
  
  # Standardize column names
  if ("start_date" %in% names(df)) df$start_date <- as.Date(df$start_date)
  if ("end_date" %in% names(df)) df$end_date <- as.Date(df$end_date)
  
  # Add metadata (use original index_type for consistency)
  df$index_type <- index_type
  df$scale <- scale
  df$index_label <- sprintf("%s_%02d", index_type, scale)   # internal label
  
  cat(sprintf("  ✓ Loaded %s: %d events\n", label, nrow(df)))
  return(df)
}

#' Identify if event overlaps with 2022-2025 period
is_recent_drought <- function(start_date, end_date) {
  start_date <= RECENT_DROUGHT_END & end_date >= RECENT_DROUGHT_START
}

####################################################################################
# PART 1: LOAD AND MERGE ALL EVENT CATALOGS
####################################################################################
cat("\n══════════════════════════════════════\n")
cat("PART 1: Loading Event Catalogs from w3\n")
cat("══════════════════════════════════════\n")

all_events <- list()

for (idx_info in INDICES_TO_ANALYZE) {
  for (sc in idx_info$scales) {
    catalog <- load_event_catalog(idx_info$index, sc)
    if (!is.null(catalog)) {
      key <- sprintf("%s_%02d", idx_info$index, sc)
      all_events[[key]] <- catalog
    }
  }
}

if (length(all_events) == 0) {
  stop("No event catalogs found. Please run w1_trend_test.R first.")
}

# Merge all catalogs
master_events <- dplyr::bind_rows(all_events)

# Filter to historical period
master_events <- master_events %>%
  filter(start_year >= HISTORICAL_START & start_year <= HISTORICAL_END)

cat(sprintf("\n✓ Total events loaded: %d (period: %d-%d)\n",
            nrow(master_events), HISTORICAL_START, HISTORICAL_END))

# --- Harmonize event columns from w1 catalogs (minimal patch) -----------------
# If the catalogs shipped a pre-computed severity called 'mean_severity',
# adopt it as 'severity'. Then back-compute mean_int if needed.

if (!"severity" %in% names(master_events) &&
    "mean_severity" %in% names(master_events)) {
  master_events$severity <- as.numeric(master_events$mean_severity)
}

# If mean intensity is missing but we have severity and duration,
# derive I = S / D (Sheffield & Wood severity definition).
if (!"mean_int" %in% names(master_events) &&
    "severity" %in% names(master_events) &&
    "duration_months" %in% names(master_events)) {
  master_events$mean_int <- as.numeric(master_events$severity) /
    pmax(as.numeric(master_events$duration_months), 1)
}

# Safety: ensure types
master_events$duration_months <- as.integer(master_events$duration_months)
master_events$start_date      <- as.Date(master_events$start_date)
master_events$end_date        <- as.Date(master_events$end_date)
################################################################################
# PART 2: CALCULATE RANKING METRICS
################################################################################
cat("\n══════════════════════════════════════\n")
cat("PART 2: Calculating Ranking Metrics\n")
cat("══════════════════════════════════════\n")

# Calculate severity
if (!"severity" %in% names(master_events)) {
  master_events$severity <- calc_severity(master_events$mean_int, master_events$duration_months)
}

# Classify duration
master_events$duration_class <- classify_duration(master_events$duration_months)

# Mark recent drought (2022-2025)
master_events$is_recent <- mapply(is_recent_drought,
                                  master_events$start_date,
                                  master_events$end_date)

# Calculate year range for return period
n_years_record <- HISTORICAL_END - HISTORICAL_START + 1

# Rank by severity (within each index-scale combination)
master_events <- master_events %>%
  group_by(index_label) %>%
  arrange(desc(severity)) %>%
  mutate(
    rank_by_severity = row_number(),
    return_period_years = calc_return_period(rank_by_severity, n_years_record),
    percentile_rank = round(100 * (1 - rank_by_severity / n()), 1)
  ) %>%
  ungroup()

# Identify the most severe recent event
recent_events <- master_events %>% filter(is_recent)
if (nrow(recent_events) > 0) {
  most_severe_recent <- recent_events %>%
    group_by(index_label) %>%
    slice_max(severity, n = 1) %>%
    ungroup()
  
  cat(sprintf("\n✓ Identified %d recent drought events (2022-2025)\n",
              nrow(recent_events)))
  cat(sprintf("  Most severe: %s (severity=%.2f, rank=%d)\n",
              most_severe_recent$index_label[1],
              most_severe_recent$severity[1],
              most_severe_recent$rank_by_severity[1]))
} else {
  cat("\n⚠ No events detected in 2022-2025 period\n")
  most_severe_recent <- NULL
}

####################################################################################
# PART 3: SAVE RANKED CATALOG
####################################################################################
cat("\n══════════════════════════════════════\n")
cat("PART 3: Saving Ranked Catalog\n")
cat("══════════════════════════════════════\n")

# Select columns for export
export_cols <- c("index_label", "start_date", "end_date", "start_year",
                 "duration_months", "duration_class", "mean_int",
                 "severity", "rank_by_severity", "return_period_years",
                 "percentile_rank", "is_recent")

ranked_catalog <- master_events[, export_cols]
write.csv(ranked_catalog,
          file.path(RANKING_DIR, "ranked_event_catalog.csv"),
          row.names = FALSE)

cat(sprintf("✓ Saved: ranked_event_catalog.csv (%d rows)\n", nrow(ranked_catalog)))

####################################################################################
# PART 3b: DURATION-CLASS RETURN PERIOD TABLE (Sheffield & Wood 2008, Table 4 style)
####################################################################################
cat("\n══════════════════════════════════════\n")
cat("PART 3b: Duration-Class Return Period Table\n")
cat("══════════════════════════════════════\n")

rp_table <- master_events %>%
  dplyr::filter(duration_class != "D1-3 (Sub-threshold)") %>%
  dplyr::group_by(index_label, duration_class) %>%
  dplyr::summarise(
    n_events         = dplyr::n(),
    mean_duration_mo = round(mean(duration_months, na.rm = TRUE), 1),
    mean_I           = round(mean(mean_int,   na.rm = TRUE), 4),  # S&W intensity
    mean_S           = round(mean(severity,        na.rm = TRUE), 3),  # S&W severity
    max_S            = round(max(severity,         na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    # Central Weibull return period
    return_period_yrs = round(n_years_record / n_events, 1),
    # Exact Poisson 95 % confidence interval (Wilson-Hilferty / chi-squared method):
    rp_lower_95 = round( (2 * n_years_record) / qchisq(0.975, 2 * (n_events + 1)), 1 ),
    rp_upper_95 = round( (2 * n_years_record) / qchisq(0.025, 2 * n_events), 1 ),
    duration_class = factor(duration_class,
                            levels = c("D4-6 (Short-term)",
                                       "D7-12 (Medium-term)",
                                       "D12+ (Long-term)"))
  ) %>%
  dplyr::arrange(index_label, duration_class) %>%
  dplyr::select(index_label, duration_class, n_events,
                return_period_yrs, rp_lower_95, rp_upper_95,
                mean_duration_mo, mean_I, mean_S, max_S)

write.csv(rp_table,
          file.path(RANKING_DIR, "return_period_by_duration_class.csv"),
          row.names = FALSE)
cat(sprintf("✓ Saved: return_period_by_duration_class.csv (%d rows)\n", nrow(rp_table)))

# --- Figure: return period dot-range chart, one panel per index ---------------
if (nrow(rp_table) > 0) {
  p_rp <- ggplot2::ggplot(
    rp_table %>% dplyr::filter(!is.na(return_period_yrs)),
    ggplot2::aes(x = duration_class, y = return_period_yrs,
                 colour = index_label, group = index_label)
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = rp_lower_95, ymax = rp_upper_95),
      width = 0.25, linewidth = 0.8, position = ggplot2::position_dodge(0.6)
    ) +
    ggplot2::geom_point(size = 3, position = ggplot2::position_dodge(0.6)) +
    ggplot2::scale_y_continuous(
      name   = "Return Period (years)",
      breaks = c(1, 2, 5, 10, 20, 50, 100),
      trans  = "log10",
      labels = scales::comma
    ) +
    ggplot2::scale_colour_brewer(palette = "Set2", name = "Index") +
    ggplot2::labs(
      title    = "Drought Return Period by Duration Class — Nechako Watershed",
      subtitle = sprintf("Weibull empirical estimate ± 95%% Poisson CI  |  Record: %d–%d",
                         HISTORICAL_START, HISTORICAL_END),
      x        = "Duration Class (Sheffield & Wood 2008)"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(face = "bold", size = 13),
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )
  
  ggplot2::ggsave(file.path(RANKING_DIR, "FigRP_Return_Period_By_Duration_Class.png"),
                  p_rp, width = 11, height = 7, dpi = 300)
  cat("✓ Saved: FigRP_Return_Period_By_Duration_Class.png\n")
  
  # ── FigRP variant: same plot but with 2022-2025 drought overlaid ──────────────
  # For each index, the most severe event that overlaps 2022-2025 is marked as a
  # diamond (shape = 18).  Its return period comes from the Weibull ranking already
  # computed in Part 2 (most_severe_recent$return_period_years).
  # This lets the reader directly compare the 2022-2025 event against the empirical
  # frequency distribution of ALL historical droughts of the same duration class.
  if (!is.null(most_severe_recent) && nrow(most_severe_recent) > 0) {
    
    # Build a data frame of the 2022-2025 event per index, keeping only the
    # duration classes that appear in the historical rp_table axes
    valid_dc <- c("D4-6 (Short-term)", "D7-12 (Medium-term)", "D12+ (Long-term)")
    
    recent_rp_overlay <- most_severe_recent %>%
      dplyr::filter(!is.na(return_period_years),
                    duration_class %in% valid_dc) %>%
      dplyr::mutate(
        duration_class = factor(duration_class, levels = valid_dc),
        # Label: show both return period and severity rank
        rp_label = sprintf("2022–2025\n~%.0f yr", return_period_years)
      )
    
    if (nrow(recent_rp_overlay) > 0) {
      
      p_rp_2022 <- p_rp +
        # ── Large diamond for each index × duration-class combination ──────────
        # shape=18 (filled diamond) distinguishes the 2022-2025 event from the
        # historical-mean dots (shape default = 19, open circle)
        ggplot2::geom_point(
          data     = recent_rp_overlay,
          ggplot2::aes(x     = duration_class,
                       y     = return_period_years,
                       colour = index_label,
                       group  = index_label),
          size     = 5,
          shape    = 18,   # filled diamond
          position = ggplot2::position_dodge(0.6),
          show.legend = FALSE
        ) +
        # ── Text label above each diamond ──────────────────────────────────────
        ggplot2::geom_label(
          data     = recent_rp_overlay,
          ggplot2::aes(x      = duration_class,
                       y      = return_period_years,
                       colour = index_label,
                       group  = index_label,
                       label  = rp_label),
          size        = 2.8,
          vjust       = -0.45,
          fontface    = "bold",
          label.size  = 0.25,
          label.padding = ggplot2::unit(0.15, "lines"),
          position    = ggplot2::position_dodge(0.6),
          show.legend = FALSE
        ) +
        ggplot2::labs(
          title    = "Drought Return Period by Duration Class — Nechako Watershed",
          subtitle = sprintf(
            paste0("Weibull empirical estimate \u00b1 95%% Poisson CI  |  Record: %d\u2013%d  ",
                   " |  \u25c6 = 2022\u20132025 event return period"),
            HISTORICAL_START, HISTORICAL_END
          )
        )
      
      ggplot2::ggsave(
        file.path(RANKING_DIR,
                  "FigRP_Return_Period_By_Duration_Class_with_2022-2025.png"),
        p_rp_2022, width = 11, height = 7, dpi = 300
      )
      cat("✓ Saved: FigRP_Return_Period_By_Duration_Class_with_2022-2025.png\n")
    } else {
      cat("  ⚠ 2022-2025 event has no valid duration class for RP overlay – skipped\n")
    }
  } else {
    cat("  ⚠ most_severe_recent is NULL – 2022-2025 RP overlay not drawn\n")
  }
}

####################################################################################
# PART 4: CREATE PRESENTATION FIGURES
####################################################################################
cat("\n══════════════════════════════════════\n")
cat("PART 4: Creating Presentation Figures\n")
cat("══════════════════════════════════════\n")

# Use SPI-12 as primary index for presentation (most common for meteorological drought)
PRIMARY_INDEX <- "spi_12"
if (!(PRIMARY_INDEX %in% unique(master_events$index_label))) {
  PRIMARY_INDEX <- unique(master_events$index_label)[1]
  cat(sprintf("⚠ %s not available, using %s instead\n", "spi_12", PRIMARY_INDEX))
}

primary_events <- master_events %>% filter(index_label == PRIMARY_INDEX)

####################################################################################
# FIGURE 1: Top 10-15 Droughts by Severity (with 2022-2025 highlighted)
####################################################################################
cat("\n  Fig 1: Top Droughts by Severity...\n")

top_n <- 15
top_events <- primary_events %>%
  arrange(desc(severity)) %>%
  head(top_n)

# Ensure recent drought is included even if not in top 15
if (any(top_events$is_recent)) {
  plot_events <- top_events
} else {
  recent_top <- primary_events %>%
    filter(is_recent) %>%
    arrange(desc(severity)) %>%
    head(1)
  if (nrow(recent_top) > 0) {
    plot_events <- dplyr::bind_rows(top_events[-top_n], recent_top) %>%
      arrange(desc(severity))
  } else {
    plot_events <- top_events
  }
}

plot_events$label <- sprintf("%s-%s\n(%d mo)",
                             format(plot_events$start_date, "%Y"),
                             format(plot_events$end_date, "%Y"),
                             plot_events$duration_months)

p1 <- ggplot2::ggplot(plot_events,
                      ggplot2::aes(x = reorder(label, severity), y = severity)) +
  ggplot2::geom_col(data = plot_events %>% filter(!is_recent),
                    fill = "#3498db", alpha = 0.7) +
  ggplot2::geom_col(data = plot_events %>% filter(is_recent),
                    fill = "#e74c3c", alpha = 0.9) +
  ggplot2::geom_text(data = plot_events %>% filter(is_recent),
                     ggplot2::aes(label = sprintf("Rank #%d\n(%.1fth percentile)",
                                                  rank_by_severity, percentile_rank)),
                     vjust = -0.5, size = 3.5, fontface = "bold", color = "#c0392b") +
  ggplot2::coord_flip() +
  ggplot2::scale_y_continuous(labels = scales::comma) +
  ggplot2::labs(
    title = sprintf("Top %d Historical Droughts in Nechako Watershed by Severity", top_n),
    subtitle = sprintf("%s | 1950-2025 | Red = 2022-2025 Drought", toupper(PRIMARY_INDEX)),
    x = "Drought Event",
    y = "Severity (Intensity × Duration)"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 14),
    plot.subtitle = ggplot2::element_text(color = "gray40", size = 11),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    panel.grid.major.y = ggplot2::element_line(color = "gray90")
  )

ggplot2::ggsave(file.path(RANKING_DIR, "Fig1_Top_Droughts_By_Severity.png"),
                p1, width = 12, height = 8, dpi = 300)
cat("    ✓ Saved: Fig1_Top_Droughts_By_Severity.png\n")

####################################################################################
# FIGURE 2: Severity Distribution with Percentile Rank (CDF)
####################################################################################
cat("\n  Fig 2: Severity Distribution (CDF)...\n")

# Calculate empirical CDF
cdf_data <- primary_events %>%
  arrange(severity) %>%
  mutate(cdf = row_number() / n())

# Get recent drought severity for annotation
if (!is.null(most_severe_recent)) {
  recent_sev <- most_severe_recent %>%
    filter(index_label == PRIMARY_INDEX) %>%
    pull(severity)
  recent_pct <- most_severe_recent %>%
    filter(index_label == PRIMARY_INDEX) %>%
    pull(percentile_rank)
} else {
  recent_sev <- NA
  recent_pct <- NA
}

p2 <- ggplot2::ggplot(cdf_data, ggplot2::aes(x = severity, y = cdf * 100)) +
  ggplot2::geom_line(color = "#2c3e50", linewidth = 1.2) +
  ggplot2::geom_point(color = "#2c3e50", size = 1.5, alpha = 0.5) +
  {if (!is.na(recent_sev))
    ggplot2::geom_vline(xintercept = recent_sev, linetype = "dashed",
                        color = "#e74c3c", linewidth = 1.5)} +
  {if (!is.na(recent_sev))
    ggplot2::annotate("segment",
                      x = recent_sev, xend = recent_sev,
                      y = 0, yend = recent_pct,
                      color = "#e74c3c", linewidth = 1, linetype = "dotted")} +
  {if (!is.na(recent_sev))
    ggplot2::annotate("point",
                      x = recent_sev, y = recent_pct,
                      color = "#e74c3c", size = 4, shape = 21, fill = "white", stroke = 1.5)} +
  ggplot2::scale_y_continuous(breaks = seq(0, 100, by = 10)) +
  ggplot2::labs(
    title = "Historical Drought Severity Distribution",
    subtitle = sprintf("%s | 1950-2025 | 2022-2025 drought at %.1fth percentile",
                       toupper(PRIMARY_INDEX), round(recent_pct, 1)),
    x = "Drought Severity",
    y = "Cumulative Frequency (%)"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 14),
    panel.grid.major = ggplot2::element_line(color = "gray90")
  )

ggplot2::ggsave(file.path(RANKING_DIR, "Fig2_Drought_Severity_Distribution.png"),
                p2, width = 10, height = 7, dpi = 300)
cat("    ✓ Saved: Fig2_Drought_Severity_Distribution.png\n")

####################################################################################
# FIGURE 3: Drought Frequency by Duration Class (Sheffield & Wood 2008 style)
####################################################################################
####################################################################################
# FIGURE 3: Drought Frequency by Duration Class (Sheffield & Wood 2008 style)
####################################################################################
cat("\n  Fig 3: Drought Frequency by Duration Class...\n")

duration_summary <- primary_events %>%
  dplyr::group_by(duration_class) %>%
  dplyr::summarise(
    n_events      = dplyr::n(),
    mean_duration = mean(duration_months, na.rm = TRUE),
    mean_severity = mean(severity,         na.rm = TRUE),  # <-- was mean(mean_S)
    max_severity  = max(severity,          na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    return_period = n_years_record / pmax(n_events, 1), # avoid divide-by-zero
    duration_class = factor(
      duration_class,
      levels = c("D1-3 (Sub-threshold)",
                 "D4-6 (Short-term)",
                 "D7-12 (Medium-term)",
                 "D12+ (Long-term)")
    )
  )

p3 <- ggplot2::ggplot(duration_summary,
                      ggplot2::aes(x = duration_class, y = n_events, fill = duration_class)) +
  ggplot2::geom_col(alpha = 0.8, color = "black", linewidth = 0.5) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("n=%d\nReturn: %.1f yrs", n_events, return_period)),
    vjust = -0.5, size = 4, fontface = "bold"
  ) +
  ggplot2::scale_fill_brewer(palette = "YlOrRd") +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.2))) +
  ggplot2::labs(
    title = "Drought Frequency by Duration Class (Sheffield & Wood 2008)",
    subtitle = sprintf("%s Basin-Averaged | %d-%d | Nechako Watershed",
                       toupper(PRIMARY_INDEX), HISTORICAL_START, HISTORICAL_END),
    x = "Duration Class",
    y = "Number of Events"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 14),
    legend.position = "none",
    panel.grid.major.y = ggplot2::element_line(color = "gray90")
  )

ggplot2::ggsave(file.path(RANKING_DIR, "Fig3_Drought_Frequency_By_Duration.png"),
                p3, width = 10, height = 7, dpi = 300)
cat("    ✓ Saved: Fig3_Drought_Frequency_By_Duration.png\n")

####################################################################################
# FIGURE 4: Combined Dashboard (all key metrics for presentation)
####################################################################################
cat("\n  Fig 4: Creating Summary Dashboard...\n")

# Summary statistics table – now built as a proper table grob
summary_stats <- data.frame(
  Metric = c("Total Drought Events (1950-2025)",
             "Short-term Droughts (4-6 mo)",
             "Medium-term Droughts (7-12 mo)",
             "Long-term Droughts (12+ mo)",
             "Most Severe Historical Drought",
             "2022-2025 Drought Rank",
             "2022-2025 Drought Percentile",
             "2022-2025 Return Period",
             "Mean Drought Duration",
             "Max Drought Duration"),
  Value = c(
    nrow(primary_events),
    sum(primary_events$duration_class == "D4-6 (Short-term)"),
    sum(primary_events$duration_class == "D7-12 (Medium-term)"),
    sum(primary_events$duration_class == "D12+ (Long-term)"),
    sprintf("%s-%s (Severity: %.1f)",
            format(max(primary_events$start_date), "%Y"),
            format(max(primary_events$end_date), "%Y"),
            max(primary_events$severity)),
    sprintf("#%d of %d",
            most_severe_recent$rank_by_severity[1] %||% NA,
            nrow(primary_events)),
    sprintf("%.1fth percentile", most_severe_recent$percentile_rank[1] %||% NA),
    sprintf("%.1f years", most_severe_recent$return_period_years[1] %||% NA),
    sprintf("%.1f months", mean(primary_events$duration_months)),
    sprintf("%d months", max(primary_events$duration_months))
  ),
  stringsAsFactors = FALSE
)

# Convert to a grid table
table_grob <- gridExtra::tableGrob(
  summary_stats,
  rows = NULL,
  theme = gridExtra::ttheme_minimal(
    base_size = 10,
    padding = unit(c(4, 2), "mm")
  )
)

# Wrap it so patchwork can handle it
p4_table <- patchwork::wrap_elements(table_grob)

# Combine figures
dashboard <- (p1 | p2) / (p3 | p4_table) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    title = "Nechako Watershed Drought Analysis for Climate Adaptation Roundtable",
    subtitle = "Prepared for April 2025 Presentation | Based on Sheffield & Wood (2008) Methodology",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 18, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(color = "gray50", size = 12, hjust = 0.5)
    )
  )

ggplot2::ggsave(file.path(RANKING_DIR, "Fig4_Presentation_Dashboard.png"),
                dashboard, width = 16, height = 12, dpi = 300)
cat("    ✓ Saved: Fig4_Presentation_Dashboard.png\n")

####################################################################################
# PART 5: SPATIAL MAP OF 2022-2025 DROUGHT (pixel‑level severity)
####################################################################################
cat("\n══════════════════════════════════════\n")
cat("PART 5: Spatial Analysis of 2022-2025 Drought\n")
cat("══════════════════════════════════════\n")

# Load basin as SpatVector for raster extraction
basin_vect <- tryCatch(
  terra::vect(BASIN_SHP),
  error = function(e) NULL
)

if (!is.null(basin_vect) && !is.null(most_severe_recent)) {
  cat("  ✓ Basin shapefile loaded\n")
  
  # Get the exact dates of the recent drought from the catalog
  recent_event <- most_severe_recent %>%
    filter(index_label == PRIMARY_INDEX)
  
  if (nrow(recent_event) == 1) {
    start_d <- recent_event$start_date[1]
    end_d   <- recent_event$end_date[1]
    cat(sprintf("  Recent drought period: %s to %s\n", start_d, end_d))
    
    # Parse primary index into type and scale
    idx_parts <- strsplit(PRIMARY_INDEX, "_")[[1]]
    idx_type <- idx_parts[1]
    idx_scale <- as.integer(idx_parts[2])
    
    # Determine which file finder to use
    if (idx_type == "swei") {
      find_fn <- function(scale) find_swei_seasonal_files(SWEI_SEAS_DIR, scale)
      data_dir <- SWEI_SEAS_DIR
    } else {
      find_fn <- function(scale) find_seasonal_nc_files(get(paste0(toupper(idx_type), "_SEAS_DIR")), idx_type, scale)
      data_dir <- get(paste0(toupper(idx_type), "_SEAS_DIR"))
    }
    
    # Locate all monthly NetCDF files for this index/scale
    all_files <- find_fn(idx_scale)
    
    if (length(all_files) == 0) {
      cat("  ⚠ No NetCDF files found for ", PRIMARY_INDEX, "\n")
    } else {
      cat(sprintf("  Found %d monthly NetCDF files\n", length(all_files)))
      
      # Determine which layers fall within the event period
      ts_list <- list()
      for (f in all_files) {
        r <- terra::rast(f)
        dates <- extract_dates_from_nc(f, terra::nlyr(r))
        # Keep only layers with date between start_d and end_d
        keep <- which(dates >= start_d & dates <= end_d)
        if (length(keep) > 0) {
          ts_list[[length(ts_list)+1]] <- r[[keep]]
        }
      }
      
      if (length(ts_list) == 0) {
        cat("  ⚠ No data found for the exact event period\n")
      } else {
        # Combine all selected layers into a single SpatRaster
        event_stack <- do.call(c, ts_list)
        cat(sprintf("  Extracted %d months of data for the event\n", terra::nlyr(event_stack)))
        
        # Mask to basin
        event_stack <- terra::mask(event_stack, basin_vect)
        
        # Compute per-pixel metrics
        r_mean <- terra::app(event_stack, fun = mean, na.rm = TRUE)
        r_min  <- terra::app(event_stack, fun = min, na.rm = TRUE)
        r_dur  <- terra::app(event_stack, fun = function(x) sum(x < DROUGHT_ONSET, na.rm = TRUE))
        
        # Convert to data frame for ggplot
        df_mean <- as.data.frame(r_mean, xy = TRUE, na.rm = FALSE)
        colnames(df_mean)[3] <- "mean_value"
        df_min  <- as.data.frame(r_min, xy = TRUE, na.rm = FALSE)
        colnames(df_min)[3] <- "min_value"
        df_dur  <- as.data.frame(r_dur, xy = TRUE, na.rm = FALSE)
        colnames(df_dur)[3] <- "months_drought"
        
        # Merge
        df_spatial <- merge(df_mean, df_min[, c("x","y","min_value")], by = c("x","y"))
        df_spatial <- merge(df_spatial, df_dur[, c("x","y","months_drought")], by = c("x","y"))
        
        # Convert basin to sf for plotting
        basin_sf <- sf::st_as_sf(basin_vect)
        
        # Create three-panel map
        # ── Fig 5 panel A: Mean index value during the event ──────────────────
        # Values are negative for drought (SPI/SPEI onset threshold ~ -0.5 to -1.0).
        # plasma with direction = -1 maps low (most-negative = most intense) to the
        # bright/yellow end and high (near-zero = mild) to dark purple. We therefore
        # use direction = 1 (default) so that MORE NEGATIVE → darker purple, which
        # is the perceptually intuitive "darker = drier" direction.
        p_mean <- ggplot2::ggplot() +
          ggplot2::geom_raster(data = df_spatial, ggplot2::aes(x = x, y = y, fill = mean_value)) +
          ggplot2::scale_fill_viridis_c(
            option    = "plasma",
            direction = 1,        # low (most negative) = darkest, 0 = yellow → darker = drier
            name      = "Mean Index"
          ) +
          ggplot2::geom_sf(data = basin_sf, fill = NA, color = "black", size = 0.6) +
          ggplot2::coord_sf() +
          ggplot2::theme_minimal() +
          ggplot2::labs(title = "Mean Index Value",
                        subtitle = sprintf("%s to %s", start_d, end_d)) +
          ggplot2::theme(legend.position = "bottom")
        
        # ── Fig 5 panel B: Minimum (peak intensity) during the event ──────────
        # Same logic as panel A: most-negative SPI/SPEI = peak drought = darkest.
        # Using inferno (direction = 1) keeps that mapping.
        p_min <- ggplot2::ggplot() +
          ggplot2::geom_raster(data = df_spatial, ggplot2::aes(x = x, y = y, fill = min_value)) +
          ggplot2::scale_fill_viridis_c(
            option    = "inferno",
            direction = 1,        # most-negative = darkest = most intense
            name      = "Minimum Index"
          ) +
          ggplot2::geom_sf(data = basin_sf, fill = NA, color = "black", size = 0.6) +
          ggplot2::coord_sf() +
          ggplot2::theme_minimal() +
          ggplot2::labs(title = "Minimum Index (Peak Intensity)") +
          ggplot2::theme(legend.position = "bottom")
        
        # ── Fig 5 panel C: Drought duration (months in drought) ───────────────
        # months_drought is a POSITIVE count (0 = never in drought, max = all months).
        # cividis default (direction = 1) goes dark-blue (low) → yellow (high), so
        # pixels with MANY drought months appear LIGHT yellow — the wrong direction.
        # Reversing with direction = -1 makes high counts dark → darker = longer drought.
        p_dur <- ggplot2::ggplot() +
          ggplot2::geom_raster(data = df_spatial, ggplot2::aes(x = x, y = y, fill = months_drought)) +
          ggplot2::scale_fill_viridis_c(
            option    = "cividis",
            direction = -1,       # flip: high month-count = dark, few months = light
            name      = "Months in Drought"
          ) +
          ggplot2::geom_sf(data = basin_sf, fill = NA, color = "black", size = 0.6) +
          ggplot2::coord_sf() +
          ggplot2::theme_minimal() +
          ggplot2::labs(title = "Drought Duration within Event") +
          ggplot2::theme(legend.position = "bottom")
        
        # Combine panels
        p_spatial <- p_mean + p_min + p_dur +
          patchwork::plot_annotation(
            title = "2022-2025 Drought: Spatial Characteristics",
            subtitle = sprintf("%s | Pixel-level metrics", toupper(PRIMARY_INDEX))
          )
        
        ggplot2::ggsave(file.path(RANKING_DIR, "Fig5_2022-2025_Spatial_Severity.png"),
                        p_spatial, width = 15, height = 6, dpi = 300)
        cat("    ✓ Saved: Fig5_2022-2025_Spatial_Severity.png\n")
      }
    }
  } else {
    cat("  ⚠ Could not isolate recent event for spatial mapping\n")
  }
} else {
  cat("  ⚠ Basin shapefile not found or recent event missing – skipping spatial map\n")
  # Fallback placeholder
  if (!is.null(basin_sf)) {
    p_spatial <- ggplot2::ggplot() +
      ggplot2::geom_sf(data = basin_sf, fill = "#ecf0f1", color = "#2c3e50", linewidth = 0.8) +
      ggplot2::annotate("text",
                        x = sf::st_coordinates(sf::st_centroid(basin_sf))[1],
                        y = sf::st_coordinates(sf::st_centroid(basin_sf))[2],
                        label = sprintf("2022-2025 Drought\nSeverity: %.1f\nRank: #%d",
                                        most_severe_recent$severity[1] %||% NA,
                                        most_severe_recent$rank_by_severity[1] %||% NA),
                        size = 5, fontface = "bold", color = "#c0392b",
                        hjust = 0.5, vjust = 0.5) +
      ggplot2::labs(
        title = "2022-2025 Drought Spatial Context (Placeholder)",
        subtitle = "Nechako Watershed | Severity ranking relative to 1950-2025"
      ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14, hjust = 0.5),
        panel.grid = ggplot2::element_blank()
      )
    
    ggplot2::ggsave(file.path(RANKING_DIR, "Fig5_2022-2025_Spatial_Context_placeholder.png"),
                    p_spatial, width = 10, height = 8, dpi = 300)
    cat("    ✓ Saved placeholder map (actual pixel data unavailable)\n")
  }
}

####################################################################################
# PART 5b: 2022-2025 PIXEL-LEVEL SEVERITY RANK MAP
####################################################################################
cat("\n══════════════════════════════════════\n")
cat("PART 5b: 2022-2025 Pixel-Level Severity Rank Map\n")
cat("══════════════════════════════════════\n")

# Helper: extract S&W events from a single pixel time-series
extract_pixel_events_sw <- function(x, onset = DROUGHT_ONSET, end_thr = DROUGHT_END) {
  x <- x[!is.na(x)]
  if (length(x) < 4L) return(data.frame(dur = integer(0), I = numeric(0), S = numeric(0)))
  in_d  <- FALSE; s_idx <- NA_integer_
  durs  <- integer(0); I_vec <- numeric(0); S_vec <- numeric(0)
  for (j in seq_along(x)) {
    v <- x[j]
    if (!in_d && !is.na(v) && v < onset) {
      in_d <- TRUE; s_idx <- j
    } else if (in_d && !is.na(v) && v >= end_thr) {
      dur   <- j - s_idx
      I_ev  <- mean(onset - x[s_idx:(j - 1L)])
      durs  <- c(durs, dur); I_vec <- c(I_vec, I_ev); S_vec <- c(S_vec, I_ev * dur)
      in_d  <- FALSE; s_idx <- NA_integer_
    }
  }
  if (in_d && !is.na(s_idx)) {
    dur  <- length(x) - s_idx + 1L
    I_ev <- mean(onset - x[s_idx:length(x)])
    durs <- c(durs, dur); I_vec <- c(I_vec, I_ev); S_vec <- c(S_vec, I_ev * dur)
  }
  data.frame(dur = durs, I = I_vec, S = S_vec)
}

# Resolve primary index type and scale
idx_parts_p <- strsplit(PRIMARY_INDEX, "_")[[1]]
idx_type_p  <- idx_parts_p[1]
idx_scale_p <- as.integer(idx_parts_p[2])

# Part 5b reads the time-series RDS written by w1_trend_test.R (not the results
# CSV, which deliberately omits raw time-series to stay lightweight).
ts_rds_file <- file.path(TREND_DIR,
                         sprintf("%s_%02d_timeseries.rds", idx_type_p, idx_scale_p))

rank_map_done <- FALSE

if (file.exists(ts_rds_file)) {
  cat(sprintf("  Loading pixel time-series matrix: %s\n", basename(ts_rds_file)))
  
  ts_obj <- tryCatch(
    readRDS(ts_rds_file),
    error = function(e) {
      cat(sprintf("  ⚠ Could not load time-series RDS: %s\n", e$message))
      NULL
    }
  )
  
  if (!is.null(ts_obj)) {
    ts_mat_p        <- ts_obj$ts_mat
    lon_p           <- ts_obj$lon
    lat_p           <- ts_obj$lat
    record_start_yr <- ts_obj$record_start_year %||% 1950L
    n_tcols         <- ncol(ts_mat_p)
    n_pix_p         <- nrow(ts_mat_p)
    
    cat(sprintf("  Pixels: %d  |  Time columns: %d  |  Record start: %d\n",
                n_pix_p, n_tcols, record_start_yr))
    
    col_2022_start  <- (2022L - record_start_yr) * 12L + 1L
    col_2025_end    <- min((2025L - record_start_yr + 1L) * 12L, n_tcols)
    
    if (col_2022_start < 1L || col_2022_start > n_tcols) {
      cat("  ⚠ 2022-2025 columns out of range for this record – skipping rank map\n")
    } else {
      cat(sprintf("  Computing per-pixel event catalogs (onset < %.1f)...\n", DROUGHT_ONSET))
      
      rank_vec  <- rep(NA_real_, n_pix_p)
      sev_2022  <- rep(NA_real_, n_pix_p)
      n_hist    <- rep(NA_integer_, n_pix_p)
      
      for (i in seq_len(n_pix_p)) {
        full_ts   <- ts_mat_p[i, ]
        hist_ts   <- full_ts
        
        idx_22_25 <- seq(col_2022_start, col_2025_end)
        ts_22_25  <- full_ts[idx_22_25]
        
        ev_all    <- extract_pixel_events_sw(hist_ts)
        ev_recent <- extract_pixel_events_sw(ts_22_25)
        
        if (nrow(ev_recent) == 0L || nrow(ev_all) == 0L) next
        
        s_recent <- max(ev_recent$S, na.rm = TRUE)
        sev_2022[i] <- s_recent
        
        n_less      <- sum(ev_all$S < s_recent, na.rm = TRUE)
        rank_vec[i] <- 100 * n_less / nrow(ev_all)
        n_hist[i]   <- nrow(ev_all)
      }
      
      # Assemble output data frame
      rank_df <- data.frame(
        lon           = lon_p,
        lat           = lat_p,
        severity_2022 = sev_2022,
        pct_rank_2022 = rank_vec,
        n_hist_events = n_hist
      )
      write.csv(rank_df,
                file.path(RANKING_DIR, "pixel_rank_2022-2025.csv"),
                row.names = FALSE)
      cat("  ✓ Saved: pixel_rank_2022-2025.csv\n")
      
      # Build rasters
      r_base <- tryCatch(terra::rast(r_tmpl_path <- file.path(
        TREND_DIR, sprintf("%s_%02d_D12p_pixel_frequency.nc", idx_type_p, idx_scale_p)
      )), error = function(e) NULL)
      
      if (is.null(r_base)) {
        r_base <- terra::rast(
          xmin = min(rank_df$lon) - 0.05, xmax = max(rank_df$lon) + 0.05,
          ymin = min(rank_df$lat) - 0.05, ymax = max(rank_df$lat) + 0.05,
          resolution = 0.1, crs = "EPSG:4326"
        )
        cat("  ℹ Using synthetic raster extent (D12p NC not found)\n")
      }
      
      fill_raster <- function(base_r, df, col_name) {
        r_out <- base_r[[1]]
        terra::values(r_out) <- NA_real_
        cell_ids <- terra::cellFromXY(r_out, cbind(df$lon, df$lat))
        valid    <- !is.na(cell_ids)
        terra::values(r_out)[cell_ids[valid]] <- df[[col_name]][valid]
        r_out
      }
      
      r_rank <- fill_raster(r_base, rank_df, "pct_rank_2022")
      r_sev  <- fill_raster(r_base, rank_df, "severity_2022")
      
      names(r_rank) <- "pct_rank_2022_2025"
      names(r_sev)  <- "severity_2022_2025"
      
      terra::writeCDF(
        c(r_rank, r_sev),
        file.path(RANKING_DIR, "pixel_rank_2022-2025.nc"),
        overwrite = TRUE
      )
      cat("  ✓ Saved: pixel_rank_2022-2025.nc\n")
      
      # Two-panel PNG map
      df_rank <- as.data.frame(r_rank, xy = TRUE, na.rm = TRUE)
      df_sev  <- as.data.frame(r_sev,  xy = TRUE, na.rm = TRUE)
      colnames(df_rank)[3] <- "pct_rank"
      colnames(df_sev)[3]  <- "severity"
      
      basin_sf_plot <- tryCatch(
        sf::st_as_sf(terra::vect(BASIN_SHP)),
        error = function(e) NULL
      )
      
      # ──────────────────────────────────────────────────────────────────────
      # WHAT DOES "PERCENTILE RANK" MEAN IN THIS MAP?
      # ──────────────────────────────────────────────────────────────────────
      # For every grid pixel, we extracted the full time-series of the primary
      # drought index (e.g. SPI-12) from HISTORICAL_START to HISTORICAL_END.
      # Using the Sheffield & Wood (2008) event-identification algorithm
      # (drought onset < DROUGHT_ONSET, termination ≥ DROUGHT_END) we built a
      # catalogue of ALL discrete drought events at that pixel over the entire
      # record.  Each event has a severity score S = mean_intensity × duration
      # (where intensity = onset_threshold − index_value, always positive).
      #
      # For the 2022-2025 period we then isolated the sub-series (time columns
      # col_2022_start : col_2025_end) and computed the maximum severity S_recent
      # of any event detected within that window.
      #
      # The PERCENTILE RANK at each pixel is defined as:
      #
      #   pct_rank = 100 × (#historical events with S < S_recent) / #total events
      #
      # Interpretation:
      #   • pct_rank = 90  → the 2022-2025 drought was more severe than 90 % of
      #                       all historical droughts at that pixel.  In other
      #                       words it is approximately a 1-in-10 event locally.
      #   • pct_rank = 99  → exceeds 99 % of the record; roughly a 1-in-75-year
      #                       event (given a ~75-year record).
      #   • pct_rank = 0   → the 2022-2025 event is the mildest on record there.
      #   • NA (grey)      → no drought event was detected at that pixel for
      #                       either the historical period or 2022-2025 (e.g.
      #                       perennially wet pixels at high elevation).
      #
      # The colour scale uses plasma with direction = -1 so that high percentile
      # ranks appear DARK (visually prominent) and near-zero ranks appear LIGHT.
      # ──────────────────────────────────────────────────────────────────────
      # WHAT DOES "SEVERITY (I × D)" MEAN IN THIS MAP?
      # ──────────────────────────────────────────────────────────────────────
      # Severity follows Sheffield & Wood (2008, J. Climate):
      #   S = Ī × D
      # where Ī = mean intensity over the event duration and D = duration in months.
      # Intensity for a single month is: I_t = onset_threshold − index_value
      # (positive when the index is below the drought-onset threshold).
      #
      # This means a pixel can have high severity because it was:
      #   (a) extremely dry (high I) even for a short time, OR
      #   (b) persistently in drought (high D) even at moderate intensity.
      # The S map therefore captures the COMBINED space-time footprint of the
      # 2022-2025 event and complements the percentile-rank map.
      #
      # The colour scale uses inferno with direction = -1 so that high severity
      # values appear DARK — matching the convention of all other maps in this
      # analysis where "darker = more extreme drought".
      # ──────────────────────────────────────────────────────────────────────
      
      make_base_map <- function(df, aes_col, fill_label, palette, title_str) {
        p <- ggplot2::ggplot() +
          ggplot2::geom_raster(data = df,
                               ggplot2::aes(x = x, y = y, fill = .data[[aes_col]])) +
          ggplot2::scale_fill_viridis_c(
            option    = palette,
            direction = -1,     # flip palette: HIGH values = DARK, low values = light.
            # For pct_rank (0–100 %): high percentile = dark = most extreme.
            # For severity  (I × D) : high severity  = dark = worst drought.
            name      = fill_label,
            na.value  = "grey90"
          ) +
          ggplot2::coord_sf() +
          ggplot2::labs(title = title_str) +
          ggplot2::theme_bw(base_size = 11) +
          ggplot2::theme(
            plot.title      = ggplot2::element_text(face = "bold"),
            legend.position = "bottom",
            axis.title      = ggplot2::element_blank()
          )
        if (!is.null(basin_sf_plot))
          p <- p + ggplot2::geom_sf(data = basin_sf_plot, fill = NA,
                                    colour = "black", linewidth = 0.6)
        p
      }
      
      p_rank_map <- make_base_map(
        df_rank, "pct_rank",
        "Percentile rank (%)",
        "plasma",
        "2022-2025 Drought: Percentile rank vs all historical events"
      )
      p_sev_map <- make_base_map(
        df_sev, "severity",
        "Severity (I × D)",
        "inferno",
        "2022-2025 Drought: Pixel-level S&W severity"
      )
      
      p_rank_combo <- p_rank_map + p_sev_map +
        patchwork::plot_annotation(
          title    = "2022-2025 Drought — Pixel-Level Rank vs Historical Record",
          subtitle = sprintf(
            "%s  |  Percentile of 2022-2025 severity relative to all events at each pixel  |  %d–%d",
            toupper(PRIMARY_INDEX), HISTORICAL_START, HISTORICAL_END
          ),
          theme = ggplot2::theme(
            plot.title    = ggplot2::element_text(face = "bold", size = 14, hjust = 0.5),
            plot.subtitle = ggplot2::element_text(colour = "grey40", size = 10, hjust = 0.5)
          )
        )
      
      ggplot2::ggsave(
        file.path(RANKING_DIR, "Fig6_2022-2025_Pixel_Rank_Map.png"),
        p_rank_combo, width = 14, height = 7, dpi = 300
      )
      cat("  ✓ Saved: Fig6_2022-2025_Pixel_Rank_Map.png\n")
      rank_map_done <- TRUE
      
      cat(sprintf(
        "  Basin-wide: median percentile = %.1f%%  |  %.1f%% of pixels at ≥ 75th pct\n",
        median(rank_vec, na.rm = TRUE),
        100 * mean(rank_vec >= 75, na.rm = TRUE)
      ))
    }
  }
} else {
  cat(sprintf("  ⚠ Time-series RDS not found: %s\n  → Run w1_trend_test.R first, then re-run w5.\n",
              basename(ts_rds_file)))
}

if (!rank_map_done)
  cat("  ℹ Pixel rank map not produced (see warnings above)\n")

####################################################################################
# PART 6: EXPORT EXCEL SUMMARY FOR PRESENTATION
####################################################################################
cat("\n══════════════════════════════════════\n")
cat("PART 6: Creating Excel Summary\n")
cat("══════════════════════════════════════\n")

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Executive_Summary")
openxlsx::addWorksheet(wb, "Ranked_Events")
openxlsx::addWorksheet(wb, "Duration_Class_Stats")
openxlsx::addWorksheet(wb, "Return_Period_Table")

# Executive Summary
exec_summary <- data.frame(
  Question = c(
    "Q1: How common are droughts?",
    "Q2: Spatial characteristics?",
    "Q3: Are droughts increasing?",
    "Q4: How does 2022-2025 compare?"
  ),
  Answer = c(
    sprintf("Short-term: %d events (1 per %.1f yrs) | Long-term: %d events (1 per %.1f yrs)",
            sum(primary_events$duration_class == "D4-6 (Short-term)"),
            n_years_record / sum(primary_events$duration_class == "D4-6 (Short-term)"),
            sum(primary_events$duration_class == "D12+ (Long-term)"),
            n_years_record / sum(primary_events$duration_class == "D12+ (Long-term)")),
    "See Fig5 for pixel-level severity map (mean, minimum, duration during 2022-2025).",
    "Refer to w1_trend_test.R for Mann-Kendall trend analysis on frequency/intensity.",
    sprintf("Rank #%d of %d events (%.1fth percentile, ~1-in-%.1f year event)",
            most_severe_recent$rank_by_severity[1] %||% NA,
            nrow(primary_events),
            most_severe_recent$percentile_rank[1] %||% NA,
            most_severe_recent$return_period_years[1] %||% NA)
  ),
  stringsAsFactors = FALSE
)

openxlsx::writeData(wb, "Executive_Summary", exec_summary)
openxlsx::addStyle(wb, "Executive_Summary",
                   openxlsx::createStyle(textDecoration = "bold", fontSize = 12),
                   rows = 1, cols = 1:2)
openxlsx::setColWidths(wb, "Executive_Summary", cols = 1:2, widths = c(30, 80))

# Ranked Events (top 20)
top_20 <- ranked_catalog %>%
  filter(index_label == PRIMARY_INDEX) %>%
  arrange(rank_by_severity) %>%
  head(20) %>%
  select(-index_label)

openxlsx::writeData(wb, "Ranked_Events", top_20)
openxlsx::setColWidths(wb, "Ranked_Events", cols = 1:ncol(top_20), widths = 15)

# Duration Class Statistics
openxlsx::writeData(wb, "Duration_Class_Stats", duration_summary)
openxlsx::setColWidths(wb, "Duration_Class_Stats", cols = 1:ncol(duration_summary), widths = 20)

# Return Period Table
if (exists("rp_table") && nrow(rp_table) > 0) {
  openxlsx::writeData(wb, "Return_Period_Table", as.data.frame(rp_table))
  openxlsx::addStyle(
    wb, "Return_Period_Table",
    openxlsx::createStyle(textDecoration = "bold", fontSize = 11),
    rows = 1, cols = seq_len(ncol(rp_table))
  )
  openxlsx::setColWidths(wb, "Return_Period_Table",
                         cols = seq_len(ncol(rp_table)),
                         widths = c(18, 24, 10, 16, 12, 12, 16, 10, 10, 10))
  cat("✓ Return_Period_Table sheet written\n")
}

excel_file <- file.path(RANKING_DIR, "Drought_Ranking_Summary.xlsx")
openxlsx::saveWorkbook(wb, excel_file, overwrite = TRUE)
cat(sprintf("✓ Saved: %s\n", basename(excel_file)))

####################################################################################
# FINAL SUMMARY
####################################################################################
cat("\n╔══════════════════════════════════════╗\n")
cat("║  w5_event_ranking.R  COMPLETE        ║\n")
cat("╚══════════════════════════════════════╝\n\n")

cat("KEY FINDINGS FOR PRESENTATION:\n")
cat(sprintf("  • Total drought events (1950-2025): %d\n", nrow(primary_events)))
cat(sprintf("  • 2022-2025 drought rank: #%d of %d\n",
            most_severe_recent$rank_by_severity[1] %||% NA,
            nrow(primary_events)))
cat(sprintf("  • Percentile: %.1fth (more severe than %.1f%% of historical droughts)\n",
            most_severe_recent$percentile_rank[1] %||% NA,
            most_severe_recent$percentile_rank[1] %||% NA))
cat(sprintf("  • Return period: ~1-in-%.1f years\n",
            most_severe_recent$return_period_years[1] %||% NA))
cat(sprintf("  • Duration class: %s\n",
            most_severe_recent$duration_class[1] %||% "N/A"))

cat("\nOUTPUT FILES:\n")
cat("  Ranked catalog    : ", file.path(RANKING_DIR, "ranked_event_catalog.csv"), "\n")
cat("  Return period tbl : ", file.path(RANKING_DIR, "return_period_by_duration_class.csv"), "\n")
cat("  Excel summary     : ", file.path(RANKING_DIR, "Drought_Ranking_Summary.xlsx"), "\n")
cat("  Figure 1          : ", file.path(RANKING_DIR, "Fig1_Top_Droughts_By_Severity.png"), "\n")
cat("  Figure 2          : ", file.path(RANKING_DIR, "Fig2_Drought_Severity_Distribution.png"), "\n")
cat("  Figure 3          : ", file.path(RANKING_DIR, "Fig3_Drought_Frequency_By_Duration.png"), "\n")
cat("  Figure 4          : ", file.path(RANKING_DIR, "Fig4_Presentation_Dashboard.png"), "\n")
cat("  Figure RP         : ", file.path(RANKING_DIR, "FigRP_Return_Period_By_Duration_Class.png"), "\n")
cat("  Figure RP+2022    : ", file.path(RANKING_DIR, "FigRP_Return_Period_By_Duration_Class_with_2022-2025.png"), "\n")
if (file.exists(file.path(RANKING_DIR, "Fig5_2022-2025_Spatial_Severity.png")))
  cat("  Figure 5          : ", file.path(RANKING_DIR, "Fig5_2022-2025_Spatial_Severity.png"), "\n")
if (file.exists(file.path(RANKING_DIR, "Fig6_2022-2025_Pixel_Rank_Map.png")))
  cat("  Figure 6          : ", file.path(RANKING_DIR, "Fig6_2022-2025_Pixel_Rank_Map.png"), "\n")
if (file.exists(file.path(RANKING_DIR, "pixel_rank_2022-2025.nc")))
  cat("  Rank raster       : ", file.path(RANKING_DIR, "pixel_rank_2022-2025.nc"), "\n")