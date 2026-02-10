##############################################
# BASIN-AVERAGED MSPI & MSPEI VISUALIZATION
# - MSPI/MSPEI are SINGLE multi-scalar indices (no timescales)
# - Generates SEPARATE PDFs for each index + comparison plot
# - Drought definition: onset < -1.0, termination >= 0.0
# - Shaded drought periods with severity labels
# - Titles EMBEDDED in plots (PDF-safe)
# - Excel summary with detailed drought event tracking
##############################################

library(openxlsx)
library(ggplot2)
library(gridExtra)
library(grid)
library(lubridate)

#---- Helper: Safe PDF creation with automatic device cleanup ----
safe_pdf <- function(filename, width = 12, height = 8) {
  while (dev.cur() > 1) try(dev.off(), silent = TRUE)
  
  tryCatch({
    pdf(filename, width = width, height = height, family = "Helvetica")
    TRUE
  }, error = function(e) {
    cat(sprintf("ERROR: Cannot create PDF '%s': %s\n", filename, e$message))
    FALSE
  })
}

#---- Paths ----
setwd("D:/Nechako_Drought/")
output_dir <- "mspi_mspei_plots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

mspi_file <- "mspi_results/mspi_basin_timeseries_1950-2025.xlsx"
mspei_file <- "mspei_results/mspei_basin_timeseries_1950-2025.xlsx"

#---- Configuration ----
drought_threshold_onset <- -1.0
drought_threshold_end <- 0.0
severe_threshold <- -1.3
extreme_threshold <- -1.6
exceptional_threshold <- -2.0
min_duration <- 2  # months

#---- DIAGNOSTIC: Check file availability ----
cat("\n--- Checking MSPI/MSPEI data files ---\n")
if (!file.exists(mspi_file)) {
  stop(sprintf("MSPI file not found: %s\n  Please run MSPI calculation script first.", mspi_file))
}
if (!file.exists(mspei_file)) {
  stop(sprintf("MSPEI file not found: %s\n  Please run MSPEI calculation script first.", mspei_file))
}
cat(sprintf("✓ MSPI file found: %s\n", basename(mspi_file)))
cat(sprintf("✓ MSPEI file found: %s\n", basename(mspei_file)))

#---- Load and prepare data ----
cat("\n--- Loading data ---\n")
mspi <- read.xlsx(mspi_file, sheet = "Basin_Timeseries")
mspei <- read.xlsx(mspei_file, sheet = "Basin_Timeseries")

# Convert date strings to Date objects
mspi$Date_formatted <- as.Date(paste0(mspi$Date, "-01"))
mspei$Date_formatted <- as.Date(paste0(mspei$Date, "-01"))

cat(sprintf("✓ MSPI data loaded: %d records from %s to %s\n", 
            nrow(mspi), 
            format(min(mspi$Date_formatted), "%Y-%m"), 
            format(max(mspi$Date_formatted), "%Y-%m")))
cat(sprintf("✓ MSPEI data loaded: %d records from %s to %s\n", 
            nrow(mspei), 
            format(min(mspei$Date_formatted), "%Y-%m"), 
            format(max(mspei$Date_formatted), "%Y-%m")))

#---- Detect drought events (enhanced version matching SPI/SPEI) ----
detect_drought_events <- function(df, onset_threshold = -1.0, termination_threshold = 0.0, min_duration = 2) {
  df <- df[order(df$Date_formatted), ]
  values <- df$Basin_Mean
  dates <- df$Date_formatted
  
  events <- data.frame(
    event_id = integer(),
    start_date = as.Date(character()),
    end_date = as.Date(character()),
    duration_months = integer(),
    min_value = numeric(),
    severity = character(),
    stringsAsFactors = FALSE
  )
  
  event_id <- 0
  i <- 1
  
  while (i <= length(values)) {
    # Check for drought onset
    if (!is.na(values[i]) && values[i] < onset_threshold) {
      event_id <- event_id + 1
      start_idx <- i
      min_val <- values[i]
      
      # Track drought until termination
      j <- i + 1
      while (j <= length(values) && !is.na(values[j]) && values[j] < termination_threshold) {
        min_val <- min(min_val, values[j], na.rm = TRUE)
        j <- j + 1
      }
      
      end_idx <- j - 1
      duration <- end_idx - start_idx + 1
      
      # Only record if meets minimum duration
      if (duration >= min_duration) {
        # Classify severity based on minimum value
        if (min_val < exceptional_threshold) {
          severity <- "Exceptional"
        } else if (min_val < extreme_threshold) {
          severity <- "Extreme"
        } else if (min_val < severe_threshold) {
          severity <- "Severe"
        } else {
          severity <- "Moderate"
        }
        
        events <- rbind(events, data.frame(
          event_id = event_id,
          start_date = dates[start_idx],
          end_date = dates[end_idx],
          duration_months = duration,
          min_value = round(min_val, 3),
          severity = severity,
          stringsAsFactors = FALSE
        ))
      }
      
      i <- end_idx + 1
    } else {
      i <- i + 1
    }
  }
  
  return(events)
}

##############################################
# PLOTTING FUNCTIONS (matching SPI/SPEI style)
##############################################

#---- Create full-period plot with drought shading (1950-2025) ----
create_timeseries_plot_full <- function(df, index_type) {
  events <- detect_drought_events(df, drought_threshold_onset, drought_threshold_end, min_duration)
  
  index_label <- toupper(index_type)
  title_text <- sprintf("Nechako Basin %s Basin-Wide Mean (1950-2025)", index_label)
  
  # Base plot
  p <- ggplot(df, aes(x = Date_formatted, y = Basin_Mean)) +
    geom_line(color = if (index_type == "mspi") "#1f78b4" else "#e31a1c", size = 0.9) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.7)
  
  # Add drought shading FIRST (so it appears behind the line)
  if (nrow(events) > 0) {
    for (i in seq_len(nrow(events))) {
      shade_color <- switch(events$severity[i],
                            "Exceptional" = "purple",
                            "Extreme" = "darkred",
                            "Severe" = "red",
                            "Moderate" = "orange")
      
      p <- p + annotate("rect",
                        xmin = events$start_date[i], 
                        xmax = events$end_date[i],
                        ymin = -Inf, ymax = Inf,
                        fill = shade_color, alpha = 0.12)
    }
  }
  
  # Add threshold lines
  p <- p +
    geom_hline(yintercept = drought_threshold_onset, linetype = "dashed", color = "orange", size = 0.8) +
    geom_hline(yintercept = severe_threshold, linetype = "dashed", color = "red", size = 0.8) +
    geom_hline(yintercept = extreme_threshold, linetype = "dashed", color = "darkred", size = 0.8) +
    annotate("text", x = as.Date("1955-01-01"), y = drought_threshold_onset, 
             label = "Moderate", hjust = 0, vjust = -0.5, size = 3.5, color = "orange", fontface = "bold") +
    annotate("text", x = as.Date("1955-01-01"), y = severe_threshold, 
             label = "Severe", hjust = 0, vjust = -0.5, size = 3.5, color = "red", fontface = "bold") +
    annotate("text", x = as.Date("1955-01-01"), y = extreme_threshold, 
             label = "Extreme", hjust = 0, vjust = -0.5, size = 3.5, color = "darkred", fontface = "bold") +
    labs(title = title_text,
         subtitle = sprintf("Drought thresholds: Moderate (< %.1f) | Severe (< %.1f) | Extreme (< %.1f) | %d drought events detected (≥%d months)", 
                            drought_threshold_onset, severe_threshold, extreme_threshold, 
                            nrow(events), min_duration),
         x = "Year", y = "Index Value") +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    ) +
    scale_x_date(date_breaks = "5 years", date_labels = "%Y", expand = expansion(mult = c(0.01, 0.01))) +
    coord_cartesian(ylim = c(-3.5, 3.5))
  
  return(p)
}

#---- Create recent-period plot with drought shading (2020-2025) ----
create_timeseries_plot_recent <- function(df, index_type) {
  df_recent <- df[df$Date_formatted >= as.Date("2020-01-01") & df$Date_formatted <= as.Date("2025-12-31"), ]
  events <- detect_drought_events(df_recent, drought_threshold_onset, drought_threshold_end, min_duration)
  
  index_label <- toupper(index_type)
  title_text <- sprintf("Nechako Basin %s Drought Events (2020-2025)", index_label)
  
  p <- ggplot(df_recent, aes(x = Date_formatted, y = Basin_Mean)) +
    geom_line(color = if (index_type == "mspi") "#1f78b4" else "#e31a1c", size = 1.1) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.7) +
    geom_hline(yintercept = drought_threshold_onset, linetype = "dashed", color = "orange", size = 0.9) +
    geom_hline(yintercept = severe_threshold, linetype = "dashed", color = "red", size = 0.9) +
    geom_hline(yintercept = extreme_threshold, linetype = "dashed", color = "darkred", size = 0.9) +
    labs(title = title_text,
         subtitle = sprintf("Drought definition: Onset < %.1f | Termination ≥ %.1f", 
                            drought_threshold_onset, drought_threshold_end),
         x = "Year", y = "Index Value",
         caption = "Shaded periods indicate drought events") +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30"),
      plot.caption = element_text(hjust = 0.5, size = 10, face = "italic", color = "gray40"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    ) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = expansion(mult = c(0.02, 0.02))) +
    coord_cartesian(ylim = c(-3.5, 3.5))
  
  # Add drought shading and severity labels
  if (nrow(events) > 0) {
    for (i in seq_len(nrow(events))) {
      shade_color <- switch(events$severity[i],
                            "Exceptional" = "purple",
                            "Extreme" = "darkred",
                            "Severe" = "red",
                            "Moderate" = "orange")
      
      p <- p + annotate("rect",
                        xmin = events$start_date[i], xmax = events$end_date[i],
                        ymin = -Inf, ymax = Inf,
                        fill = shade_color, alpha = 0.18)
      
      # Add severity label
      mid_date <- events$start_date[i] + (events$end_date[i] - events$start_date[i]) / 2
      label_y <- 3.2
      
      p <- p + annotate("text", 
                        x = mid_date, y = label_y,
                        label = sprintf("%s\n(%d mo)", events$severity[i], events$duration_months[i]),
                        size = 3.3,
                        color = shade_color,
                        fontface = "bold",
                        hjust = 0.5, vjust = 0.5)
    }
  }
  
  return(p)
}

#---- Create MSPI vs MSPEI comparison plot ----
create_comparison_plot <- function(mspi_df, mspei_df) {
  # Merge data
  merged <- merge(mspi_df[, c("Date_formatted", "Basin_Mean")], 
                  mspei_df[, c("Date_formatted", "Basin_Mean")],
                  by = "Date_formatted", suffixes = c("_MSPI", "_MSPEI"))
  
  # Reshape for ggplot
  merged_long <- data.frame(
    Date = rep(merged$Date_formatted, 2),
    Index = rep(c("MSPI", "MSPEI"), each = nrow(merged)),
    Value = c(merged$Basin_Mean_MSPI, merged$Basin_Mean_MSPEI)
  )
  
  # Calculate correlation
  corr_val <- cor(merged$Basin_Mean_MSPI, merged$Basin_Mean_MSPEI, use = "complete.obs")
  
  # Detect drought events for both indices
  mspi_events <- detect_drought_events(mspi_df, drought_threshold_onset, drought_threshold_end, min_duration)
  mspei_events <- detect_drought_events(mspei_df, drought_threshold_onset, drought_threshold_end, min_duration)
  
  title_text <- "Nechako Basin MSPI vs MSPEI Comparison (1950-2025)"
  
  # Create plot
  p <- ggplot(merged_long, aes(x = Date, y = Value, color = Index)) +
    geom_line(size = 0.75, alpha = 0.85) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.6) +
    geom_hline(yintercept = drought_threshold_onset, linetype = "dashed", color = "gray50", size = 0.5, alpha = 0.7) +
    geom_hline(yintercept = severe_threshold, linetype = "dashed", color = "gray50", size = 0.5, alpha = 0.7) +
    geom_hline(yintercept = extreme_threshold, linetype = "dashed", color = "gray50", size = 0.5, alpha = 0.7) +
    scale_color_manual(values = c("MSPI" = "#1f78b4", "MSPEI" = "#e31a1c")) +
    theme_minimal(base_size = 13) +
    labs(
      title = title_text,
      subtitle = sprintf("Correlation: r = %.3f | MSPI events: %d | MSPEI events: %d", 
                         corr_val, nrow(mspi_events), nrow(mspei_events)),
      x = "Year",
      y = "Index Value",
      color = "Index Type",
      caption = "MSPI = precipitation only | MSPEI = precipitation minus evapotranspiration"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30"),
      axis.title = element_text(face = "bold"),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 11),
      panel.grid.minor = element_blank(),
      plot.caption = element_text(size = 9, hjust = 0.5, face = "italic", color = "gray40")
    ) +
    scale_x_date(date_breaks = "5 years", date_labels = "%Y", expand = expansion(mult = c(0.01, 0.01))) +
    coord_cartesian(ylim = c(-3.5, 3.5))
  
  return(p)
}

##############################################
# MAIN PROCESSING
##############################################

cat("\n============================================================\n")
cat("BASIN-AVERAGED MSPI & MSPEI VISUALIZATION\n")
cat("DROUGHT DEFINITION: Onset < -1.0, Termination >= 0.0\n")
cat(sprintf("MINIMUM DURATION: %d months\n", min_duration))
cat("============================================================\n")

# Close any lingering graphics devices
while (dev.cur() > 1) try(dev.off(), silent = TRUE)

for (index_type in c("mspi", "mspei")) {
  index_label <- toupper(index_type)
  cat(sprintf("\n>>>>> Processing %s <<<<<\n", index_label))
  
  # Get the appropriate dataset
  df <- if (index_type == "mspi") mspi else mspei
  
  cat(sprintf("  - Data range: %s to %s (n=%d)\n", 
              format(min(df$Date_formatted), "%Y-%m"), 
              format(max(df$Date_formatted), "%Y-%m"), 
              nrow(df)))
  
  # Save basin average as CSV (matching SPI/SPEI format)
  write.csv(df, file.path(output_dir, sprintf("%s_basin_average_1950_2025.csv", index_type)), 
            row.names = FALSE)
  cat(sprintf("  + Saved: %s_basin_average_1950_2025.csv\n", index_type))
  
  #---- FULL PERIOD (1950-2025) ----
  p_full <- create_timeseries_plot_full(df, index_type)
  pdf_file <- file.path(output_dir, sprintf("%s_full_period_1950_2025.pdf", toupper(index_type)))
  
  if (safe_pdf(pdf_file, 12, 8)) {
    print(p_full)
    dev.off()
    cat(sprintf("  + Saved: %s\n", basename(pdf_file)))
  }
  
  #---- RECENT PERIOD (2020-2025) with drought events ----
  df_recent <- df[df$Date_formatted >= as.Date("2020-01-01") & df$Date_formatted <= as.Date("2025-12-31"), ]
  events <- detect_drought_events(df_recent, drought_threshold_onset, drought_threshold_end, min_duration)
  
  if (nrow(events) > 0) {
    write.csv(events, file.path(output_dir, sprintf("%s_drought_events_2020_2025.csv", index_type)), 
              row.names = FALSE)
    cat(sprintf("  + Detected %d drought events (min duration %d months)\n", nrow(events), min_duration))
  } else {
    cat(sprintf("  + No drought events detected in 2020-2025 period (min duration %d months)\n", min_duration))
  }
  
  p_recent <- create_timeseries_plot_recent(df, index_type)
  pdf_file <- file.path(output_dir, sprintf("%s_recent_period_2020_2025.pdf", toupper(index_type)))
  
  if (safe_pdf(pdf_file, 12, 8)) {
    print(p_recent)
    dev.off()
    cat(sprintf("  + Saved: %s\n", basename(pdf_file)))
  }
  
  # Clean up graphics device
  while (dev.cur() > 1) try(dev.off(), silent = TRUE)
}

#---- COMPARISON PLOT ----
cat("\n>>>>> Creating MSPI vs MSPEI Comparison Plot <<<<<\n")
p_comparison <- create_comparison_plot(mspi, mspei)
pdf_file <- file.path(output_dir, "MSPI_MSPEI_comparison_1950_2025.pdf")

if (safe_pdf(pdf_file, 14, 9)) {
  print(p_comparison)
  dev.off()
  cat(sprintf("  + Saved: %s\n", basename(pdf_file)))
}

##############################################
# GENERATE SUMMARY STATISTICS EXCEL FILE
##############################################

cat("\n--- Generating Summary Statistics ---\n")

# Detect drought events for full period
mspi_droughts_full <- detect_drought_events(mspi, drought_threshold_onset, drought_threshold_end, min_duration)
mspei_droughts_full <- detect_drought_events(mspei, drought_threshold_onset, drought_threshold_end, min_duration)

# Recent period
mspi_recent <- mspi[mspi$Date_formatted >= as.Date("2020-01-01"), ]
mspei_recent <- mspei[mspei$Date_formatted >= as.Date("2020-01-01"), ]
mspi_droughts_recent <- detect_drought_events(mspi_recent, drought_threshold_onset, drought_threshold_end, min_duration)
mspei_droughts_recent <- detect_drought_events(mspei_recent, drought_threshold_onset, drought_threshold_end, min_duration)

summary_wb <- createWorkbook()

#---- Summary Statistics Sheet ----
stats_data <- data.frame(
  Index = c("MSPI", "MSPEI"),
  Mean = c(mean(mspi$Basin_Mean, na.rm = TRUE), mean(mspei$Basin_Mean, na.rm = TRUE)),
  Median = c(median(mspi$Basin_Mean, na.rm = TRUE), median(mspei$Basin_Mean, na.rm = TRUE)),
  StdDev = c(sd(mspi$Basin_Mean, na.rm = TRUE), sd(mspei$Basin_Mean, na.rm = TRUE)),
  Min = c(min(mspi$Basin_Mean, na.rm = TRUE), min(mspei$Basin_Mean, na.rm = TRUE)),
  Max = c(max(mspi$Basin_Mean, na.rm = TRUE), max(mspei$Basin_Mean, na.rm = TRUE)),
  Drought_Events_Full_Period = c(nrow(mspi_droughts_full), nrow(mspei_droughts_full)),
  Drought_Events_2020_2025 = c(nrow(mspi_droughts_recent), nrow(mspei_droughts_recent))
)

addWorksheet(summary_wb, "Summary_Statistics")
writeData(summary_wb, "Summary_Statistics", stats_data, rowNames = FALSE)
addStyle(summary_wb, "Summary_Statistics", 
         style = createStyle(textDecoration = "bold", fgFill = "#D3D3D3"), 
         rows = 1, cols = 1:ncol(stats_data))

#---- MSPI Drought Events (Full Period) ----
if (nrow(mspi_droughts_full) > 0) {
  addWorksheet(summary_wb, "MSPI_Droughts_Full")
  writeData(summary_wb, "MSPI_Droughts_Full", mspi_droughts_full)
  addStyle(summary_wb, "MSPI_Droughts_Full", 
           style = createStyle(textDecoration = "bold", fgFill = "#D3D3D3"), 
           rows = 1, cols = 1:ncol(mspi_droughts_full))
}

#---- MSPEI Drought Events (Full Period) ----
if (nrow(mspei_droughts_full) > 0) {
  addWorksheet(summary_wb, "MSPEI_Droughts_Full")
  writeData(summary_wb, "MSPEI_Droughts_Full", mspei_droughts_full)
  addStyle(summary_wb, "MSPEI_Droughts_Full", 
           style = createStyle(textDecoration = "bold", fgFill = "#D3D3D3"), 
           rows = 1, cols = 1:ncol(mspei_droughts_full))
}

#---- MSPI Drought Events (2020-2025) ----
if (nrow(mspi_droughts_recent) > 0) {
  addWorksheet(summary_wb, "MSPI_Droughts_2020_2025")
  writeData(summary_wb, "MSPI_Droughts_2020_2025", mspi_droughts_recent)
  addStyle(summary_wb, "MSPI_Droughts_2020_2025", 
           style = createStyle(textDecoration = "bold", fgFill = "#D3D3D3"), 
           rows = 1, cols = 1:ncol(mspi_droughts_recent))
}

#---- MSPEI Drought Events (2020-2025) ----
if (nrow(mspei_droughts_recent) > 0) {
  addWorksheet(summary_wb, "MSPEI_Droughts_2020_2025")
  writeData(summary_wb, "MSPEI_Droughts_2020_2025", mspei_droughts_recent)
  addStyle(summary_wb, "MSPEI_Droughts_2020_2025", 
           style = createStyle(textDecoration = "bold", fgFill = "#D3D3D3"), 
           rows = 1, cols = 1:ncol(mspei_droughts_recent))
}

excel_file <- file.path(output_dir, "MSPI_MSPEI_Summary_Statistics.xlsx")
saveWorkbook(summary_wb, excel_file, overwrite = TRUE)
cat(sprintf("+ Saved summary statistics: %s\n", basename(excel_file)))

##############################################
# FINAL SUMMARY
##############################################

cat("\n*** ALL VISUALIZATIONS COMPLETE ***\n")
cat(sprintf("Output directory: %s\n", normalizePath(output_dir)))
cat(sprintf("Generated:\n"))
cat(sprintf("  - 5 PDF files (2 indices × 2 periods + 1 comparison)\n"))
cat(sprintf("  - 2 CSV basin average files\n"))
cat(sprintf("  - %d CSV drought event files\n", 
            sum(nrow(mspi_droughts_recent) > 0, nrow(mspei_droughts_recent) > 0)))
cat(sprintf("  - 1 Excel summary statistics file\n"))
cat(sprintf("\nDrought Event Summary:\n"))
cat(sprintf("  MSPI  - Full period: %d events | Recent (2020-2025): %d events\n", 
            nrow(mspi_droughts_full), nrow(mspi_droughts_recent)))
cat(sprintf("  MSPEI - Full period: %d events | Recent (2020-2025): %d events\n", 
            nrow(mspei_droughts_full), nrow(mspei_droughts_recent)))