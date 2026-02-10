# Load libraries
library(tidyverse)
library(lubridate)
library(scales)

# 1. Configuration
setwd("D:/Nechako_Drought/")

data_folder <- "./Hydrology/SWE_data/"
output_folder <- "./Hydrology/SWE_plots/"

if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

files <- list.files(path = data_folder, pattern = "\\.csv$", full.names = TRUE)

# 2. Helper: Detect Missing Gaps
get_data_gaps <- function(df) {
  if (nrow(df) < 2) return(list(gaps = NULL, full_data = df))
  
  dt_diff <- as.numeric(difftime(df$Date[-1], df$Date[-nrow(df)], units = "days"))
  median_diff <- median(dt_diff, na.rm = TRUE)
  seq_freq <- if (median_diff > 25) "month" else "day"
  
  full_dates <- seq(min(df$Date), max(df$Date), by = seq_freq)
  
  df_full <- tibble(Date = full_dates) %>%
    left_join(df, by = "Date")
  
  is_missing <- is.na(df_full$SWE_mm)
  if (!any(is_missing)) return(list(gaps = NULL, full_data = df_full))
  
  runs <- rle(is_missing)
  end_idx <- cumsum(runs$lengths)
  start_idx <- c(1, head(end_idx, -1) + 1)
  missing_run_indices <- which(runs$values == TRUE)
  
  gaps <- tibble(
    xmin = df_full$Date[start_idx[missing_run_indices]],
    xmax = df_full$Date[end_idx[missing_run_indices]]
  )
  
  gaps <- gaps %>% mutate(xmax = if_else(xmin == xmax, xmax + hours(12), xmax))
  return(list(gaps = gaps, full_data = df_full))
}

# 3. Main Plotting Function
process_and_plot_SWE <- function(file_path) {
  # Read with manual column names and skip first row
  raw_data <- read_csv(file_path, col_names = c("Date", "SWE_mm", "Grade"), skip = 1, show_col_types = FALSE)
  if (nrow(raw_data) == 0) return(NULL)
  print(paste("Reading:", file_path))
  print(head(raw_data))
  
  file_name <- basename(file_path)
  station_name <- str_remove(file_name, "\\.csv$")
  
  df_clean <- raw_data %>%
    mutate(Date = parse_date_time(Date, orders = c("ymd", "mdy", "dmy", "Ymd HMS"))) %>%
    filter(!is.na(Date)) %>%
    arrange(Date)
  
  
  res <- get_data_gaps(df_clean)
  df_plot <- res$full_data
  gaps_df <- res$gaps
  
  # Grade colors
  grade_colors <- c("Good" = "#1b9e77", "Estimated" = "#d95f02", "Poor" = "#7570b3", "Missing" = "gray60")
  
  p <- ggplot() +
    {if (!is.null(gaps_df)) geom_rect(data = gaps_df, 
                                      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                                      fill = "red", alpha = 0.20)} +
    geom_line(data = df_plot, aes(x = Date, y = SWE_mm, color = Grade), linewidth = 0.6) +
    scale_color_manual(values = grade_colors, na.value = "gray50") +
    theme_minimal() +
    scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
    labs(
      title = paste(station_name),
      subtitle = "Snow Water Equivalent (mm) | Red shading = Missing Data | Color = Data Grade",
      y = "SWE (mm)",
      x = "Year",
      color = "Grade"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "bottom"
    )
  
  clean_name <- str_replace_all(station_name, "[^[:alnum:]]", "_")
  plot_filename <- paste0(output_folder, clean_name, "_SWE_graded.png")
  ggsave(filename = plot_filename, plot = p, width = 12, height = 6)
  
  message(paste("Processed:", station_name))
}

# 4. Run for all files
walk(files, function(f) {
  tryCatch({
    process_and_plot_SWE(f)
  }, error = function(e) {
    message(paste("Error in file:", f, "|", e$message))
  })
})
