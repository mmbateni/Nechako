# 1. Load Required Libraries
library(tidyverse)
library(lubridate)
library(scales)

# 2. Configuration & Paths
setwd("D:/Nechako_Drought/")

met_data_folder   <- "./Meteorology/data_retrievalECCCanada/climate_daily_downloads/" 
output_folder     <- "./Meteorology/primary_plots/" 
met_metadata_file <- "./Meteorology/data_retrievalECCCanada/Stations.csv" 

if(!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# 3. Load Station Metadata
met_metadata <- read_csv(met_metadata_file, show_col_types = FALSE) %>%
  select(STATION_ID, STATION_NAME) %>%
  mutate(STATION_ID = as.character(STATION_ID)) %>%
  distinct()

# 4. Generate File List
met_files <- list.files(path = met_data_folder, pattern = "_CLIMATE_DAILY\\.csv$", full.names = TRUE)

# 5. Helper Function: Detect Data Gaps
get_met_gaps <- function(df, check_col) {
  if(nrow(df) < 2) return(list(gaps = NULL, full_data = df))
  full_dates <- seq(min(df$Date), max(df$Date), by = "day")
  df_full <- tibble(Date = full_dates) %>% left_join(df, by = "Date")
  is_missing <- is.na(df_full[[check_col]])
  if(!any(is_missing)) return(list(gaps = NULL, full_data = df_full))
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

# 6. Main Processing Function
process_met_station <- function(file_path, metadata_df) {
  
  raw_met <- read_csv(file_path, show_col_types = FALSE)
  if(nrow(raw_met) == 0) return(NULL)
  
  file_name <- basename(file_path)
  stn_id <- str_extract(file_name, "^[0-9]+") 
  stn_info <- metadata_df %>% filter(STATION_ID == stn_id)
  stn_name <- if(nrow(stn_info) > 0) stn_info$STATION_NAME[1] else stn_id
  clean_stn_name <- str_replace_all(stn_name, "[^[:alnum:]]", "_")
  
  df <- raw_met %>% mutate(Date = as.Date(Date)) %>% arrange(Date)
  
  # Determine axis breaks based on record length
  # If the record is very long (>30 years), show every 2 years. Otherwise, show every 1 year.
  year_range <- as.numeric(difftime(max(df$Date), min(df$Date), units = "days")) / 365
  break_interval <- if(year_range > 30) "2 years" else "1 year"
  
  # --- PLOT 1: TEMPERATURES ---
  df_temp <- df %>% select(Date, Tmax_C, Tmin_C, Tmean_C)
  res_temp <- get_met_gaps(df_temp, "Tmean_C")
  df_temp_long <- res_temp$full_data %>%
    pivot_longer(cols = starts_with("T"), names_to = "Variable", values_to = "Temp")
  
  p_temp <- ggplot() +
    {if(!is.null(res_temp$gaps)) geom_rect(data = res_temp$gaps, 
                                           aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                                           fill = "red", alpha = 0.15)} +
    geom_line(data = df_temp_long, aes(x = Date, y = Temp, color = Variable), alpha = 0.8) +
    scale_color_manual(values = c("Tmax_C" = "firebrick", "Tmin_C" = "royalblue", "Tmean_C" = "black")) +
    scale_x_date(date_breaks = break_interval, date_labels = "%Y", expand = c(0.01, 0.01)) + 
    theme_minimal() +
    labs(title = paste(stn_id, "-", stn_name), subtitle = "Temperature Records (°C) | Red shading = Gaps", 
         y = "Temperature (°C)", x = "Year") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
          panel.grid.major.x = element_line(color = "gray90"),
          legend.position = "bottom")
  
  ggsave(paste0(output_folder, stn_id, "_", clean_stn_name, "_TEMPERATURE.png"), p_temp, width = 16, height = 8)
  
  # --- PLOT 2: PRECIPITATION BREAKDOWN ---
  df_precip_long <- df %>%
    select(Date, Rain_mm, Snow_cm) %>%
    pivot_longer(cols = c(Rain_mm, Snow_cm), names_to = "Type", values_to = "Amount") %>%
    mutate(Type = factor(Type, levels = c("Snow_cm", "Rain_mm"), labels = c("Snow (cm)", "Rain (mm)")))
  
  res_precip <- get_met_gaps(df %>% select(Date, Total_Precip_mm), "Total_Precip_mm")
  
  p_precip <- ggplot() +
    {if(!is.null(res_precip$gaps)) geom_rect(data = res_precip$gaps, 
                                             aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf),
                                             fill = "red", alpha = 0.1)} +
    geom_col(data = df_precip_long, aes(x = Date, y = Amount, fill = Type), width = 1, alpha = 0.9) +
    scale_fill_manual(values = c("Rain (mm)" = "dodgerblue4", "Snow (cm)" = "skyblue")) +
    scale_y_continuous(trans = pseudo_log_trans(base = 10),
                       breaks = c(0, 1, 5, 10, 25, 50, 100, 200),
                       labels = c("0", "1", "5", "10", "25", "50", "100", "200")) +
    scale_x_date(date_breaks = break_interval, date_labels = "%Y", expand = c(0.01, 0.01)) + 
    theme_minimal() +
    labs(title = paste(stn_id, "-", stn_name), 
         subtitle = "Daily Precipitation Breakdown | Dark Blue = Rain (mm), Light Blue = Snow (cm)",
         y = "Amount (Log Scale)", x = "Year", fill = "Type") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
          panel.grid.major.x = element_line(color = "gray90"),
          legend.position = "bottom")
  
  ggsave(paste0(output_folder, stn_id, "_", clean_stn_name, "_PRECIPITATION.png"), p_precip, width = 16, height = 8)
  
  message(paste("Saved plots for:", stn_name))
}

# 7. Execution Logic
if (exists("met_files") && length(met_files) > 0) {
  walk(met_files, function(f) {
    tryCatch({
      process_met_station(f, met_metadata)
    }, error = function(e) {
      message(paste("Error in file:", f, "| Details:", e$message))
    })
  })
} else {
  # If the variable was lost, try to recreate it once more before failing
  met_files <- list.files(path = met_data_folder, pattern = "_CLIMATE_DAILY\\.csv$", full.names = TRUE)
  if(length(met_files) > 0) {
    walk(met_files, ~try(process_met_station(.x, met_metadata)))
  } else {
    message("No files found. Please verify your folder path.")
  }
}