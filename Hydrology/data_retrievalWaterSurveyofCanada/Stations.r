library(tidyhydat)
library(dplyr)
library(readr) # Recommended for cleaner CSV writing

# Download the database if you haven't already
# download_hydat()

# Get all BC stations and filter for Nechako basin codes (08J...)
nechako_gauges <- hy_stations(prov_terr_state_loc = "BC") %>%
  filter(grepl("^08J", STATION_NUMBER))

# Write the list to a CSV file called "stations.csv"
# We exclude row names to keep the file clean
write_csv(nechako_gauges, "stations.csv")

# Print a preview to the console
print(head(nechako_gauges[, c("STATION_NUMBER", "STATION_NAME")]))
