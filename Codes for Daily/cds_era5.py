import cdsapi
import calendar
import os

# Settings
# dataset = "reanalysis-era5-single-levels"
dataset = "reanalysis-era5-land"
variable_name = "snowfall"
#"snow_depth_water_equivalent"
#
#"surface_thermal_radiation_downwards"
#"10m_v_component_of_wind"
#"surface_net_solar_radiation"
#"total_precipitation"
#"2m_temperature"
#"2m_dewpoint_temperature"
#"10m_u_component_of_wind"
#
#
#
#
#"total_cloud_cover"
 
# "total_cloud_cover"  # Set the variable to download
#base_download_directory = r"D:\ERA5 download-Blackwater"  # For blackwater ERA5
# base_download_directory = r"D:\ERA5 download"  # For Nechako ERA5
base_download_directory = r"D:\ERA5-Land download-N and B"  # For Nechako & Blackwater ERA5-Land
output_directory = os.path.join(base_download_directory, variable_name)

# Create output directory if it doesn't exist
os.makedirs(output_directory, exist_ok=True)

start_year = 2024
end_year = 2024

# Initialize the CDS API client
client = cdsapi.Client()

# Loop over years and months
for year in range(start_year, end_year + 1):
    for month_num in range(1, 11):
        num_days = calendar.monthrange(year, month_num)[1]
        days_in_month = [f"{day:02d}" for day in range(1, num_days + 1)]
        month_str = f"{month_num:02d}"

        # Define the request
        request = {
            "product_type": "reanalysis",
            "variable": variable_name,
            "year": str(year),
            "month": month_str,
            "day": days_in_month,
            "time": [
                "00:00", "01:00", "02:00", "03:00", "04:00", "05:00",
                "06:00", "07:00", "08:00", "09:00", "10:00", "11:00",
                "12:00", "13:00", "14:00", "15:00", "16:00", "17:00",
                "18:00", "19:00", "20:00", "21:00", "22:00", "23:00"
            ],
            "format": "netcdf",
            #"area": [56.5, -128.5, 52.5, -122.0] # [North, West, South, East]Nechako
            #"area": [53.9, -126.1, 51.9, -122.0] # [North, West, South, East]Blackwater
            "area": [56.5, -128.5, 51.9, -122.0] # [North, West, South, East]Nechako&Blackwater
        }

        output_filename = os.path.join(output_directory, f"era5_{variable_name}_{year}_{month_str}.nc")

        print(f"Requesting data for {year}-{month_str}")
        print(f"Saving to: {output_filename}")

        try:
            client.retrieve(dataset, request, output_filename)
            print(f"Successfully downloaded {year}-{month_str}")
        except Exception as e:
            print(f"Error downloading {year}-{month_str}: {e}")

print("Download complete.")
