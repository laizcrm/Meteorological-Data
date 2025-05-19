#############################################
#    CLIMATOLOGY DATA DOWNLOAD SCRIPT
#    Dataset: reanalysis-era5-land
#    Source: Copernicus Climate Data Store (CDS)
#    API documentation: https://cds.climate.copernicus.eu/api-how-to
#############################################

import cdsapi
import os
import calendar

def download_climatology_data(year, month):
    """
    Downloads hourly ERA5-Land reanalysis data (temperature, wind, and precipitation)
    for a given year and month from the Copernicus Climate Data Store (CDS),
    and saves them in a structured folder format.
    
    Parameters:
    - year (int): The year of the data to be downloaded.
    - month (int): The month of the data to be downloaded.
    """
    
    # Define base directory to save the data
    base_dir = f"/home/fieldpro/Climatology/data/{year}/{month:02d}"
    os.makedirs(base_dir, exist_ok=True)

    # Get number of days in the month
    num_days = calendar.monthrange(year, month)[1]

    # Loop through all days and all 24 hours
    for day in range(1, num_days + 1):
        for time in [f"{h:02d}:00" for h in range(24)]:

            # Create directory for each day
            day_dir = os.path.join(base_dir, str(day).zfill(2))
            os.makedirs(day_dir, exist_ok=True)

            # Set request parameters for CDS API
            request = {
                "variable": [
                    "2m_temperature",
                    "10m_u_component_of_wind",
                    "10m_v_component_of_wind",
                    "total_precipitation"
                ],
                "year": str(year),
                "month": str(month).zfill(2),
                "day": [str(day).zfill(2)],
                "time": [time],
                "data_format": "grib",
                "download_format": "unarchived",
                # Define region of interest: [North, West, South, East]
                "area": [12, -80, -40, -20],  
                "format": "grib"
            }

            # Define output file name and path
            file_name = f"{year}_{month:02d}_{str(day).zfill(2)}_{time.replace(':', '')}.grib"
            file_path = os.path.join(day_dir, file_name)

            # Skip if file already exists
            if os.path.exists(file_path):
                print(f"File already exists: {file_name}, skipping download.")
            else:
                print(f"Downloading: {file_name}...")
                client = cdsapi.Client()
                client.retrieve("reanalysis-era5-land", request, file_path)
                print(f"Downloaded: {file_name}")

# Download data from 2020 to 2021
for year in range(2020, 2022):  
    for month in range(1, 13):  
        download_climatology_data(year, month)
