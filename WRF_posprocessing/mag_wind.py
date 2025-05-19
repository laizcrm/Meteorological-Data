"""
Script to extract and save wind speed data at 750 hPa from WRF model output files.

Description:
- Loads multiple WRF output NetCDF files from a specified directory.
- Matches file timestamps with a defined date range.
- Extracts staggered wind speed and direction variable.
- Destaggers the wind speed data to the model grid.
- Selects wind speed at a specific location (lat/lon indices) and pressure level (750 hPa).
- Stores extracted wind speed magnitudes over time in a DataFrame.
- Saves the resulting time series to a CSV file.
"""
import os
from glob import glob
from netCDF4 import Dataset
import pandas as pd
import wrf as wrf

wrf_files = sorted(glob("run2/wrfout_d01_*"))

# Define the expected timestamps for the files (hourly from 2016-07-13 00:00 to 2016-07-15 00:00)
expected_times = pd.date_range('2016-07-13T00:00:00', '2016-07-15T00:00:00', freq='1H').strftime('%Y-%m-%d_%H:%M:%S').tolist()

# Filter files matching the expected timestamps
matched_files = []
for file in wrf_files:
    file_time = os.path.basename(file)[16:]  # Extract timestamp substring from filename
    if file_time in expected_times:
        matched_files.append(file)

wind_speed_values = []

# Loop through matched files and extract wind speed at 750 hPa and specific lat/lon indices
for file in matched_files:
    with Dataset(file) as ds:
        # Extract staggered wind speed and direction at all times
        wspd_wdir_staggered = wrf.getvar(ds, "wspd_wdir", units="m/s")[0, :]  # Select first time index

        # Destagger wind speed from WRF grid to mass grid (axis=1 corresponds to west-east staggering)
        wspd = wrf.destagger(wspd_wdir_staggered, axis=1, meta=True)

        # Select wind speed at specific indices (south_north=290, west_east=205, bottom_top=28)
        # These correspond roughly to lat/lon and vertical level (750 hPa)
        wind_speed_point = wspd.sel(south_north=290, west_east=205, bottom_top=28).values.item()

        wind_speed_values.append(wind_speed_point)

# Save extracted wind speed time series to CSV
df_wind_speed = pd.DataFrame(wind_speed_values, columns=['WindSpeed_750hPa'])
df_wind_speed.to_csv('wrf_wind_speed_750hPa.csv', index=False)
