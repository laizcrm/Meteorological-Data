"""
This script performs numerical analysis and visualization of atmospheric data based on weather model outputs. 
It loads meteorological variables such as temperature, humidity, and wind speed from NetCDF files, processes the data 
to compute derived parameters, and generates visual plots for exploratory analysis and reporting. The script is designed 
to be modular and extensible for various meteorological datasets and includes functionality to handle time series 
and spatial data slices efficiently.
"""

import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Point
import cartopy.crs as ccrs

# ===============================================================
#                      Local Weather Stations Data
# ===============================================================

local_stations = pd.read_csv('/home/data_provider/monthly_accumulated/2024-10-01_2024-10-31_precipitation_data.csv')

geometry = [Point(xy) for xy in zip(local_stations['lon'], local_stations['lat'])]
gdf_local_stations = gpd.GeoDataFrame(local_stations, geometry=geometry)
gdf_local_stations.set_crs("EPSG:4326", inplace=True)
gdf_local_stations.to_crs(shape.crs, inplace=True)

merged_local = gpd.sjoin(gdf_local_stations, shape, how="inner")

print('Number of local stations =', f'{local_stations.shape[0]} stations')


# ===============================================================
#                      CEMADEN Stations Data
# ===============================================================

cemaden_data = pd.read_csv('/home/cemaden/data/cemaden_stations.csv')

geometry_cemaden = [Point(xy) for xy in zip(cemaden_data['longitude'], cemaden_data['latitude'])]
gdf_cemaden = gpd.GeoDataFrame(cemaden_data, geometry=geometry_cemaden)
gdf_cemaden.set_crs("EPSG:4326", inplace=True)
gdf_cemaden.to_crs(shape.crs, inplace=True)

merged_cemaden = gpd.sjoin(gdf_cemaden, shape, how="inner")

print('Number of CEMADEN stations =', f'{cemaden_data.shape[0]} stations')


# ===============================================================
#                         Data Comparison Plot
# ===============================================================

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10, 10))
ax.coastlines()

# Plot local stations
ax.scatter(merged_local['lon'], merged_local['lat'], color='blue', label='Local Stations', s=10, transform=ccrs.PlateCarree())

# Plot CEMADEN stations
ax.scatter(merged_cemaden['longitude'], merged_cemaden['latitude'], color='orange', label='CEMADEN', s=10, transform=ccrs.PlateCarree())

ax.set_title("Local Stations and CEMADEN Stations", fontsize=14)
ax.legend()

plt.show()


# ===============================================================
#                  Save figures and data as needed
# ===============================================================

plt.savefig(f'./figures/precipitation_{row["station_id"]}.png')
