#################################################
# Monthly Climatology of ERA5-Land Reanalysis
# Dataset source: https://cds.climate.copernicus.eu
# Note: ERA5T data (latest 3 months) may differ slightly from ERA5.

# Note:
# - ERA5T data (latest 3 months) may differ slightly from ERA5.
# - Precipitation data from ERA5-Land is not ideal for representing 
#   rainfall extremes in Brazil, as regional models tend to 
#   underestimate heavy rainfall events. These datasets should be 
#   interpreted as initial estimates and not as precise values 
#   for extreme precipitation characterization.
##################################################

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.colors as mcolors
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from datetime import datetime, timedelta
import glob
import metpy.calc as mpcalc
import geopandas as gpd
import calendar
import os

# Open ERA5-Land reanalysis files (assumes multiple .nc files)
files = xr.open_mfdataset('data/era5/ecmwf_seas5_*.nc', concat_dim="valid_time", combine="nested")

# Preprocessing variables
files['tp'] = files.tp * 10000  # Convert to mm
files['mag'] = np.sqrt(files.u10**2 + files.v10**2)  # Wind magnitude
files['t2m'] = files.t2m - 273.15  # Convert temperature to Celsius

# Monthly aggregation
clima_tp = files.tp.groupby('valid_time.month').max(dim='valid_time')
clima_t2m_mean = files.t2m.groupby('valid_time.month').mean(dim='valid_time')
clima_t2m_max = files.t2m.groupby('valid_time.month').max(dim='valid_time')
clima_t2m_min = files.t2m.groupby('valid_time.month').min(dim='valid_time')
clima_mag = files.mag.groupby('valid_time.month').mean(dim='valid_time')

# Path to shapefile and map extent
shapefile_path = 'shapefiles/BR_UF_2019.shp'  # Update path as needed
extent = [-80.0, -40.0, -30.0, 10.0]  # [min_lon, min_lat, max_lon, max_lat]

# Custom color palette for precipitation
rain_colors = [
    "#b4f0f0", "#96d2fa", "#78b9fa", "#3c95f5", "#1e6deb", "#1463d2",
    "#0fa00f", "#28be28", "#50f050", "#72f06e", "#b3faaa", "#fff9aa",
    "#ffe978", "#ffc13c", "#ffa200", "#ff6200", "#ff3300", "#ff1500",
    "#c00100", "#a50200", "#870000", "#653b32"
]
cmap_rain = mcolors.ListedColormap(rain_colors)
cmap_rain.set_over('#000000')
cmap_rain.set_under('#aaaaaa')

# Define contour levels
data_min = 500
data_max = 2000
interval = 150
levels = np.arange(data_min, data_max + interval, interval)

# Generate meshgrid from coordinates
xx, yy = np.meshgrid(files.longitude.values, files.latitude.values)

##################################################
# Quick Plot: Single Month Test (e.g., January)
##################################################

fig = plt.figure(figsize=(10, 10))
ax = plt.axes(projection=ccrs.PlateCarree())

ax.set_extent([extent[0], extent[2], extent[1], extent[3]], crs=ccrs.PlateCarree())

shapefile = list(shpreader.Reader(shapefile_path).geometries())
ax.add_geometries(shapefile, crs=ccrs.PlateCarree(), edgecolor='gray', facecolor='none', linewidth=0.3)

ax.coastlines(resolution='10m', color='black', linewidth=0.8)
ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=0.5)

gl = ax.gridlines(draw_labels=True, linestyle='--', linewidth=0.25, color='gray', alpha=1.0,
                  xlocs=np.arange(-180, 180, 5), ylocs=np.arange(-90, 90, 5))
gl.top_labels = False
gl.right_labels = False

S1 = ax.contourf(xx, yy, clima_tp[0, :, :].values,
                 transform=ccrs.PlateCarree(),
                 cmap=cmap_rain,
                 extend='max',
                 alpha=0.6,
                 levels=levels)

cb2 = plt.colorbar(S1, format='%.0f', ticks=levels, pad=0.009, shrink=0.8)
plt.title('Monthly Maximum Precipitation - January', fontsize=14)
plt.show()

##################################################
# Full Monthly Maps
##################################################

fig, axs = plt.subplots(3, 4, figsize=(15, 12), subplot_kw={'projection': ccrs.PlateCarree()})
fig.subplots_adjust(hspace=0.2, wspace=0.2)

for i, ax in enumerate(axs.flatten()):
    month = clima_tp.sel(month=i + 1)

    ax.set_extent([extent[0], extent[2], extent[1], extent[3]], ccrs.PlateCarree())

    shapefile = list(shpreader.Reader(shapefile_path).geometries())
    ax.add_geometries(shapefile, ccrs.PlateCarree(), 
                      edgecolor='gray', facecolor='none', linewidth=0.3)

    ax.coastlines(resolution='10m', color='black', linewidth=0.8)
    ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=0.5)

    gl = ax.gridlines(draw_labels=True, linestyle='--', linewidth=0.25, color='gray', alpha=1.0,
                      xlocs=np.arange(-180, 180, 5), ylocs=np.arange(-90, 90, 5))
    gl.top_labels = False
    gl.right_labels = False

    S1 = ax.contourf(xx, yy, month.values,
                     transform=ccrs.PlateCarree(),
                     cmap=cmap_rain,
                     extend='max',
                     alpha=0.6,
                     levels=levels)

    month_name = calendar.month_abbr[i + 1]
    ax.set_title(f'{month_name}')

# Add shared colorbar
cbar_ax = fig.add_axes([0.2, 0.05, 0.6, 0.02])
cb = plt.colorbar(S1, cax=cbar_ax, orientation='horizontal',
                  label='Precipitation (mm)', format='%.0f', ticks=levels)

plt.suptitle('Monthly Maximum Precipitation - ERA5-Land Climatology', fontsize=16, y=1.02)
plt.savefig('monthly_precip_era5_land.png', bbox_inches='tight', dpi=300)
plt.show()
