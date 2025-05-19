"""
Script to visualize Convective Inhibition Energy (CINe) from WRF model output.

This script reads WRF output data, extracts meteorological variables including
pressure, geopotential height, and CINe. It interpolates the geopotential height
at 500 hPa and plots the spatial distribution of CINe over a defined geographic
region using Cartopy for mapping. The plot includes administrative boundaries,
coastlines, and gridlines, with a colorbar to indicate CINe intensity.
"""
import datetime
from netCDF4 import num2date

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from netCDF4 import Dataset
import pandas as pd
import xarray as xr

from glob import glob
import sys
import os


import cartopy
import cartopy.crs as crs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader


import wrf as wrf
from wrf import getvar, latlon_coords, get_cartopy, to_np, smooth2d, interplevel, interpline, CoordPair, geo_bounds


ds = Dataset('run2/wrfout_d01_2016-07-14_22:00:00')

# Convert time variable to datetime object
time_var = getvar(ds, "times")
year = time_var.dt.year
day = time_var.dt.day
hour = time_var.dt.hour
month = time_var.dt.month

# Extract variables
pressure = getvar(ds, "pressure")
cine = getvar(ds, 'cape_2d')[1, :]  # Convective Inhibition Energy (CINe)
geopotential_height = getvar(ds, "z", units="dm")

# Interpolate geopotential height at 500 hPa level
height_500hpa = interplevel(geopotential_height, pressure, 500)

# Get latitude and longitude points
lats, lons = latlon_coords(height_500hpa)

# Get cartopy projection object
cart_proj = get_cartopy(height_500hpa)

# Define the map extent (west_lon, south_lat, east_lon, north_lat)
extent = [-80, -41, -38, -6]
interval = 5

min_lon, max_lon = extent[0], extent[2]
min_lat, max_lat = extent[1], extent[3]

# Create figure and axis with PlateCarree projection
fig = plt.figure(figsize=(10, 10))
ax = plt.axes(projection=crs.PlateCarree())
ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs.PlateCarree())

# Add background land feature with light gray shading
land = ax.add_feature(cfeature.LAND, facecolor='gray', alpha=0.4)

# Add shapefile (admin boundaries)
shapefile = list(shpreader.Reader('shapefile/ne_10m_admin_1_states_provinces.shp').geometries())
ax.add_geometries(shapefile, crs.PlateCarree(), edgecolor='gray', facecolor='none', linewidth=0.3)

# Add coastlines and borders
ax.coastlines(resolution='10m', color='black', linewidth=1)
ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=0.5)

# Add gridlines with labels (only left and bottom)
gridlines = ax.gridlines(crs=crs.PlateCarree(), color='gray', alpha=1.0, linestyle='--', linewidth=0.25,
                         xlocs=np.arange(min_lon, max_lon, interval),
                         ylocs=np.arange(min_lat, max_lat, interval),
                         draw_labels=True)
gridlines.top_labels = False
gridlines.right_labels = False

# Define contour levels for CINe shading
min_cine = np.min(cine)
max_cine = np.max(cine)
contour_interval = 50
levels = np.arange(0, 600, contour_interval)

# Plot filled contour for CINe
cine_contour = ax.contourf(to_np(lons), to_np(lats), to_np(cine), levels=levels, cmap='RdPu', extend='max',
                           transform=crs.PlateCarree())

# Add colorbar
cbar = plt.colorbar(cine_contour, ticks=np.arange(0, 600, contour_interval), format='%.0f', pad=0.01, fraction=0.04)
cbar.set_label('CINe (J/kg)', fontsize=12)

# Plot titles
plt.title('Maximum Convective Inhibition (CINe)', loc='center', fontsize=14, y=1.05)
plt.title('WRF Model - 9 Km Resolution', loc='left', fontsize=12)
plt.title(f'{pd.to_datetime(str(time_var.values), format="%Y-%m-%d %H")} UTC', fontsize=9, loc='right')

# Save figure
plt.savefig('cine_map.png', dpi=300)
plt.close()
