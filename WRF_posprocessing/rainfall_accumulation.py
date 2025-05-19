"""
Script to plot accumulated precipitation from WRF output using Python.

- Loads WRF output NetCDF file.
- Extracts accumulated rainfall variables (RAINC + RAINNC).
- Defines geographic extent for the plot.
- Uses Cartopy to plot precipitation with state boundaries and coastlines.
- Applies a custom colormap and contour levels for precipitation shading.
- Saves the figure with timestamp in the filename.
"""

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
import pandas as pd
from wrf import getvar, ALL_TIMES, latlon_coords, to_np, extract_vars
import cartopy.io.shapereader as shpreader

# Load WRF output file
wrfin = Dataset('run2/wrfout_d01_2016-07-15_12:00:00')

# Extract precipitation variables for all times
precip_vars = extract_vars(wrfin, ALL_TIMES, ("RAINC", "RAINNC"))

# Extract time variable
times = getvar(wrfin, "times")
year = times.dt.year
day = times.dt.day
hour = times.dt.hour
month = times.dt.month

# Calculate total accumulated rainfall (convective + non-convective)
accum_rain = precip_vars["RAINC"] + precip_vars["RAINNC"]

# Get latitude and longitude arrays
lats, lons = latlon_coords(accum_rain)

# Define map extent: [min_lon, max_lon, min_lat, max_lat]
extent = [-80, -41, -38, -6]

# Grid interval for longitude and latitude gridlines
grid_interval = 5
min_lon, max_lon, min_lat, max_lat = extent[0], extent[1], extent[2], extent[3]

# Create figure and axis with PlateCarree projection
fig = plt.figure(figsize=(10, 10))
ax = plt.axes(projection=ccrs.PlateCarree())

# Set plot extent
ax.set_extent(extent, crs=ccrs.PlateCarree())

# Add state boundaries from shapefile
shapefile = list(shpreader.Reader('shapefile/ne_10m_admin_1_states_provinces.shp').geometries())
ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor='gray', facecolor='none', linewidth=0.3)

# Add coastlines and country borders
ax.coastlines(resolution='10m', color='black', linewidth=0.8)
ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=0.5)

# Add gridlines with labels only on left and bottom
gl = ax.gridlines(crs=ccrs.PlateCarree(), color='gray', linestyle='--', linewidth=0.25,
                  xlocs=np.arange(min_lon, max_lon, grid_interval),
                  ylocs=np.arange(min_lat, max_lat, grid_interval),
                  draw_labels=True)
gl.top_labels = False
gl.right_labels = False

# Define contour levels and corresponding colors for precipitation
levels = [0.7, 6.0, 10, 20, 30, 40, 50, 60, 70, 150, 250, 300]
colors = ["#b4f0f0", "#96d2fa", "#78b9fa", "#3c95f5", "#1e6deb", "#1463d2", 
          "#0fa00f", "#28be28", "#50f050", "#72f06e", "#b3faaa", "#fff9aa"]
cmap = matplotlib.colors.ListedColormap(colors)
cmap.set_over('#294f47')
cmap.set_under('#533e27')
norm = matplotlib.colors.BoundaryNorm(levels, cmap.N)

# Plot accumulated precipitation
cs = plt.contourf(to_np(lons), to_np(lats), accum_rain, levels=levels, norm=norm, cmap=cmap, extend='max')

# Add colorbar with label
cbar = plt.colorbar(cs, ticks=levels, format="%.0f", pad=0.01, fraction=0.04)
cbar.set_label('Precipitation (mm)', fontsize=14)

# Add titles
plt.title('Accumulated Precipitation', fontsize=16, loc='center', y=1.05)
plt.title('WRF - 9 km Resolution', fontsize=12, loc='left')
plt.title(f'{pd.to_datetime(str(times.values)).strftime("%Y-%m-%d %H UTC")}', fontsize=9, loc='right')

# Save the figure
plt.savefig(f'./accumulated_precipitation/{pd.to_datetime(str(times.values)).strftime("%Y%m%d%H")}_rain9km.png', dpi=300)
