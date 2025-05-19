"""
WRF Output Visualization Script

This script processes WRF model output files for a specific date/time range,
extracting sea level pressure (SLP) and wind components at 1000 hPa.
It smooths the SLP field, interpolates wind to the 1000 hPa level,
and plots these variables over a specified geographic region
using Cartopy for map projection and visualization.

The output is saved as PNG images showing wind vectors and sea level pressure contours,
with geographic features such as coastlines, borders, and state boundaries for context.
"""

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from netCDF4 import Dataset
import wrf
from wrf import getvar, latlon_coords, get_cartopy, to_np, smooth2d, interplevel
import pandas as pd
from glob import glob
from scipy.ndimage.filters import gaussian_filter
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader

# Load and sort WRF output files
files = sorted(glob("run2/wrfout_d01_*"))

# Define datetime strings to match files
dates = list(pd.date_range('2016-07-13T00:00:00', '2016-07-13T00:00:00', freq='1H').strftime('%Y-%m-%d_%H:%M:%S'))

# Select matching files based on datetime strings
selected_files = [f for f in files if f[16:] in dates]

for file in selected_files:
    ds = Dataset(file)
    time_var = getvar(ds, "times")
    
    # Extract date/time components
    year = time_var.dt.year
    month = time_var.dt.month
    day = time_var.dt.day
    hour = time_var.dt.hour
    
    # Get variables
    slp = getvar(ds, "slp")
    pressure = getvar(ds, "pressure")
    u = getvar(ds, "ua", units="m/s")
    v = getvar(ds, "va", units="m/s")
    
    # Interpolate winds to 1000 hPa
    ua_1000 = interplevel(u, pressure, 1000)
    va_1000 = interplevel(v, pressure, 1000)
    
    # Smooth sea level pressure
    slp_smooth = smooth2d(slp, 10, cenweight=4)
    slp_smooth = gaussian_filter(slp_smooth, sigma=3)
    
    # Get lat/lon coordinates and projection
    lats, lons = latlon_coords(ua_1000)
    projection = get_cartopy(ua_1000)
    
    # Define geographic extent [min_lon, min_lat, max_lon, max_lat]
    extent = [-80, -38, -41, -6]
    
    # Create figure and axis with PlateCarree projection
    fig = plt.figure(figsize=(10,10))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    
    # Add background land with transparency
    ax.add_feature(cfeature.LAND, facecolor='gray', alpha=0.4)
    
    # Add state boundaries from shapefile
    shapefile = list(shpreader.Reader('shapefile/ne_10m_admin_1_states_provinces.shp').geometries())
    ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor='gray', facecolor='none', linewidth=0.3)
    
    # Add coastlines and borders
    ax.coastlines(resolution='10m', color='gray', linewidth=0.8)
    ax.add_feature(cfeature.BORDERS, edgecolor='gray', linewidth=0.5)
    
    # Add gridlines with labels
    gl = ax.gridlines(crs=ccrs.PlateCarree(), color='gray', linestyle='--', linewidth=0.25,
                      xlocs=np.arange(extent[0], extent[2], 5),
                      ylocs=np.arange(extent[1], extent[3], 5),
                      draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False
    
    # Contour sea level pressure
    contour_levels = np.arange(np.min(slp_smooth), np.max(slp_smooth), 2)
    slp_contours = ax.contour(to_np(lons), to_np(lats), to_np(slp_smooth), levels=contour_levels,
                              cmap='jet', linewidths=1.6, alpha=0.6, transform=ccrs.PlateCarree())
    plt.clabel(slp_contours, inline=True, fontsize=12, fmt='%1.0f', colors='darkred')
    
    # Plot wind vectors (every 10th point)
    skip = 10
    quiver = ax.quiver(to_np(lons[::skip, ::skip]), to_np(lats[::skip, ::skip]),
                       to_np(ua_1000[::skip, ::skip]), to_np(va_1000[::skip, ::skip]),
                       color='black', transform=ccrs.PlateCarree())
    ax.quiverkey(quiver, 0.8, 0.12, 12, '12 m/s', labelpos='E', coordinates='figure')
    
    # Titles
    plt.title('Wind (m/s) and Sea Level Pressure (hPa)', fontsize=14, y=1.05)
    plt.title('WRF - 9 km', loc='left', fontsize=12)
    plt.title(str(pd.to_datetime(str(time_var.values))) + ' UTC', fontsize=9, loc='right')
    
    # Save figure
    output_filename = f"./surface/{pd.to_datetime(str(time_var.values)).strftime('%Y-%m-%d_%H')}pnmm.png"
    plt.savefig(output_filename, dpi=300)
    plt.close(fig)
    
    print('---------------------------------------')
    print('Script executed for:', time_var.values)
