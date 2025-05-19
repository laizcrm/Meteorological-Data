"""
Script to calculate and visualize vertically integrated moisture flux between 900-700 hPa from WRF model output.

This script reads multiple WRF output files, extracts meteorological variables such as temperature, dew point,
wind components, and pressure levels. It computes the specific humidity and moisture flux vectors at several pressure
levels between 900 and 700 hPa, integrates them vertically, and plots the resulting moisture flux magnitude and vectors
over a specified geographic region in South America. The plot includes administrative boundaries, coastlines, gridlines,
and a colorbar indicating the magnitude of moisture flux.
"""

import datetime 
import scipy.ndimage as ndimage
from siphon.ncss import NCSS
from scipy.ndimage.filters import gaussian_filter
import scipy.integrate as integrate


import numpy as np
from numpy import linspace, meshgrid


import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from netCDF4 import Dataset
import pandas as pd
import xarray as xr

from glob import glob
import os


import cartopy.crs as crs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader


import metpy.calc as mpcalc
from metpy.units import units


import wrf as wrf
from wrf import getvar, latlon_coords, get_cartopy, to_np, interplevel

# -----------------------
# File and time selection
# -----------------------

# List all WRF output files
wrf_files = sorted(glob("run2/wrfout_d01_*"))

# Define desired datetime range and format
date_range = list(pd.date_range('2016-07-13T00:00:00', '2016-07-15T12:00:00', freq='1H').strftime('%Y-%m-%d_%H:%M:%S'))

# Filter files matching the date range exactly
selected_files = []
for idx in range(len(date_range)):
    filename_date = wrf_files[idx][16:]
    if filename_date == date_range[idx]:
        selected_files.append(wrf_files[idx])

# -----------------------
# Process each selected file
# -----------------------
for file_path in selected_files:

    ds = Dataset(file_path)

    # Extract time variable
    time_var = getvar(ds, "times")
    
    # Extract meteorological variables
    pressure = getvar(ds, "pressure")
    temp_c = getvar(ds, "tc")        # Temperature in Celsius
    dewpoint_c = getvar(ds, 'td', units='degC')
    u_wind = getvar(ds, "ua", units="m/s")
    v_wind = getvar(ds, "va", units="m/s")
    rh = getvar(ds, 'rh2')
    surface_pressure = getvar(ds, "pressure")

    # Interpolate temperature and dewpoint at specific pressure levels
    pressure_levels = [900, 850, 800, 750, 700]
    temp_levels = interplevel(temp_c, pressure, pressure_levels)
    dewpoint_levels = interplevel(dewpoint_c, pressure, pressure_levels)
    u_levels = interplevel(u_wind, pressure, pressure_levels)
    v_levels = interplevel(v_wind, pressure, pressure_levels)
    pz = interplevel(surface_pressure, pressure, pressure_levels)

    # Constants
    g = 9.8          # Gravity acceleration (m/s^2)
    epsilon = 0.622  # Ratio of molecular weights of water vapor/dry air

    # Calculate vapor pressure e (hPa)
    e = 6.11 * np.exp((17.67 * dewpoint_levels) / (dewpoint_levels + 243.5))

    # Calculate specific humidity q (g/kg)
    q = (epsilon * e) / (pz + e * (epsilon - 1)) * 1e3

    # Calculate moisture flux components
    qx = q * u_levels
    qy = q * v_levels

    # Vertically integrate moisture flux (average over levels)
    q_vec_u = np.sum(qx, axis=0) / len(pressure_levels)
    q_vec_v = np.sum(qy, axis=0) / len(pressure_levels)

    # Magnitude of moisture flux vector
    q_mag = np.sqrt(q_vec_u**2 + q_vec_v**2)

    # Get lat/lon coordinates
    lats, lons = latlon_coords(u_levels)

    # Set map extent
    extent = [-80, -41, -38, -6]#(west_lon, east_lon, south_lat, north_lat)
    interval = 5  # Gridline interval degrees

    fig = plt.figure(figsize=(10, 10))
    ax = plt.axes(projection=crs.PlateCarree())
    ax.set_extent([extent[0], extent[1], extent[2], extent[3]], crs=crs.PlateCarree())

    # Add background land with transparency
    ax.add_feature(cfeature.LAND, facecolor='gray', alpha=0.4)

    # Add shapefile geometries
    shapefile = list(shpreader.Reader('shapefile/ne_10m_admin_1_states_provinces.shp').geometries())
    ax.add_geometries(shapefile, crs.PlateCarree(), edgecolor='gray', facecolor='none', linewidth=0.3)

    # Add coastlines, borders and gridlines
    ax.coastlines(resolution='10m', color='black', linewidth=0.8)
    ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=0.5)
    gl = ax.gridlines(crs=crs.PlateCarree(), color='gray', linestyle='--', linewidth=0.25,
                      xlocs=np.arange(extent[0], extent[1], interval),
                      ylocs=np.arange(extent[2], extent[3], interval),
                      draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False

    # Create contour levels for plotting
    levels = np.arange(100, 550, 50)

    # Define color map
    colors_list = ['#66aad4', '#0072b8', '#a9d281', '#6fb42d', '#f7ba73',
                   '#f28c16', '#f070b3', '#e61082']
    cmap = mpl.colors.ListedColormap(colors_list)
    cmap.set_over('#b80d68')

    # Plot filled contours for moisture flux magnitude
    contourf = ax.contourf(to_np(lons), to_np(lats), to_np(q_mag), levels=levels,
                          cmap=cmap, extend='max', transform=crs.PlateCarree())

    # Add colorbar
    cbar = plt.colorbar(contourf, ticks=levels, format='%.0f', pad=0.01, fraction=0.04)
    cbar.set_label('Moisture Flux Magnitude (g m$^{-2}$ s$^{-1}$)', fontsize=14)

    # Plot moisture flux vectors (quiver)
    step = 15  # Vector spacing
    q = ax.quiver(to_np(lons[::step, ::step]), to_np(lats[::step, ::step]),
                  to_np(q_vec_u[::step, ::step]), to_np(q_vec_v[::step, ::step]),
                  color='black', transform=crs.PlateCarree())

    ax.quiverkey(q, 0.80, 0.12, 100, '100 g m$^{-2}$ s$^{-1}$', labelpos='E', coordinates='figure')

    plt.title('Vertically Integrated Moisture Flux between 900-700 hPa', fontsize=16, y=1.05)
    plt.title(f'WRF Output Time: {pd.to_datetime(str(time_var.values))} UTC', fontsize=9, loc='right')
    plt.title('WRF - 9 km resolution', fontsize=12, loc='left')

    output_dir = './850hpa/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    output_filename = f"{output_dir}{pd.to_datetime(str(time_var.values)).strftime('%Y-%m-%d_%H')}_moisture_flux.png"
    plt.savefig(output_filename, dpi=300)
    plt.close()

    print('---------------------------------------')
    print('Processed file:', file_path)
    print('Output saved as:', output_filename)

