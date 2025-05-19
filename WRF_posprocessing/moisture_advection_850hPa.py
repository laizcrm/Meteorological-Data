"""
Script to calculate and plot moisture advection at 850 hPa from WRF output files.

- Loads multiple WRF output NetCDF files.
- Extracts meteorological variables (temperature, wind components, precipitable water, geopotential height).
- Interpolates variables to 850 hPa level.
- Calculates moisture advection using MetPy.
- Plots precipitable water with state boundaries, coastlines, and gridlines using Cartopy.
- Saves the figure with timestamp in filename.
"""

import os
from glob import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as crs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
import wrf
from wrf import getvar, interplevel, smooth2d, latlon_coords, to_np, get_cartopy
import metpy.calc as mpcalc

# Define input files and sort
wrf_files = sorted(glob("run2/wrfout_d01_*"))

# Define expected timestamps for filtering files
expected_times = pd.date_range('2016-07-13T00:00:00', '2016-07-15T12:00:00', freq='1H').strftime('%Y-%m-%d_%H:%M:%S').tolist()

# Select files matching expected timestamps
selected_files = []
for f in wrf_files:
    file_time = os.path.basename(f)[16:]
    if file_time in expected_times:
        selected_files.append(f)

# Geographic extent for plotting [lon_min, lon_max, lat_min, lat_max]
extent = [-80, -41, -38, -6]
contour_interval = 5  # For precipitable water shading

for file in selected_files:
    with wrf.Dataset(file) as ds:
        # Extract time variable and convert to pandas datetime
        time = getvar(ds, "times")
        dt = pd.to_datetime(str(time.values))

        # Extract meteorological variables
        slp = getvar(ds, "slp")
        p = getvar(ds, "pressure")
        pw = getvar(ds, 'pw')  # Precipitable water
        z = getvar(ds, "z", units="dm")
        tc = getvar(ds, "tc")  # Temperature (Celsius)
        u = getvar(ds, "ua", units="m/s")
        v = getvar(ds, "va", units="m/s")
        wspd_dir = getvar(ds, "wspd_wdir", units="m/s")[0, :]

        # Interpolate variables to 850 hPa level
        temp_850 = interplevel(tc, p, 850)
        ua_850 = interplevel(u, p, 850)
        va_850 = interplevel(v, p, 850)
        hgt_850 = interplevel(z, p, 850)
        wspd_850 = interplevel(wspd_dir, p, 850)

        # Smooth height field for plotting
        hgt_smoothed = smooth2d(hgt_850, 1)

        # Get lat/lon coordinates
        lats, lons = latlon_coords(wspd_850)

        # Calculate grid spacing in meters for advection calculation
        dx, dy = mpcalc.lat_lon_grid_deltas(lons.values, lats.values)

        # Calculate moisture advection (units scaled by 1e4 for convenience)
        adv = mpcalc.advection(temp_850.values, ua_850.values, va_850.values, dx=dx, dy=dy) * 1e4

        # Setup Cartopy projection
        cart_proj = get_cartopy(wspd_850)

        # Plotting setup
        fig = plt.figure(figsize=(10, 10))
        ax = plt.axes(projection=crs.PlateCarree())
        ax.set_extent([extent[0], extent[1], extent[2], extent[3]], crs.PlateCarree())

        # Add state boundaries from shapefile
        shapefile = list(shpreader.Reader('shapefile/ne_10m_admin_1_states_provinces.shp').geometries())
        ax.add_geometries(shapefile, crs.PlateCarree(), edgecolor='gray', facecolor='none', linewidth=0.3)

        # Add coastlines and borders
        ax.coastlines(resolution='10m', color='gray', linewidth=0.8)
        ax.add_feature(cfeature.BORDERS, edgecolor='gray', linewidth=0.5)

        # Add gridlines
        gl = ax.gridlines(crs=crs.PlateCarree(), color='gray', alpha=1.0, linestyle='--', linewidth=0.25,
                          xlocs=np.arange(extent[0], extent[1], 5),
                          ylocs=np.arange(extent[2], extent[3], 5),
                          draw_labels=True)
        gl.top_labels = False
        gl.right_labels = False

        # Prepare contour levels for precipitable water
        pw_min = np.min(pw)
        pw_max = np.max(pw)
        levels = np.arange(pw_min, pw_max, contour_interval)

        # Plot precipitable water shaded contours
        cf = ax.contourf(to_np(lons), to_np(lats), to_np(pw), levels=levels, cmap='Greens', extend='max', transform=crs.PlateCarree())
        cb = plt.colorbar(cf, ticks=levels, format='%.0f', pad=0.01, fraction=0.04)
        cb.set_label('mm')

        # Titles
        plt.title('Precipitable Water', loc='center', fontsize=14, y=1.05)
        plt.title('WRF Model - 9 km Resolution', loc='left', fontsize=12)
        plt.title(f'{dt.strftime("%Y-%m-%d %H UTC")}', loc='right', fontsize=9)

        # Save figure
        plt.savefig(f'pw_9km_{dt.strftime("%Y%m%d_%H%M")}.png', dpi=300)
        plt.close(fig)

        print(f'Processed and saved figure for time: {dt}')
