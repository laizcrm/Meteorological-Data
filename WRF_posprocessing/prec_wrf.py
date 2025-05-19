"""
This script reads output data from the WRF model (NetCDF file), calculates 24-hour accumulated precipitation 
from the RAINC and RAINNC variables, and generates contour maps of this precipitation using the Cartopy library. 
For each 24-hour interval, the script plots and saves a georeferenced figure showing the accumulated rainfall, 
including state boundaries and a colorbar representing values in millimeters. The PlateCarree projection is used, 
and the script loops over all available time steps in the dataset.
"""


import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from datetime import datetime
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Open WRF output file
f = Dataset("wrfout_d01_2019-01-20_00:00:00.nc")

# Calculate precipitation
rain = f.variables["RAINC"][:] + f.variables["RAINNC"][:]
rain_1h = np.zeros(rain.shape)
rain_24h = np.zeros(rain.shape)

rain_1h = rain[1:, :, :] - rain[:-1, :, :]
rain_24h = rain[24:, :, :] - rain[:-24, :, :]

# Location variables (do not change with time)
lat = f.variables["XLAT"][0]
lon = f.variables["XLONG"][0]

# Loop through available time indices (one every 24 hours)
for idx_time in range(24, len(rain_24h[:, 0, 0]), 24):
    
    # Extract and decode time information
    time_str = f.variables["Times"][idx_time].tobytes().decode("UTF-8").replace("\x00", "")
    time_obj = datetime.strptime(time_str, "%Y-%m-%d_%H:%M:%S")
    print(time_obj)

    # Set up projection
    crs = ccrs.PlateCarree()
    fig = plt.figure(dpi=300, facecolor='none', edgecolor='k')
    ax = plt.axes(projection=crs, facecolor='#f0edf0ff')

    # Plot title
    title = "24h Accumulated Precipitation (mm) on {}".format(time_obj)
    plt.title(title, fontsize=9)

    # Define map extent based on lat/lon limits
    ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()])

    # Add state/province borders
    states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces',
        scale='50m',
        facecolor="none"
    )
    ax.add_feature(states_provinces, edgecolor='black', linewidth=0.5, alpha=0.5)

    # Define contour levels and plot precipitation
    clevels = [0.7, 4, 6, 8, 10, 15, 20, 25, 30,
               35, 40, 45, 50, 55, 60, 65, 70, 80]
    p = plt.contourf(lon, lat, rain_24h[idx_time], clevels, cmap=plt.get_cmap('terrain'), transform=crs)

    # Add colorbar
    cbar = plt.colorbar(p, ax=ax, orientation="horizontal")
    cbar.ax.tick_params(labelsize=5)

    # Save figure
    plt.savefig("{:%Y%m%d%H}_rain24h.png".format(time_obj))
