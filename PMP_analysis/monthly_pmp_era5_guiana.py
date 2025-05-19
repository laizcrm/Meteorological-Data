##################################################
""" 
Probable Maximum Precipitation (PMP) 
monthly for the Guiana region
based on precipitation data from ERA5 reanalysis.

The World Meteorological Organization (WMO) defines PMP as the theoretically greatest depth
of precipitation for a given duration that is physically possible over a drainage area
at a certain time of the year.

Here, hourly data are used to calculate the moving maximum accumulated precipitation
for different durations (1h, 12h, and 24h), and the maximum values are grouped by month
to obtain the monthly PMP for each duration. This aims to analyze the climatological
behavior of extreme precipitation in the region.
"""



import numpy as np
import matplotlib.pyplot as plt
# %matplotlib inline
from netCDF4 import Dataset
import pandas as pd
import xarray as xr
import matplotlib.colors as colors
import calendar

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from mpl_toolkits.axes_grid1 import make_axes_locatable

import matplotlib
from matplotlib.cm import get_cmap
import matplotlib.colors
import matplotlib as mpl
import matplotlib.cm as cm

import shapely.geometry as sgeom

import shapefile as shp

import datetime as dt
import salem
import glob
# import warnings
# warnings.filterwarnings("ignore")
# from osgeo import gdal # Import the GDAL library

files = glob.glob('/guiana/reanalysis_era5*.nc')
files_list = sorted(files)

# Open multiple NetCDF files as one dataset
ds = xr.open_mfdataset(files_list, combine='by_coords')

# Correct units of precipitation data => m to mm
ds['tp'] = ds.tp * 1000  # [mm]

tp = ds.tp.sel(expver=1)

# Maximum 1-hour precipitation (rolling max over 1 hour)
max_1h = tp.rolling(time=1).max()
pmp_1h = max_1h.groupby('time.month').sum()

# Maximum moving 12-hour precipitation
max_12h = tp.rolling(time=12).max()
pmp_12h = max_12h.groupby('time.month').sum()

# Maximum moving 24-hour precipitation
max_24h = tp.rolling(time=24).max()
pmp_24h = max_24h.groupby('time.month').sum()

pmp_1h.sel(month=1).plot()

# Color palette for precipitation (rain)
rain_colors = [
    "#b4f0f0", "#96d2fa", "#78b9fa", "#3c95f5", "#1e6deb", "#1463d2", 
    "#0fa00f", "#28be28", "#50f050", "#72f06e", "#b3faaa", "#fff9aa", 
    "#ffe978", "#ffc13c", "#ffa200", "#ff6200", "#ff3300", "#ff1500", 
    "#c00100", "#a50200", "#870000", "#653b32"
]
cmap = matplotlib.colors.ListedColormap(rain_colors)
cmap.set_over('#000000')
cmap.set_under('#aaaaaa')

# Map extent: [min_lon, min_lat, max_lon, max_lat]
extent = [-62.0, 1.1, -56.0, 9.0]
min_lon = extent[0]
max_lon = extent[2]
min_lat = extent[1]
max_lat = extent[3]

# Regions with coordinates and names
regions_labels = {
    "R1": (-59.813, 7.6902, "Barima–Waini"),
    "R2": (-60.040, 6.300, "Cuyuni–Mazaruni"),
    "R3": (-57.980, 6.530, "Demerara–Mahaica"),
    "R4": (-57.358, 5.218, "Berbice Oriental"),
    "R5": (-58.200, 6.6114, "Essequibo Islands"),
    "R6": (-57.700, 6.245, "Mahaica–Berbice"),
    "R7": (-58.733, 7.184, "Pomeroon–Supenaam"),
    "R8": (-59.28, 5.00, "Potaro–Siparuni"),
    "R9": (-58.1802, 5.40, "Upper Tacutu"),
    "R10": (-58.96, 2.85, "Upper Demerara"),
}

# Markers and labels for legend
handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='black', markersize=3) for _ in regions_labels]
labels = [f"{region} - {name}" for region, (_, _, name) in regions_labels.items()]

# Create figure with 3 rows and 4 columns, map projection PlateCarree
fig, axs = plt.subplots(3, 4, figsize=(15, 12), subplot_kw={'projection': ccrs.PlateCarree()})
fig.subplots_adjust(hspace=0.2, wspace=0.2)

# Loop through each month and plot PMP 1h
for i, ax in enumerate(axs.flatten()):
    # Select PMP for the month (1 to 12)
    month_data = pmp_1h.sel(month=i + 1)

    # Set map extent
    ax.set_extent([extent[0], extent[2], extent[1], extent[3]], ccrs.PlateCarree())

    # Add ocean feature with gray color
    ax.add_feature(cfeature.OCEAN, facecolor='gray')

    # Add shapefile geometry for Guyana admin boundaries
    shapefile_geoms = list(shpreader.Reader('/guiana/shapefile/guy_admbnd_adm1_2021.shp').geometries())
    ax.add_geometries(shapefile_geoms, ccrs.PlateCarree(), edgecolor='black', facecolor='none', linewidth=0.9)

    # Add coastlines
    ax.coastlines(resolution='10m', color='black', linewidth=0.3)

    # Add gridlines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), color='black', alpha=0.3, linestyle='--', linewidth=0.25,
                      xlocs=np.arange(-180, 180, 1), ylocs=np.arange(-90, 90, 1), draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False

    # Plot contourf of precipitation data
    cont = ax.contourf(month_data.longitude, month_data.latitude, month_data, transform=ccrs.PlateCarree(),
                       cmap=cmap, extend='both', alpha=0.6)

    # Month abbreviation as title
    month_name = calendar.month_abbr[i + 1]
    ax.set_title(f'{month_name}')

    # Plot region labels
    for region, (lon, lat, name) in regions_labels.items():
        ax.text(lon, lat, f"{region}", color='red', size=6, ha='right', va='bottom', transform=ccrs.PlateCarree())

# Add legend and colorbar
legend = fig.legend(handles, labels, loc='center right', bbox_to_anchor=(0.12, 0.2), title="Regions", fontsize=9)
plt.setp(legend.get_title(), fontsize=9)

cbar_ax = fig.add_axes([0.2, 0.05, 0.6, 0.02])
cb = plt.colorbar(cont, cax=cbar_ax, orientation='horizontal', label='Precipitation (mm)')

# Overall figure title
plt.suptitle('Monthly Precipitation - Maximum Accumulated Precipitation of 1H', fontsize=16, y=1.02)

# Save figure
plt.savefig('climatology_precipitation_guiana_1h.png')

plt.show()
