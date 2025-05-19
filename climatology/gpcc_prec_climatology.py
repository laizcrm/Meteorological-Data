"""
Data source: GPCC (Global Precipitation Climatology Centre)
URL: https://psl.noaa.gov/data/gridded/data.gpcc.html

Description:
The GPCC dataset provides monthly total precipitation on a 2.5° x 2.5° global grid,
derived from rain gauge measurements collected worldwide. It offers long-term climate
records of precipitation and is widely used for climate monitoring, model validation,
and hydrological studies. The version used here covers monthly precipitation totals
aggregated from station observations and interpolated to a regular grid, representing
climatological averages for each month.
"""


import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import numpy as np

# Shapefile for Brazilian states boundaries
shapefile_path = '/home/fieldpro/shapefile/BR_UF_2019.shp'

# Map extent [min_lon, max_lon, min_lat, max_lat]
extent = [-80.0, -40.0, -30.0, 10.0]

lon = file.variables['lon'][:]
lat = file.variables['lat'][:]
xx, yy = np.meshgrid(lon, lat)

# Create figure and axis with PlateCarree projection
fig, axs = plt.subplots(3, 4, figsize=(15, 12), subplot_kw={'projection': ccrs.PlateCarree()})
fig.subplots_adjust(hspace=0.2, wspace=0.2)

# Loop through each month and plot climatology
for i, ax in enumerate(axs.flatten()):
    month_data = clima_tp.sel(month=i + 1)
    
    # Set map extent
    ax.set_extent([extent[0], extent[1], extent[2], extent[3]], ccrs.PlateCarree())
    
    # Add shapefile geometry (Brazil states borders)
    shapefile = list(shpreader.Reader(shapefile_path).geometries())
    ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor='gray', facecolor='none', linewidth=0.3)
    
    # Add coastlines and borders
    ax.coastlines(resolution='10m', color='black', linewidth=0.8)
    ax.add_feature(cartopy.feature.BORDERS, edgecolor='black', linewidth=0.5)
    
    # Add gridlines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), color='gray', alpha=1.0, linestyle='--', linewidth=0.25,
                      xlocs=np.arange(-180, 180, 5), ylocs=np.arange(-90, 90, 5), draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False
    
    # Plot precipitation monthly mean using contourf
    contour = ax.contourf(xx, yy, month_data, 
                          transform=ccrs.PlateCarree(),
                          cmap='rainbow',
                          extend='max',
                          alpha=0.6)
    
    # Add month title on each subplot
    ax.set_title(calendar.month_name[i+1], fontsize=10)

# Add a single colorbar below all plots
cbar_ax = fig.add_axes([0.2, 0.05, 0.6, 0.02])
cbar = plt.colorbar(contour, cax=cbar_ax, spacing='proportional', 
                    orientation='horizontal', label='Precipitation (mm)', format='%.1f')

plt.suptitle('Monthly Mean Precipitation Climatology - GPCC', fontsize=16, y=1.02)

plt.savefig('monthly_precip_climatology_GPCC.png', bbox_inches='tight', dpi=300)
plt.show()
