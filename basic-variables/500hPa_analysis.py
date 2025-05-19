#------- INSTALL BASEMAP AND DEPENDENCIES -----------------------

# At the time this script was written, Cartopy had compatibility issues
# with shapefiles and certain map projections. As a workaround, Basemap was used
# despite being deprecated. Today, Cartopy is preferred for new projects.
#!sudo apt-get install python-grib 
#----------------------------------------------------------------
#       NOTES
#!sudo pip3 install -U git+https://github.com/matplotlib/basemap.git
#https://jswhit.github.io/pygrib/docs/--- GRIB files
#-------------------------------------------------------------------------

#Instalar o Pyproj
!pip install -q pyproj==2.3.0

#Install Basemap 
!apt-get install -q libgeos-dev

!pip install -q https://github.com/matplotlib/basemap/archive/master.zip

!pip install netCDF4

!pip install cartopy

!pip install metpy==1.0


# ----- MAP SETUP AND RESOURCES ------------------------------------
!pip install cartopy
!pip install shapely --no-binary shapely --force
print('\n')

# Download the shapefile of Brazilian states from IBGE
!wget -c https://geoftp.ibge.gov.br/organizacao_do_territorio/malhas_territoriais/malhas_municipais/municipio_2019/Brasil/BR/br_unidades_da_federacao.zip
print('\n')

# Download shapefile of world states/provinces from Natural Earth
!wget -c https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/cultural/ne_10m_admin_1_states_provinces.zip
print('\n')

!unzip -o ne_10m_admin_1_states_provinces.zip
print('\n')

#-------------------- LIBRARIES --------------------
# Core libraries
import os
import glob
import numpy as np
import pandas as pd
import datetime as dt

# Plotting
import matplotlib.pyplot as plt
%matplotlib inline
import matplotlib.colors as colors
import matplotlib.cm as cm


# NetCDF and data handling
from netCDF4 import Dataset, num2date
import xarray as xr

# Meteorological tools
import metpy.calc as mpcalc
from metpy.units import units
from metpy.interpolate import cross_section

# Remote data access
from siphon.ncss import NCSS

# Mapping
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
#from osgeo import gdal


# Load the NetCDF file containing 500 hPa atmospheric data 
file = xr.open_dataset("/content/nivel_500.nc")
print(file)


# Time range
# timestamps every 6 hours between July 13 and July 15, 2016
date = list(pd.date_range('2016-07-13T00:00:00.00','2016-07-15T00:00:00.00', freq='6H').strftime('%Y-%m-%d %H'))

# Define map extent 
extent = [-95.0,-60.0,-20.0,10.0] #(west_lon, east_lon, south_lat, north_lat)
min_lon = extent[0]; max_lon = extent[2]; min_lat = extent[1]; max_lat = extent[3]

#----------Vorticity (500hPa)--------------------------------
for time_idx in range(len(date)):
    # Select data for the current time step
    ds = file.sel(time=date[time_idx]).metpy.parse_cf()

    # Extract variables
    u = ds['u'].values
    v = ds['v'].values
    vo = ds['vo'].values
    vo_c = vo * 1e5  # Convert vorticity units to 10^-5 s^-1 scale
    z = ds['z'].values

    # Convert wind components to meters per second with units
    u_wind_m = units('m/s') * u.squeeze()
    v_wind_m = units('m/s') * v.squeeze()

    # Get latitude and longitude arrays
    lat = ds['latitude'].values.squeeze()
    lon = ds['longitude'].values.squeeze()

    # Create 2D grid of lat/lon for plotting
    lon_2d, lat_2d = np.meshgrid(lon, lat)

    # Smooth geopotential height with Gaussian filter
    Z_500 = ndimage.gaussian_filter(z, sigma=3, order=0) * units.meter

    #-----------------Figure--------------
    fig, ax = plt.subplots(figsize=(12, 9))

    # Set up Basemap with rectangular projection and defined extent
    bmap = Basemap(llcrnrlon=extent[0], llcrnrlat=extent[1],
                   urcrnrlon=extent[2], urcrnrlat=extent[3],
                   epsg=4326, resolution='i', ax=ax)


    bmap.readshapefile('/content/BR_UF_2019', 'BRA_UF_2019', linewidth=0.3, color='black')
    bmap.readshapefile('/content/South_America', 'South_America', linewidth=0.3, color='black')


    bmap.drawparallels(np.arange(-90, 90, 5), linewidth=0.2, dashes=[4, 4], color='white',
                       labels=[True, False, False, True], fmt='%g', labelstyle="+/-", size=12)
    bmap.drawmeridians(np.arange(0, 360, 5), linewidth=0.2, dashes=[4, 4], color='white',
                       labels=[True, False, False, True], fmt='%g', labelstyle="+/-", size=12)

    # Prepare latitude and longitude arrays for contour plotting (rectangular grid)
    x = np.linspace(extent[0], extent[2], vo_c.shape[1])
    y = np.linspace(extent[3], extent[1], vo_c.shape[0])  # note reversed order for latitudes
    xx, yy = np.meshgrid(x, y)

    # Define contour levels for vorticity
    contour_levels = np.arange(-10, 0, 1)

    # vorticity relative 
    vort_plot = ax.contourf(xx, yy, vo_c, cmap='Blues_r', levels=contour_levels, extend='min')

    # Add colorbar for vorticity
    cbar = plt.colorbar(vort_plot, ax=ax, ticks=contour_levels, format='%.0f')
    cbar.set_label('Relative Vorticity ($10^{-5}$ s$^{-1}$)')

    # Plot wind barbs (every 15th point to avoid clutter)
    ax.barbs(lon_2d[::15, ::15], lat_2d[::15, ::15],
             u[::15, ::15], v[::15, ::15],
             length=5, fill_empty=False,
             sizes=dict(emptybarb=0., spacing=0.2, height=0.3))

    # Plot geopotential height contours (in decameters)
    geo_contours = ax.contour(xx, yy, Z_500.magnitude / 10, colors='darkred',
                              linestyles='dashed', linewidths=1.5)
    ax.clabel(geo_contours, inline=True, fontsize=10, fmt='%1.0f', colors='darkred')

    # Titles
    ax.set_title('Relative Vorticity ($10^{-5}$ s$^{-1}$) and Geopotential Height (dam)',
                 fontsize=12, y=1.05)
    ax.text(0.01, 0.97, 'ERA-Interim 0.25° x 0.25°', transform=ax.transAxes, fontsize=9, ha='left')
    ax.text(0.99, 0.97, pd.to_datetime(str(ds.time.values)).strftime('%Y-%m-%d %H UTC'),
            transform=ax.transAxes, fontsize=9, ha='right')

    # Save figure with timestamp in filename
    timestamp = pd.to_datetime(str(ds.time.values)).strftime('%Y%m%d_%H%M')
    plt.savefig(f'{timestamp}_500_vorticity.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

#---------------------Relative Vorticity Advection(500 hPa)-----------------------------------
for time_idx in range(len(date)):
    # Select dataset for current time
    ds = file.sel(time=date[time_idx]).metpy.parse_cf()

    # Extract variables
    u = ds['u'].values
    v = ds['v'].values
    vo = ds['vo'].values
    z = ds['z'].values

    # Convert winds to knots for plotting wind barbs
    u_wind_knots = u * units.knots
    v_wind_knots = v * units.knots

    # Create 1D latitude and longitude arrays
    lat = ds['latitude'].values.squeeze()
    lon = ds['longitude'].values.squeeze()

    # Create 2D meshgrid of lon and lat for plotting
    lon_2d, lat_2d = np.meshgrid(lon, lat)

    # Calculate grid spacing (meters) for advection calculation
    dx, dy = mpcalc.lat_lon_grid_deltas(lon, lat)

    # Calculate vorticity advection (units: s^-2), scale by 1e9 for better visualization
    adv = mpcalc.advection(vo, u, v) * 1e9

    # Smooth geopotential height and advection fields with Gaussian filter
    Z_500 = ndimage.gaussian_filter(z, sigma=3) * units.meter
    adv = ndimage.gaussian_filter(adv, sigma=3)
    #--------------------Figure----------------------------------------
    fig, ax = plt.subplots(figsize=(12, 9))

    bmap = Basemap(llcrnrlon=extent[0], llcrnrlat=extent[1],
                   urcrnrlon=extent[2], urcrnrlat=extent[3],
                   epsg=4326, resolution='i', ax=ax)

  
    bmap.readshapefile('/content/BR_UF_2019', 'BRA_UF_2019', linewidth=0.3, color='black')
    bmap.readshapefile('/content/South_America', 'South_America', linewidth=0.3, color='black')


    x = np.linspace(extent[0], extent[2], adv.shape[1])
    y = np.linspace(extent[3], extent[1], adv.shape[0])  # note reversed latitude order
    xx, yy = np.meshgrid(x, y)

    # Define contour levels for advection, centered around zero
    levels = np.arange(-10, 10, 1)

    # Prepare custom colormap: blue-white-red with white near zero
    ncolors = list(plt.cm.bwr_r(np.linspace(0, 1, len(levels))))
    ncolors[9] = "white"  # Make level near zero white
    ncolors[10] = "white"
    cmap = colors.ListedColormap(ncolors, name="custom_bwr", N=len(ncolors))
    norm = colors.BoundaryNorm(levels, ncolors=len(ncolors))

    # Contours and filled contours
    cf = ax.contourf(xx, yy, adv, cmap=cmap, levels=levels, norm=norm, extend="both")
    plt.colorbar(cf, ax=ax, ticks=np.arange(-10, 10, 1), format='%.1f', label='Vorticity Advection ($10^{-9} s^{-2}$)')

    # Plot wind barbs every 15th point to avoid clutter
    w = 15
    ax.barbs(lon_2d[::w, ::w], lat_2d[::w, ::w],
             u_wind_knots[::w, ::w], v_wind_knots[::w, ::w],
             length=5, fill_empty=False, sizes=dict(emptybarb=0., spacing=0.2, height=0.3))

    # Plot geopotential height contours (in decameters)
    cs = ax.contour(xx, yy, Z_500.magnitude / 10, colors='black', linestyles='dashed', linewidths=1.5)
    ax.clabel(cs, inline=True, fontsize=10, fmt='%1.0f')

    # Titles and annotations
    ax.set_title('Relative Vorticity Advection ($10^{-9}$ $s^{-2}$) and Geopotential Height (dam)',
                 fontsize=12, y=1.05)
    ax.text(0.01, 0.97, 'ERA-Interim 0.25° x 0.25°', transform=ax.transAxes, fontsize=9, ha='left')
    time_str = pd.to_datetime(str(ds.time.values)).strftime('%Y-%m-%d %H UTC')
    ax.text(0.99, 0.97, time_str, transform=ax.transAxes, fontsize=9, ha='right')

    # Save figure with timestamp
    filename = f"{pd.to_datetime(str(ds.time.values)).strftime('%Y%m%d_%H%M')}_vort_adv.png"
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close(fig)
