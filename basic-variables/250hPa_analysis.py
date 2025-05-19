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

# Pyproj
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


# Load the NetCDF file containing 250 hPa atmospheric data 
file = xr.open_dataset("/content/nivel_250.nc")
print(file)

# Time range
# timestamps every 6 hours between July 13 and July 15, 2016
date = list(pd.date_range('2016-07-13T00:00:00.00','2016-07-15T00:00:00.00', freq='6H').strftime('%Y-%m-%d %H'))

# Define map extent
extent = [-95.0,-60.0,-20.0,10.0] # (west_lon, east_lon, south_lat, north_lat)
min_lon = extent[0]; max_lon = extent[2]; min_lat = extent[1]; max_lat = extent[3]

#----------Plotting upper-level Jet Stream data (250 hPa)--------------------------
for idx_time in range(len(date)):
    # Select the dataset for the current time step, parsing CF metadata
    ds = file.sel(time=date[idx_time]).metpy.parse_cf()

    # Extract wind components at 250 hPa
    u_250 = ds['u'][:].values
    v_250 = ds['v'][:].values

    # Extract geopotential height variable
    z = ds['z']

    # Convert wind components to knots using MetPy units
    u_wind_knots = ds['u'].values * units.knots
    v_wind_knots = ds['v'].values * units.knots

    # Also convert wind components to meters per second
    u_wind_mps = units('m/s') * ds.variables['u'][:].squeeze()
    v_wind_mps = units('m/s') * ds.variables['v'][:].squeeze()

    # Extract latitude and longitude arrays
    lat = ds.variables['latitude'][:].squeeze()
    lon = ds.variables['longitude'][:].squeeze()

    # Calculate wind speed magnitude
    wind_speed = np.sqrt(u_250**2 + v_250**2)

    # Create 2D meshgrid of lon/lat coordinates
    lon_2d, lat_2d = np.meshgrid(lon, lat)

    # Smooth geopotential height with Gaussian filter (sigma=3)
    Z_250_smoothed = ndimage.gaussian_filter(z, sigma=3, order=0) * units.meter

    # -------------------Figure-------------------------------------------
    fig, ax = plt.subplots(figsize=(12, 9))

    bmap = Basemap(llcrnrlon=extent[0], llcrnrlat=extent[1],
                   urcrnrlon=extent[2], urcrnrlat=extent[3],
                   epsg=4326, resolution='i')


    bmap.readshapefile('/content/BR_UF_2019', 'BRA_UF_2019',
                       linewidth=0.3, color='black')
    bmap.readshapefile('/content/South_America', 'South_America',
                       linewidth=0.3, color='black')

    # Define evenly spaced grid points for plotting wind vectors
    x = np.linspace(min_lon, max_lon, u_250.shape[1])
    y = np.linspace(max_lat, min_lat, u_250.shape[0])
    
  
  # Create 2D grid from lat/lon
    xx, yy = np.meshgrid(x, y)
    

    levels = np.arange(40, np.max(mag_v), 5)
    
    # Filled contour
    contour = plt.contourf(xx, yy, mag_v, cmap='viridis_r', levels=levels, extend='max')
    
    cbar = plt.colorbar(contour, ticks=levels, format='%.1f')
    cbar.set_label('Wind Speed (m $s^{-1}$)', fontsize=14)
    
    # Plot wind barbs (vectors)
    plt.barbs(lon_2d[::14, ::14], lat_2d[::14, ::14], 
              u_wind[::14, ::14], v_wind[::14, ::14], 
              length=6, fill_empty=False, sizes=dict(emptybarb=0., spacing=0.2, height=0.3))
    
    # Geopotential height contours (in dam)
    geo_contours = plt.contour(xx, yy, Z_250 / 10, colors='darkred', 
                               linestyles='dashed', linewidths=1.5, alpha=0.9)
    plt.clabel(geo_contours, inline=True, fontsize=10, fmt='%1.0f', colors='darkred')
    
    # Titles
    plt.title('Wind Speed Magnitude (m $s^{-1}$) & Geopotential Height (dam) at 250 hPa', fontsize=12, y=1.05)
    plt.title('ERA', loc='left', fontsize=9)
    plt.title(pd.to_datetime(str(ds.time.values)).strftime('%Y-%m-%d %H UTC'), loc='right', fontsize=9)
    
    # Save figure with formatted filename
    plt.savefig('{}_250_wind.png'.format(pd.to_datetime(str(ds.time.values)).strftime('%Y-%m-%d_%H')), dpi=300)
    plt.close(fig)

#----------Plotting upper-level Divergence (250 hPa)--------------------------
for idx_time in range(len(date)):
    # Select time slice and parse metadata
    ds = file.sel(time=date[idx_time]).metpy.parse_cf()

    # Variables
    u_250 = ds['u'][:].values
    v_250 = ds['v'][:].values
    z = ds['z']
    divergence = ds['d'][:].values * 1e5  # divergence scaled (x10^-5 s^-1)

    # Winds converted to units
    u_wind_knots = ds['u'].values * units.knots
    v_wind_knots = ds['v'].values * units.knots

    u_wind_m = units('m/s') * ds.variables['u'][:].squeeze()
    v_wind_m = units('m/s') * ds.variables['v'][:].squeeze()

    lat = ds.variables['latitude'][:].squeeze()
    lon = ds.variables['longitude'][:].squeeze()

    # Calculate wind magnitude
    mag_v = np.sqrt(u_250**2 + v_250**2)

    # Create 2D meshgrid of longitude and latitude
    lon_2d, lat_2d = np.meshgrid(lon, lat)

    # Smooth geopotential height field
    Z_250 = ndimage.gaussian_filter(z, sigma=3, order=0) * units.meter

    #---------------------Figure---------------------------------------
    fig, ax = plt.subplots(figsize=(12, 9))

    # Basemap setup for rectangular projection
    bmap = Basemap(llcrnrlon=extent[0], llcrnrlat=extent[1],
                   urcrnrlon=extent[2], urcrnrlat=extent[3], 
                   epsg=4326, resolution='i')


    bmap.readshapefile('/content/BR_UF_2019', 'BRA_UF_2019', linewidth=0.3, color='black')
    bmap.readshapefile('/content/South_America', 'South_America', linewidth=0.3, color='black')

    
    bmap.drawparallels(np.arange(-90, 90, 5), linewidth=0.2, dashes=[4,4], color='white', 
                      labels=[True, False, False, True], fmt='%g', labelstyle="+/-", size=12)
    bmap.drawmeridians(np.arange(0, 360, 5), linewidth=0.2, dashes=[4,4], color='white', 
                      labels=[True, False, False, True], fmt='%g', labelstyle="+/-", size=12)

    # Define grid for plotting
    x = np.linspace(min_lon, max_lon, divergence.shape[1])
    y = np.linspace(max_lat, min_lat, divergence.shape[0])
    xx, yy = np.meshgrid(x, y)

    # Contour levels for divergence
    divergence_levels = np.arange(0, 20, 2)

    # Filled contour 
    cf = plt.contourf(xx, yy, divergence, cmap='Blues', levels=divergence_levels, extend='max')

    cbar_div = plt.colorbar(cf, ticks=divergence_levels, format='%.1f', 
                           orientation='horizontal', fraction=0.040, pad=0.05)
    cbar_div.set_label('Divergence ($10^{-5}$ s$^{-1}$)', fontsize=14)
    
  
  # Streamplot for wind vectors colored by magnitude
    mag_levels = np.arange(10, 80, 10)
  
    strm = plt.streamplot(lon_2d, lat_2d, u_250, v_250, density=3, linewidth=1,
                         color=mag_v, cmap='magma_r', arrowsize=0.9,
                         arrowstyle='-|>', minlength=0.15, maxlength=4.0)


    cbar_wind = plt.colorbar(strm.lines, format='%.0f', extend='both', pad=0.03)
    cbar_wind.set_label('Wind Speed (m s$^{-1}$)', fontsize=14)


    plt.title('Divergence ($10^{-5}$ s$^{-1}$) and Streamlines (m s$^{-1}$) at 250 hPa',
              loc='center', fontsize=12, y=1.05)
    plt.title('ERA', loc='left', fontsize=9)
    plt.title(f"{pd.to_datetime(str(ds.time.values)).strftime('%Y-%m-%d %H')} UTC",
              loc='right', fontsize=9)


    filename = pd.to_datetime(str(ds.time.values)).strftime('%Y-%m-%d_%H') + "_250_div.png"
    plt.savefig(filename, dpi=300)
    plt.close(fig)
