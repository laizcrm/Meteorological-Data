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

#-------------------- STANDARD LIBRARIES --------------------
import os
import glob
import numpy as np
import pandas as pd
import datetime as dt
from datetime import datetime

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
from osgeo import gdal


# Load the NetCDF file containing 850 hPa atmospheric data 
file = xr.open_dataset("/content/nivel_850.nc")
print(file)


# Time range
# timestamps every 6 hours between July 13 and July 15, 2016
date = list(pd.date_range('2016-07-13T00:00:00.00','2016-07-15T00:00:00.00', freq='6H').strftime('%Y-%m-%d %H'))


extent = [-95.0,-60.0,-20.0,10.0]
min_lon = extent[0]; max_lon = extent[2]; min_lat = extent[1]; max_lat = extent[3]


#----------------------- TEMPERATURE ADVECTION PLOTTING -----------------------

for idx_tempo in range(len(date)):

    # Select dataset for the current time
    ds = file.sel(time=date[idx_tempo]).metpy.parse_cf()

    # Variables
    lat = ds['latitude'][:].squeeze()
    lon = ds['longitude'][:].squeeze()
    u = ds['u']
    v = ds['v']
    t = ds['t']
    z = ds['z'][:]

    # Convert wind to knots and m/s
    u_wind = u.values * units.knots
    v_wind = v.values * units.knots
    u_wind_m = u.values * units('m/s')
    v_wind_m = v.values * units('m/s')

    # Create 2D lat/lon grid
    lon_2d, lat_2d = np.meshgrid(lon, lat)

    # Flip flag for wind barbs in Southern Hemisphere
    flip_flag = np.zeros((u.shape[0], v.shape[1]))  
    flip_flag[lat_2d < 0] = 1 

    # Calculate grid deltas
    dx, dy = mpcalc.lat_lon_grid_deltas(lon, lat)

    # Calculate temperature advection (K -> °C conversion inside function)
    adv = mpcalc.advection(t - 273.15, u, v) * 1e4

    # Smooth fields
    Z_850 = ndimage.gaussian_filter(z, sigma=3) * units.meter
    adv = ndimage.gaussian_filter(adv, sigma=3)

    #-------------figure-------------------------------------
    fig = plt.figure(figsize=(10, 10))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([extent[0], extent[2], extent[1], extent[3]], crs=ccrs.PlateCarree())

    shapefile = list(shpreader.Reader('ne_10m_admin_1_states_provinces.shp').geometries())
    ax.add_geometries(shapefile, crs=ccrs.PlateCarree(), edgecolor='gray', facecolor='none', linewidth=0.3)

    ax.coastlines(resolution='10m', color='black', linewidth=0.8)
    ax.add_feature(cartopy.feature.BORDERS, edgecolor='black', linewidth=0.5)


    gl = ax.gridlines(crs=ccrs.PlateCarree(), color='gray', linestyle='--', linewidth=0.25,
                      xlocs=np.arange(-180, 180, 5), ylocs=np.arange(-90, 90, 5), draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False


    x = linspace(min_lon, max_lon, adv.shape[1])
    y = linspace(max_lat, min_lat, adv.shape[0])
    xx, yy = meshgrid(x, y)

    # Advection color levels
    levels = np.arange(-3, 3, 0.5)
    cmap = plt.cm.bwr(np.linspace(0, 1, len(levels)))
    cmap[5:7] = "white"  # highlight 0 value
    cmap = colors.ListedColormap(cmap)
    norm = colors.BoundaryNorm(levels, len(levels))

    # Contours and filled contours
    contour = ax.contour(xx, yy, Z_850 / 10, colors='black', linestyles='dashed', linewidths=1.0)
    plt.clabel(contour, inline=True, fontsize=12, fmt='%1.0f', colors='black')

    shaded = ax.contourf(xx, yy, adv, levels=levels, cmap=cmap, norm=norm, extend='both')
    cbar = plt.colorbar(shaded, ticks=levels, format='%.1f', pad=0.009, shrink=0.8)
    cbar.set_label('Advection (°C)', fontsize=10)

    # Wind vectors (quiver)
    ax.quiver(lon_2d[::8, ::8], lat_2d[::8, ::8], u.values[::8, ::8], v.values[::8, ::8], color='gray')


    time_str = pd.to_datetime(str(ds.time.values)).strftime('%Y-%m-%d %H UTC')
    plt.title('Temperature Advection (°C) and Geopotential Height (dam)', loc='center', fontsize=13, y=1.05)
    plt.title('ERA5 0.25°', loc='left', fontsize=9)
    plt.title(time_str, loc='right', fontsize=9)


    filename = f"{time_str.replace(':', '').replace(' ', '_')}_adv_t.png"
    plt.savefig(filename, dpi=150)
    plt.close()



# --------- SPECIFIC HUMIDITY ADVECTION  ---------


for time_idx in range(len(date)):
    # Select dataset for the given time
    ds = file.sel(time=date[time_idx]).metpy.parse_cf()

    # Variables
    lat = ds.variables['latitude'][:].squeeze()
    lon = ds.variables['longitude'][:].squeeze()
    u = ds['u']
    v = ds['v']
    q = ds['q'][:] * 1000  # Convert specific humidity from kg/kg to g/kg
    z = ds['z'][:]
    t = ds['t']

    # Convert winds to knots (if needed for other purposes)
    u_wind_knots = ds['u'].values * units.knots
    v_wind_knots = ds['v'].values * units.knots

    # Convert winds to m/s
    u_wind_ms = units('m/s') * ds.variables['u'][:].squeeze()
    v_wind_ms = units('m/s') * ds.variables['v'][:].squeeze()

    # Create 2D lat/lon grid
    lon_2d, lat_2d = np.meshgrid(lon, lat)

    # Flag for flipping barbs below equator
    flip_flag = np.zeros((u_wind_knots.shape[0], v_wind_knots.shape[1]))
    flip_flag[lat_2d < 0] = 1

    # Calculate grid spacing in meters between lat/lon points
    dx, dy = mpcalc.lat_lon_grid_deltas(lon, lat)

    # Calculate specific humidity advection (units scaled)
    adv_q = mpcalc.advection(q, u, v) * 1e4  # Scaling factor

    # Smooth geopotential height and advection fields
    Z_850_smoothed = ndimage.gaussian_filter(z, sigma=3, order=0) * units.meter
    adv_q_smoothed = ndimage.gaussian_filter(adv_q, sigma=3, order=0)

    #-------------figure-------------------------------------
    fig = plt.figure(figsize=(10, 10))
    ax = plt.axes(projection=ccrs.PlateCarree())

    img_extent = [extent[0], extent[2], extent[1], extent[3]]
    ax.set_extent([extent[0], extent[2], extent[1], extent[3]], crs=ccrs.PlateCarree())


    shapefile = list(shpreader.Reader('ne_10m_admin_1_states_provinces.shp').geometries())
    ax.add_geometries(shapefile, crs=ccrs.PlateCarree(), edgecolor='gray', facecolor='none', linewidth=0.3)

    # Add coastlines, borders, and gridlines
    ax.coastlines(resolution='10m', color='black', linewidth=0.8)
    ax.add_feature(cartopy.feature.BORDERS, edgecolor='black', linewidth=0.5)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), color='gray', alpha=1.0, linestyle='--', linewidth=0.25,
                      xlocs=np.arange(-180, 180, 5), ylocs=np.arange(-90, 90, 5), draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False

    # Generate regularly spaced grid for contour plotting
    x = np.linspace(min_lon, max_lon, adv_q_smoothed.shape[1])
    y = np.linspace(max_lat, min_lat, adv_q_smoothed.shape[0])
    xx, yy = np.meshgrid(x, y)

    # Define contour levels and colors
    min_val = np.min(adv_q_smoothed)
    max_val = np.max(adv_q_smoothed)
    interval = 0.5
    levels = np.arange(-2.5, 2.5, interval)

    # Create a custom colormap with white near zero values for emphasis
    ncolors = list(plt.cm.RdYlGn(np.linspace(0, 1, len(levels))))
    ncolors[4] = "white"
    ncolors[5] = "white"
    cmap = colors.ListedColormap(ncolors, "", len(ncolors))
    norm = colors.BoundaryNorm(levels, len(levels), cmap.N)

    # Contour geopotential height (divided by 10 to get decameters)
    geopot_contours = ax.contour(xx, yy, Z_850_smoothed / 10, colors='black',
                                linewidths=1.0, linestyles='dashed', alpha=0.8)
    plt.clabel(geopot_contours, inline=1, inline_spacing=8, fontsize=12, fmt='%1.0f', colors='black')

    # Filled contour
    adv_contourf = ax.contourf(xx, yy, adv_q_smoothed, cmap=cmap, levels=levels,
                               norm=norm, extend='both')
    cbar = plt.colorbar(adv_contourf, ticks=np.arange(-2.5, 2.5, interval), format='%.1f',
                       pad=0.009, shrink=0.8)
    cbar.set_label('Specific Humidity Advection (g/kg)')

    # Wind vectors (quiver) 
    plt.quiver(lon_2d[::8, ::8], lat_2d[::8, ::8], u[::8, ::8], v[::8, ::8], color='gray')


    plt.title('Specific Humidity Advection (g/kg) and Geopotential Height (dam)',
              fontsize=12, y=1.05)
    plt.title('ERA5 0.25°', loc='left', fontsize=9)
    plt.title(f"{pd.to_datetime(str(ds.time.values), format='%Y-%m-%d %H')} UTC", fontsize=9, loc='right')

    filename = f"{pd.to_datetime(str(ds.time.values), format='%Y-%m-%d %H')}_adv_q.png"
    plt.savefig(filename, dpi=DPI)
    plt.close()



# --------- WIND MAGNITUDE AND GEOPOTENTIAL HEIGHT PLOTTING ---------

for time_idx in range(len(date)):
    # Select dataset for current time step
    ds = file.sel(time=date[time_idx]).metpy.parse_cf()

    # Variables
    lat = ds.variables['latitude'][:].squeeze()
    lon = ds.variables['longitude'][:].squeeze()
    u = ds['u']
    v = ds['v']
    # q = ds['q'][:] * 1000  # Specific humidity (not used here)
    z = ds['z'][:]
    t = ds['t']

    # Calculate wind magnitude (magnitude of the horizontal wind vector)
    wind_speed = np.sqrt(u**2 + v**2)

    # Convert winds to knots (for barbs plotting)
    u_wind_knots = ds['u'].values * units.knots
    v_wind_knots = ds['v'].values * units.knots

    # Convert winds to m/s (if needed elsewhere)
    u_wind_ms = units('m/s') * ds.variables['u'][:].squeeze()
    v_wind_ms = units('m/s') * ds.variables['v'][:].squeeze()

    # Create 2D latitude-longitude grid
    lon_2d, lat_2d = np.meshgrid(lon, lat)

    # Create flag array for flipping barbs below the equator
    flip_flag = np.zeros(u_wind_knots.shape)
    flip_flag[lat_2d < 0] = 1

    # Calculate grid spacing in meters (not used explicitly here, but useful)
    dx, dy = mpcalc.lat_lon_grid_deltas(lon, lat)

    # Smooth geopotential height field
    Z_850_smoothed = ndimage.gaussian_filter(z, sigma=3, order=0) * units.meter

 #-------------figure-------------------------------------
    fig = plt.figure(figsize=(10, 10))
    ax = plt.axes(projection=ccrs.PlateCarree())

    ax.set_extent([extent[0], extent[2], extent[1], extent[3]], crs=ccrs.PlateCarree())

    shapefile = list(shpreader.Reader('ne_10m_admin_1_states_provinces.shp').geometries())
    ax.add_geometries(shapefile, crs=ccrs.PlateCarree(), edgecolor='gray', facecolor='none', linewidth=0.3)

    ax.coastlines(resolution='10m', color='black', linewidth=0.8)
    ax.add_feature(cartopy.feature.BORDERS, edgecolor='black', linewidth=0.5)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), color='gray', alpha=1.0, linestyle='--', linewidth=0.25,
                      xlocs=np.arange(-180, 180, 5), ylocs=np.arange(-90, 90, 5), draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False

    x = np.linspace(min_lon, max_lon, wind_speed.shape[1])
    y = np.linspace(max_lat, min_lat, wind_speed.shape[0])
    xx, yy = np.meshgrid(x, y)

    min_speed = np.min(wind_speed)
    max_speed = np.max(wind_speed)
    interval = 2
    levels = np.arange(12, max_speed, interval)

    # Filled contour 
    speed_contourf = ax.contourf(xx, yy, wind_speed, cmap='viridis_r', levels=levels, extend='max')
    cbar = plt.colorbar(speed_contourf, ticks=np.arange(12, max_speed, interval), format='%.1f', pad=0.009, shrink=0.8)
    cbar.set_label('Wind Speed (m/s)')

    # Contour lines of geopotential height (in decameters)
    geopot_contours = ax.contour(xx, yy, Z_850_smoothed / 10, colors='darkred',
                                linewidths=1.2, linestyles='solid', alpha=0.8)
    plt.clabel(geopot_contours, inline=1, fontsize=12, fmt='%1.0f', colors='darkred')

    # Plot wind barbs 
    plt.barbs(lon_2d[::14, ::14], lat_2d[::14, ::14],
              u_wind_knots[::14, ::14], v_wind_knots[::14, ::14],
              flip_barb=flip_flag[::14, ::14],
              length=5.5, fill_empty=False,
              sizes=dict(emptybarb=0., spacing=0.2, height=0.3),
              alpha=0.8)

    # Optional: plot wind vectors or streamlines (commented out)
    # plt.quiver(lon_2d[::8, ::8], lat_2d[::8, ::8], u[::8, ::8], v[::8, ::8], color='gray')
    # plt.streamplot(lon_2d, lat_2d, u, v, density=2, linewidth=1.5, color='b', arrowsize=1)


    plt.title('Wind Speed (m/s) and Geopotential Height (dam)', fontsize=12, y=1.05)
    plt.title('ERA5 0.25°', loc='left', fontsize=9)
    plt.title(f"{pd.to_datetime(str(ds.time.values), format='%Y-%m-%d %H')} UTC", fontsize=9, loc='right')

  
    filename = f"{pd.to_datetime(str(ds.time.values), format='%Y-%m-%d %H')}_wind.png"
    plt.savefig(filename, dpi=DPI)
    plt.close()


