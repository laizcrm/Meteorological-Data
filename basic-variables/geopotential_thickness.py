import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pygrib
from glob import glob
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader

# Get list of files matching pattern and sort them
files = glob("../gfsanl_4_*")
file_list = sorted(files)



# Define map extent (west_lon, east_lon, south_lat, north_lat)
extent = [-85.0, -60.0, -20.0, 10.0]
min_lon, max_lon = extent[0], extent[1]
min_lat, max_lat = extent[2], extent[3]


for file_path in file_list:
    # Open GRIB file using pygrib
    grib = pygrib.open(file_path)
    
    # Extract valid time from GRIB messages
    for g in grib:
        time = g.validDate
    
    # Select relevant variables from GRIB file
    pres = grib.select(name='Pressure reduced to MSL')[0]  # Mean sea level pressure (Pa)
    sup = grib.select(name='Surface pressure')[0]          # Surface pressure (Pa)
    ht_500 = grib.select(name='Geopotential Height', typeOfLevel='isobaricInhPa', level=500)[0]     # Geopotential height at 500 hPa (gpm)
    ht_1000 = grib.select(name='Geopotential Height', typeOfLevel='isobaricInhPa', level=1000)[0]   # Geopotential height at 1000 hPa (gpm)
    u_wind = grib.select(name='U component of wind', typeOfLevel='isobaricInhPa', level=250)[0]     # U wind component at 250 hPa (m/s)
    v_wind = grib.select(name='V component of wind', typeOfLevel='isobaricInhPa', level=250)[0]     # V wind component at 250 hPa (m/s)
    
    # Extract data over specified lat/lon extent (note: longitudes shifted by +360)
    msl, lats, lons = pres.data(min_lat, max_lat, min_lon + 360, max_lon + 360)
    u = u_wind.data(min_lat, max_lat, min_lon + 360, max_lon + 360)[0]
    v = v_wind.data(min_lat, max_lat, min_lon + 360, max_lon + 360)[0]
    surface_pressure = sup.data(min_lat, max_lat, min_lon + 360, max_lon + 360)[0]
    ht_500_data = ht_500.data(min_lat, max_lat, min_lon + 360, max_lon + 360)[0]
    ht_1000_data = ht_1000.data(min_lat, max_lat, min_lon + 360, max_lon + 360)[0]
    
    # Calculate wind speed magnitude at 250 hPa
    wspd = np.sqrt(u**2 + v**2)
  
    # Calculate geopotential thickness (in decameters) between 1000 hPa and 500 hPa
    thickness = (ht_500_data - ht_1000_data) / 10
  
    #Convert MSL pressure to hPa
    msl_hpa = msl / 100
    
    #-----------------Figure----------------------------------------------------
    fig = plt.figure(figsize=(10,10))
    ax = plt.axes(projection=ccrs.PlateCarree())
    
    # Add land feature 
    ax.add_feature(cfeature.LAND, facecolor='gray', alpha=0.2)
    ax.set_extent([extent[0], extent[1], extent[2], extent[3]], crs=ccrs.PlateCarree())
    
    shapefile = list(shpreader.Reader('shapefile/ne_10m_admin_1_states_provinces.shp').geometries())
    ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor='gray', facecolor='none', linewidth=0.3)
    
    # Add coastlines, borders and gridlines
    ax.coastlines(resolution='10m', color='gray', linewidth=0.8)
    ax.add_feature(cfeature.BORDERS, edgecolor='gray', linewidth=0.5)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), color='gray', alpha=1.0, linestyle='--', linewidth=0.25,
                      xlocs=np.arange(-180, 180, 5), ylocs=np.arange(-90, 90, 5), draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False

  
    # Define contour levels
    wind_min = np.min(wspd)
    wind_max = np.max(wspd)
    contour_interval = 3
    levels = np.arange(40, wind_max, contour_interval)
    
    # Plot geopotential thickness contours (dashed brown lines)
    thickness_contours = ax.contour(lons, lats, thickness, colors='brown', linestyles='dashed', alpha=0.6)
    plt.clabel(thickness_contours, inline=True, inline_spacing=3, fontsize=12, fmt='%1.0f', colors='darkred')
    
    # Plot mean sea level pressure contours (solid black lines)
    mslp_contours = ax.contour(lons, lats, msl_hpa, colors='black', linestyles='solid', alpha=0.6)
    plt.clabel(mslp_contours, inline=True, inline_spacing=8, fontsize=12, fmt='%1.0f', colors='black')
    
    # Plot wind speed as shaded contour plot
    wind_speed_plot = ax.contourf(lons, lats, wspd, cmap='BuPu', levels=levels, extend='max')
    

    cbar = plt.colorbar(wind_speed_plot, ticks=np.arange(40, wind_max, contour_interval), format='%.0f',
                        pad=0.01, fraction=0.04)
    cbar.set_label('(m $s^{-1}$)', fontsize=14)
    
    # titles
    plt.title('Wind Speed (m $s^{-1}$) at 250 hPa, Thickness (dam) 1000-500 hPa, and MSL Pressure (hPa)', 
              loc='center', fontsize=14, y=1.05)
    plt.title('GFS 0.25°', loc='left', fontsize=13)
    plt.title(f'{pd.to_datetime(time).strftime("%Y-%m-%d %H")} UTC', fontsize=13, loc='right')
    

    plt.savefig(f'./esp/{pd.to_datetime(time).strftime("%Y-%m-%d_%H")}_thickness.png')
    
    print('---------------------------------------')
    print('------ Script executed for file ------')
    print('File time:', time)
  
