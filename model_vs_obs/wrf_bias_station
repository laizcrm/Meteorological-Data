import datetime as dt
import wrf
from netCDF4 import num2date, Dataset
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pygrib

#---WRF Libraries----
import wrf as wrf
from wrf import getvar, ALL_TIMES, latlon_coords,get_cartopy, to_np, smooth2d,cartopy_xlim, cartopy_ylim,interplevel, extract_vars


#-----------WRF Data------------------------------
file_wrf = Dataset('wrfout_d01_2022-07-01_00:00:00')


#Extract variables
wrfin = extract_vars(file_wrf, ALL_TIMES,("T2", 'PSFC' ))

#Extract time variable from WRF file across all time steps using concatenation method
#here are 49 time steps, hourly intervals from 2022-07-01 00 UTC to 2022-07-03 00 UTC

time=getvar(file_wrf,"times",timeidx=ALL_TIMES, method='cat')
#print(time)

# 2-meter temperature
t2 = wrfin['T2']

# Surface pressure
PSFC = wrfin['PSFC']

#-----------------------------------------
# Locate the nearest grid point to a given latitude and longitude (Curitiba)
# and extract the temperature time series at that location
#-----------------------------------------

lats, lons = latlon_coords(PSFC)

#select the coordinates for Curitiba
# Curitiba coordinates: [-25.44, -49.23] approximately correspond to grid indices [269, 255]

lat_est = -25.44
lon_est = -49.23

# Find the nearest grid points 
abs_lat = np.abs(np.array(lats) - lat_est)
abs_lon = np.abs(np.array(lons) - lon_est)

c = np.maximum(abs_lat, abs_lon)

# Find the index of the minimum difference (closest grid point)
latlon_idx = np.argmin(c)

# Convert the flat index back to 2D grid indices (x, y)
x, y = np.where(c == np.min(c))

# Extract the temperature time series at the nearest grid point
grid_temp = t2[:, x[0], y[0]]
# Example indices: x=128, y=228 (grid point closest to Curitiba)

# Convert temperature from Kelvin to Celsius
temp_local = grid_temp - 273.15

#------ Plot WRF temperature data over time ------

fig = plt.figure(figsize=(10,10))


plt.plot(time,temp_local,label='WRF',linewidth=1.5)

plt.title('2-meter Temperature',loc='center', fontsize=14,y=1.05)


plt.xlabel('Time [hours]',fontsize=14)
plt.ylabel('Temperature [°C]',fontsize=14)
plt.legend(loc='upper left')


plt.title('WRF- Exp1', loc='left', fontsize=12) 
#plt.title('{}'.format(str(time.values) +' UTC'), fontsize=9,loc='right')
  
plt.savefig("t2m_exp1.png",dpi=150)
plt.legend(loc='upper left')


#-------------- GFS Data --------------

path = ('../../gfs/forecasting/01/00/recorte_gfs.t00z.pgrb2.0p25.f000')


grib= pygrib.open(path)

temp = grib.select(name='2 metre temperature')[0]

extent = [-60, -35, -40, -18]  # Simepar domain region

min_lon = extent[0]; max_lon = extent[2]; min_lat = extent[1]; max_lat = extent[3]

# Extract temperature data and corresponding latitudes and longitudes within the domain
temp, lats_gfs, lons_gfs = temp.data(min_lat, max_lat, min_lon+360, max_lon+360)

lons_gfs = lons_gfs - 360

# Find the nearest grid points 
abs_lat_gfs = np.abs(np.array(lats_gfs) - lat_est)
abs_lon_gfs = np.abs(np.array(lons_gfs) - lon_est)


c = np.maximum(abs_lat_gfs, abs_lon_gfs)

latlon_idx = np.argmin(c)

# Convert the flat index back to 2D grid indices (xg, yg)
xg, yg = np.where(c == np.min(c))

# Index for the location (you can uncomment these lines to print lat/lon of the grid point)
# print(lats_gfs[xg[0], yg[0]])
# print(lons_gfs[xg[0], yg[0]])

# Select GFS files by date using a loop and extract 2-meter temperature for Curitiba

time_gfs = []  # List to store datetime objects for each forecast hour
t2_gfs =   []  # List to store 2-meter temperature values at Curitiba grid point

# Initial datetime (start of forecast)
di = dt.datetime(2022, 7, 1, 0)




# Loop through 49 forecast hours (from 0 to 48)

for i in range(49):
    # Calculate the current forecast datetime by adding i hours to the initial datetime
    d = di + dt.timedelta(hours=i)
    time_gfs.append(d)  # Append current datetime to time list
    
    path = ('../../gfs/forecasting/01/00/recorte_gfs.t00z.pgrb2.0p25.f0%02d' % i)

    # Open the GRIB file
    grib = pygrib.open(path)
  
    temp = grib.select(name='2 metre temperature')[0]
    
    temp, lats_gfs, lons_gfs = temp.data(min_lat, max_lat, min_lon + 360, max_lon + 360)



# Append the temperature at the Curitiba grid point to the list
t2_gfs.append(temp[xg[0], yg[0]])

#convert to array
t2_gfs = np.array(t2_gfs)-273.15



#--------- Plot time series of 2m temperature from WRF and GFS ---------

g = plt.figure(figsize=(10,10))

# Plot WRF temperature time series
plt.plot(time, temp_local, label='WRF', linewidth=1.5)

# Plot GFS temperature time series
plt.plot(time_gfs, t2_gfs, label='GFS', linewidth=1.5)

# Set axis labels
plt.xlabel('Time [hours]', fontsize=14)
plt.ylabel('Temperature [°C]', fontsize=14)

# Show legend at upper left
plt.legend(loc='upper left')

# Main title centered above plot
plt.title('2m Temperature from WRF and GFS', loc='center', fontsize=14, y=1.05)

# Subtitle with experiment names on the left
plt.title('WRF-Exp1 and GFS', loc='left', fontsize=12)

# Save the figure with specified dpi
plt.savefig("t2m_WRFEGFS_exp1.png", dpi=150)


#---------------- Plot Bias (WRF - GFS) of 2m Temperature -----------------

fig = plt.figure(figsize=(10,10))

#BIAS time series
plt.plot(time, temp_local - t2_gfs, label='Bias (WRF - GFS)', linewidth=1.5)

# Line for reference
plt.plot(time, np.full(len(t2_gfs), 0), linestyle='--', color='gray')


plt.xlabel('Time [hours]', fontsize=14)
plt.ylabel('Temperature [°C]', fontsize=14)

plt.legend(loc='upper left')


plt.title('Bias of 2m Temperature between WRF and GFS', loc='center', fontsize=14, y=1.05)
plt.title('WRF-Exp1 and GFS', loc='left', fontsize=12)
plt.savefig("t2m_bias_exp1.png", dpi=150)
