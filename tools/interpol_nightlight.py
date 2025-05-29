"""
GOES-16 Satellite Visualization Tool with Night Lights Overlay

This script processes GOES-16 satellite data (Channels 07 and 13) to create an RGB grayscale composite.
It remaps the infrared data to geographic coordinates (lat/lon), overlays night lights from the 
VIIRS Black Marble dataset, and interpolates the lights to match the satellite grid.

Outputs:
- 'nightlights.png' : visualization of night lights and cold cloud tops
- 'night_light.npy' : interpolated night lights array
"""

import matplotlib.pyplot as plt
plt.switch_backend('agg')

from netCDF4 import Dataset
import numpy as np
from pyproj import Proj
from remap import remap
from osgeo import gdal
from scipy.interpolate import griddata
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')
import time as t

start = t.time()

# -----------------------------
# Input GOES-16 NetCDF files
# A GOES-16 file with half day and half night (visível + IR)
#the CONUS domain
# -----------------------------
file_ch13 = 'dados/13/OR_ABI-L2-CMIPF-M6C13_G16_s20221522100206_e20221522109525_c20221522110015.nc'
file_ch07 = 'dados/07/OR_ABI-L2-CMIPF-M6C07_G16_s20221522050206_e20221522059525_c20221522059594.nc'
extent = [-130.0, -55.0, -25.0, 70.0]  # [lon_min, lon_max, lat_min, lat_max]
resolution = 1.0

# -----------------------------
# Read Channel 13 (IR) and projection info
# -----------------------------
nc = Dataset(file_ch13)
proj = nc.variables['goes_imager_projection']
H = proj.perspective_point_height
sat_lon = proj.longitude_of_projection_origin
sat_sweep = proj.sweep_angle_axis

x1 = nc.variables['x_image_bounds'][0, 0] * H
x2 = nc.variables['x_image_bounds'][0, 1] * H
y1 = nc.variables['y_image_bounds'][0, 1] * H
y2 = nc.variables['y_image_bounds'][0, 0] * H

print('Remapping Channel 13...')
grid = remap(file_ch13, extent, resolution, x1, y1, x2, y2)
data_ch13 = grid.ReadAsArray()

# The projection x and y coordinates equals- the scanning angle (in radians) multiplied by the satellite height
X = nc.variables['x'][:] * H
Y = nc.variables['y'][:] * H

# Convert map points to latitude and longitude with the magic provided by Pyproj :p
XX, YY = np.meshgrid(X, Y)
p = Proj(proj='geos', h=H, lon_0=sat_lon, sweep=sat_sweep)
lons, lats = p(XX, YY, inverse=True)

print('Projection and coordinates set.')
#print('lon_0= ')
#print(lats)
# -----------------------------
# Read and process Channel 07
# -----------------------------
print('Remapping Channel 07...')
nc07 = Dataset(file_ch07)
grid07 = remap(file_ch07, extent, resolution, x1, y1, x2, y2)
data_ch07 = grid07.ReadAsArray()

#Convert bright to reflectance 
# Normalize the channel between a range.
#cleanIR = (cleanIR-minimumValue)/(maximumValue-minimumValue)
cleanIR07 = (data_ch07 - 90) / (313 - 90)
cleanIR07 = np.clip(cleanIR07, 0, 1) # Apply range limits to make sure values are between 0 and 1

cleanIR07 = 1 - cleanIR07
cleanIR07 = cleanIR07 / 1.4
RGB_cleanIR07 = np.dstack([cleanIR07]*3)

print('Channel 07 processed for grayscale RGB.')

# -----------------------------
# Load and crop night lights
# -----------------------------
print('Opening night lights raster...')
#crs = raster.crs
#EPSG: 4326

#bounds = raster.bounds
#boundingBox(left=-180, bottom=-90, right=180,top=90)

raster = gdal.Open('BlackMarble_2016_3km_gray_geo.tif')
ulx, xres, _, uly, _, yres = raster.GetGeoTransform()
lrx = ulx + (raster.RasterXSize * xres)
# lrx = 180
lry = uly + (raster.RasterYSize * yres)
# lry = 90


min_lon, max_lat, max_lon, min_lat = extent[0], extent[3], extent[1], extent[2]
raster = gdal.Translate('teste.tif', raster, projWin=[min_lon, max_lat, max_lon, min_lat])
array = raster.ReadAsArray()

r, g, b = array[0].astype(float), array[1].astype(float), array[2].astype(float)

geo = raster.GetGeoTransform()
xmin, xres = geo[0], geo[1]
ymax, yres = geo[3], geo[5]
xmax = xmin + (xres * raster.RasterXSize)
ymin = ymax + (yres * raster.RasterYSize)

lons_n = np.arange(xmin, xmax, xres)
lats_n = np.arange(ymax, ymin, yres)
lons_n, lats_n = np.meshgrid(lons_n, lats_n)
#print('LATS e LONS do NL', lats_n.shape,lons_n.shape)  #shape(2250,2063)
#print('R:',r.shape)    #type: float64
#print('					')

length = len(lats_n) * len(lons_n)
print('Night light grid shape:', lats_n.shape, lons_n.shape)

points = np.empty((length, 2))
points[:, 1] = lats_n.flatten()
points[:, 0] = lons_n.flatten()

#print('points',points.shape)

print('Interpolating night lights to satellite grid...')
interpol = griddata(points, r.flatten(), (lons, lats), method='cubic')

np.save('night_light.npy', interpol) #how to use > np.load('luzes.npy')

plt.imsave('nightlights.png', interpol, cmap='gray')

print('--- Finished! ---')
print('--- Execution Time:', t.time() - start, 'seconds ---')
