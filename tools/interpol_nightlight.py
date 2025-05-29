"""
GOES-16 Satellite Visualization Tool with Night Lights Overlay

This script processes GOES-16 satellite data (Channels 07 and 13) to create an RGB grayscale composite.
It remaps the infrared data to geographic coordinates (lat/lon), overlays night lights from the 
VIIRS Black Marble dataset, and interpolates the lights to match the satellite grid.

Output:
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

# ----------------------------------------
# Input GOES-16 NetCDF files (Ch.13 = IR, Ch.07 = clean IR)
# ----------------------------------------

file_ch13 = 'dados/13/OR_ABI-L2-CMIPF-M6C13_G16_s20221522100206_e20221522109525_c20221522110015.nc'
file_ch07 = 'dados/07/OR_ABI-L2-CMIPF-M6C07_G16_s20221522050206_e20221522059525_c20221522059594.nc'

# Domain extent [lon_min, lon_max, lat_min, lat_max]
extent = [-130.0, -55.0, -25.0, 70.0]
resolution = 1.0

# ----------------------------------------
# Read Channel 13 (IR) and extract projection info
# ----------------------------------------

nc = Dataset(file_ch13)
proj = nc.variables['goes_imager_projection']
H = proj.perspective_point_height
sat_lon = proj.longitude_of_projection_origin
sat_sweep = proj.sweep_angle_axis

# Define remapping boundaries (in meters)
x1 = nc.variables['x_image_bounds'][0, 0] * H
x2 = nc.variables['x_image_bounds'][0, 1] * H
y1 = nc.variables['y_image_bounds'][0, 1] * H
y2 = nc.variables['y_image_bounds'][0, 0] * H

# Remap channel 13
grid = remap(file_ch13, extent, resolution, x1, y1, x2, y2)
data_ch13 = grid.ReadAsArray()

# Get lat/lon from projection
X = nc.variables['x'][:] * H
Y = nc.variables['y'][:] * H
XX, YY = np.meshgrid(X, Y)
p = Proj(proj='geos', h=H, lon_0=sat_lon, sweep=sat_sweep)
lons, lats = p(XX, YY, inverse=True)

# ----------------------------------------
# Process Channel 07 for grayscale RGB image
# ----------------------------------------

nc07 = Dataset(file_ch07)
grid07 = remap(file_ch07, extent, resolution, x1, y1, x2, y2)
data_ch07 = grid07.ReadAsArray()

# Normalize and convert to grayscale (invert cold cloud brightness)
cleanIR07 = (data_ch07 - 90) / (313 - 90)
cleanIR07 = np.clip(cleanIR07, 0, 1)
cleanIR07 = 1 - cleanIR07
cleanIR07 = cleanIR07 / 1.4
RGB_cleanIR07 = np.dstack([cleanIR07]*3)

# ----------------------------------------
# Load and crop night lights (Black Marble)
# ----------------------------------------

raster = gdal.Open('BlackMarble_2016_3km_gray_geo.tif')
ulx, xres, _, uly, _, yres = raster.GetGeoTransform()
lrx = ulx + (raster.RasterXSize * xres)
lry = uly + (raster.RasterYSize * yres)
min_lon, max_lat, max_lon, min_lat = extent[0], extent[3], extent[1], extent[2]

# Crop raster to domain
raster = gdal.Translate('teste.tif', raster, projWin=[min_lon, max_lat, max_lon, min_lat])
array = raster.ReadAsArray()
r, g, b = array[0].astype(float), array[1].astype(float), array[2].astype(float)

# Get night lights coordinates
geo = raster.GetGeoTransform()
xmin, xres = geo[0], geo[1]
ymax, yres = geo[3], geo[5]
xmax = xmin + (xres * raster.RasterXSize)
ymin = ymax + (yres * raster.RasterYSize)

lons_n = np.arange(xmin, xmax, xres)
lats_n = np.arange(ymax, ymin, yres)
lons_n, lats_n = np.meshgrid(lons_n, lats_n)

# Interpolate night lights to satellite grid
points = np.column_stack((lons_n.flatten(), lats_n.flatten()))
interpol = griddata(points, r.flatten(), (lons, lats), method='cubic')

# Save outputs
np.save('night_light.npy', interpol)

plt.imsave('nightlights.png', interpol, cmap='gray')

print('--- Finished! ---')
