"""
This script processes GOES-16 satellite data to generate a false color composite image
combining visible and infrared channels over the CONUS domain.

It reads multiple netCDF files corresponding to different spectral bands,
extracts reflectance and brightness temperature data, and applies
contrast enhancement and a false color look-up table (LUT) for visualization.

The script also computes geolocation coordinates (latitude and longitude)
from satellite projection data to enable mapping and spatial analysis.
"""
#import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt
plt.switch_backend('agg')

from netCDF4 import Dataset
 # Import the Basemap toolkit
import numpy as np # Import the Numpy package

from matplotlib.colors import LinearSegmentedColormap # Linear interpolation for color maps

from cpt_convert import loadCPT # Import the CPT convert function
from pyproj import Proj
from remap import remap # Import the Remap function
from pyorbital import astronomy
from osgeo import gdal

#---Cartopy 
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader

import glob, os
import subprocess

import locale


# A GOES-16 file with half day and half night (visível mais IR)

from datetime import datetime, timedelta
from skimage import exposure #histogram equalization

import metpy  # noqa: F401
import numpy as np
import xarray as xr

#Import do script ára fazer a conversão para lat/lon
# Miscellaneous operating system interfaces
from utilities import download_CMI       # Our own utilities
from utilities import geo2grid, latlon2xy, convertExtent2GOESProjection  
import warnings

warnings.filterwarnings('ignore')

#---------------- PARTE I----------------------------- 
#Abrir o arquivo for the CONUS domain
print("Reading Files...")

f13 = ('dados/13/OR_ABI-L2-CMIPF-M6C13_G16_s20221522100206_e20221522109525_c20221522110015.nc')

f7 = ('dados/07/OR_ABI-L2-CMIPF-M6C07_G16_s20221522050206_e20221522059525_c20221522059594.nc')


#-------RGB-------------------
#Conversão para gridArray
#-----------------------------
f1 = ('dados/01/OR_ABI-L2-CMIPF-M6C01_G16_s20221522050206_e20221522059514_c20221522059597.nc')


f2 = ('dados/02/OR_ABI-L2-CMIPF-M6C02_G16_s20221522040206_e20221522049514_c20221522049596.nc')


f3 = ('dados/03/OR_ABI-L2-CMIPF-M6C03_G16_s20221522050206_e20221522059514_c20221522059591.nc')


canal1 = Dataset(f1)
canal2 = Dataset(f2)
canal3 = Dataset(f3)
canal7 = Dataset(f7)
canal13 = Dataset(f13)

#  RGB: 

B=canal1.variables['CMI'][0][::2,::2]
R=canal2.variables['CMI'][0][::2,::2]
G=canal3.variables['CMI'][0][::2,::2]

#IR
data7=canal7.variables['CMI'][0][::1,::1]
data13=canal13.variables['CMI'][0][::1,::1]



print("Loading RGB ...")



# The RGB array for the true color image
RGB = np.dstack([R, G, B])
#print(RGB.shape) #(linha,coluna,canal)=(5424,5424,3)

# Load the Channel 02 contrast curve
contrast_curve = [0.00000, 0.02576, 0.05148, 0.07712, 0.10264, 0.12799, 0.15313, 0.17803, 0.20264, 0.22692, 0.25083, 0.27432, 0.29737, 0.31991, 0.34193, 0.36336,
0.38418, 0.40433, 0.42379, 0.44250, 0.46043, 0.47754, 0.49378, 0.50911, 0.52350, 0.53690, 0.54926, 0.56055, 0.57073, 0.57976, 0.58984, 0.59659,
0.60321, 0.60969, 0.61604, 0.62226, 0.62835, 0.63432, 0.64016, 0.64588, 0.65147, 0.65694, 0.66230, 0.66754, 0.67267, 0.67768, 0.68258, 0.68738,
0.69206, 0.69664, 0.70112, 0.70549, 0.70976, 0.71394, 0.71802, 0.72200, 0.72589, 0.72968, 0.73339, 0.73701, 0.74055, 0.74399, 0.74736, 0.75065,
0.75385, 0.75698, 0.76003, 0.76301, 0.76592, 0.76875, 0.77152, 0.77422, 0.77686, 0.77943, 0.78194, 0.78439, 0.78679, 0.78912, 0.79140, 0.79363,
0.79581, 0.79794, 0.80002, 0.80206, 0.80405, 0.80600, 0.80791, 0.80978, 0.81162, 0.81342, 0.81518, 0.81692, 0.81862, 0.82030, 0.82195, 0.82358,
0.82518, 0.82676, 0.82833, 0.82987, 0.83140, 0.83292, 0.83442, 0.83592, 0.83740, 0.83888, 0.84036, 0.84183, 0.84329, 0.84476, 0.84623, 0.84771,
0.84919, 0.85068, 0.85217, 0.85368, 0.85520, 0.85674, 0.85829, 0.85986, 0.86145, 0.86306, 0.86469, 0.86635, 0.86803, 0.86974, 0.87149, 0.87326,
0.87500, 0.87681, 0.87861, 0.88038, 0.88214, 0.88388, 0.88560, 0.88730, 0.88898, 0.89064, 0.89228, 0.89391, 0.89552, 0.89711, 0.89868, 0.90023,
0.90177, 0.90329, 0.90479, 0.90627, 0.90774, 0.90919, 0.91063, 0.91205, 0.91345, 0.91483, 0.91620, 0.91756, 0.91890, 0.92022, 0.92153, 0.92282,
0.92410, 0.92536, 0.92661, 0.92784, 0.92906, 0.93027, 0.93146, 0.93263, 0.93380, 0.93495, 0.93608, 0.93720, 0.93831, 0.93941, 0.94050, 0.94157,
0.94263, 0.94367, 0.94471, 0.94573, 0.94674, 0.94774, 0.94872, 0.94970, 0.95066, 0.95162, 0.95256, 0.95349, 0.95441, 0.95532, 0.95622, 0.95711,
0.95799, 0.95886, 0.95973, 0.96058, 0.96142, 0.96225, 0.96307, 0.96389, 0.96469, 0.96549, 0.96628, 0.96706, 0.96783, 0.96860, 0.96936, 0.97010,
0.97085, 0.97158, 0.97231, 0.97303, 0.97374, 0.97445, 0.97515, 0.97584, 0.97653, 0.97721, 0.97789, 0.97856, 0.97922, 0.97988, 0.98053, 0.98118,
0.98182, 0.98246, 0.98309, 0.98372, 0.98435, 0.98497, 0.98559, 0.98620, 0.98681, 0.98741, 0.98802, 0.98862, 0.98921, 0.98980, 0.99039, 0.99098,
0.99157, 0.99215, 0.99273, 0.99331, 0.99389, 0.99446, 0.99503, 0.99561, 0.99618, 0.99675, 0.99732, 0.99788, 0.99845, 0.99902, 0.99959, 1.000000]

# Convert the contrast curve to a numpy array
curve = np.array(contrast_curve)
#print(np.shape(contrast_curve))

#Junta os três canais em media 
RGB_mean=np.mean(RGB,axis=2) #tons de cinza

rgb_convert=RGB_mean
 
  
# Minimum and maximum reflectances
VISmin = 0.0
VISmax = 1.0

# Anything that is below the min or greater than max, keep min and max
rgb_convert[rgb_convert>VISmax] = VISmax
# Convert to 0 - 255
rgb_convert = ((rgb_convert - VISmin) / (VISmax - VISmin)) * 255


# Convert to int
rgb_convert = rgb_convert.astype(int)
# Apply the contrast curve

# Apply the contrast curve
rgb_convert = curve[rgb_convert]* 255 

rgb_convert = rgb_convert.astype(int)

print('rgb_convert:  ',rgb_convert.shape)


print("Loading Infrared ...")

#------------Clean IR--------

#Converter brilho para refletância 
cleanIR13 = data13 #serve para compor o RGB



IR13=data13 #plot separado de IR



# Minimum and maximum brightness temperatures
IRmin = 89.62
IRmax = 341.27
# Anything that is below the min or greater than max, keep min and max
cleanIR13[cleanIR13>IRmax] = IRmax
# Convert to 0 - 255
cleanIR13 = ((cleanIR13 - IRmin) / (IRmax - IRmin)) * 255
# Convert to int
cleanIR13 = cleanIR13.astype(int)


print('cleanIR13:  ', cleanIR13.shape)


#----Canal 07-----------------



#-----------------PARTE II-------------------------------



print("Reading the False Color Look Up Table...")
 
import matplotlib.image as mpimg
# Open the False Color Look Up Table
img = mpimg.imread('wx-star.com_GOES-R_ABI_False-Color-LUT_sat20.png')
# Flip the array (for some reason, is necessary to flip the LUT horizontally)
img = np.fliplr(img)

print('---IMG shape:--- ',img.shape)

def is_int(n):
  try :
     if np.isnan(n) == False:
      n=int(n)
      return True
     else:
      return False
  except:
     return False


final_image=[]

for jx in range(rgb_convert.shape[0]):
  for jy in range(rgb_convert.shape[1]):
    #print(data1[jx,jy])
    #print(data2[jx,jy])
    if is_int(rgb_convert[jx,jy])==True and is_int(cleanIR13[jx,jy])==True:
      r = img[int(rgb_convert[jx,jy]),int(cleanIR13[jx,jy]),0:3]
    else:
      r=[0,0,0] #cor preta 
    final_image.append(r)


final_image=np.reshape(final_image,(rgb_convert.shape[0],rgb_convert.shape[1],3))

print('Final Image:  ',final_image.shape)

print("Calculating the lons and lats...")
# Satellite height
sat_h = canal7.variables['goes_imager_projection'].perspective_point_height
# Satellite longitude
sat_lon = canal7.variables['goes_imager_projection'].longitude_of_projection_origin
# Satellite sweep
sat_sweep = canal7.variables['goes_imager_projection'].sweep_angle_axis
# The projection x and y coordinates equals
# the scanning angle (in radians) multiplied by the satellite height (http://proj4.org/projections/geos.html)
X = canal7.variables['x'][:][::1] * sat_h
Y = canal7.variables['y'][:][::1] * sat_h
# map object with pyproj
p = Proj(proj='geos', h=sat_h, lon_0=sat_lon, sweep=sat_sweep, a=6378137.0)
# Convert map points to latitude and longitude with the magic provided by Pyproj
XX, YY = np.meshgrid(X, Y)
lons, lats = p(XX, YY, inverse=True)
 
print('LATS e LONS   ', lons.shape) 
 
print("                                 ")
print("Calculating the sun zenith angle...")


utc_time = datetime(2022, 6, 1, 21, 00)
sun_zenith = np.zeros((rgb_convert.shape[0], rgb_convert.shape[1]))
sun_zenith = astronomy.sun_zenith_angle(utc_time, lons, lats)
 
print("sun_zenith",sun_zenith.shape) 
 
print("Putting zeroes in the areas without sunlight")


final_image[sun_zenith > 80] = 0
mask = (final_image == [0.0,0.0,0.0]).all(axis=2)
#apply the mask to overwrite the pixels
final_image[mask] = [0,0,0]
alpha = ~np.all(final_image == 0, axis=2)
final_image_alpha = np.dstack((final_image, alpha))
print(np.shape(final_image))
print(np.shape(final_image_alpha))   

raster = gdal.Open('BlackMarble_2016_3km_gray_geo.tif') 
ulx, xres, xskew, uly, yskew, yres = raster.GetGeoTransform()
lrx = ulx + (raster.RasterXSize * xres)
lry = uly + (raster.RasterYSize * yres)
corners = [ulx, lry, lrx, uly]
#extent = [-135.0, -75.0, -14.0, 72.0]
extent = [-85.0,-50.0,-30.0,10.0]
min_lon = extent[0]; max_lon = extent[2]; min_lat = extent[1]; max_lat = extent[3]
raster = gdal.Translate('teste.tif', raster, projWin = [min_lon, max_lat, max_lon, min_lat])
array2 = raster.ReadAsArray()
r = array2[0,:,:].astype(float)
g = array2[1,:,:].astype(float)
b = array2[2,:,:].astype(float)
#r[r==4] = 1
#g[g==5] = 4
#b[b==15] = 30
r[r==4] = 0
g[g==5] = 0
b[b==15] = 0
geo = raster.GetGeoTransform()
xres = geo[1]
yres = geo[5]
xmin = geo[0]
xmax = geo[0] + (xres * raster.RasterXSize)
ymin = geo[3] + (yres * raster.RasterYSize)
ymax = geo[3]
lons_n = np.arange(xmin,xmax,xres)
lats_n = np.arange(ymax,ymin,yres)
lons_n, lats_n = np.meshgrid(lons_n,lats_n)
color_tuples = (np.array([r[:-1,:-1].flatten(), g[:-1,:-1].flatten(), b[:-1,:-1].flatten()]).transpose())/255

 
 # Define the size of the saved picture=========================================
DPI = 150
fig = plt.figure(figsize=(rgb_convert.shape[1]/float(DPI), rgb_convert.shape[0]/float(DPI)), frameon=False, dpi=DPI)
#fig = plt.figure(figsize=(10,6), frameon=False, dpi=DPI)
ax = plt.Axes(fig, [0., 0., 1., 1.])
ax.set_axis_off()
fig.add_axes(ax)
ax = plt.axis('off')
#==============================================================================

# Converts a CPT file to be used in Python
cpt = loadCPT('IR4AVHRR6.cpt')
IR4AVHRR6 = LinearSegmentedColormap('cpt', cpt)



# Use the Geostationary projection in cartopy
ax = plt.axes(projection=ccrs.Geostationary(central_longitude=-75.0, satellite_height=35786023.0))
img_extent = (-5434894.67527,5434894.67527,-5434894.67527,5434894.67527)
 
 
 


print("First layer...")# TEMPERATURA DE BRILHO



IR07 = data7 #plot separado de IR

IR07[sun_zenith < 80] = np.nan
#-----------Para plotar 07+13--------------------

# Apply range limits for clean IR channel 13:

IR13 = np.maximum(IR13, 90)
IR13 = np.minimum(IR13, 313)
# Normalize the channel between a range
IR13 = (IR13-90)/(313-90)
# Invert colors
IR13 = 1 - IR13

# Apply range limits for clean IR channel 07:

IR07 = np.maximum(IR07, 90)
IR07 = np.minimum(IR07, 313)
# Normalize the channel between a range
IR07 = (IR07-90)/(313-90)
# Invert colors
IR07 = 1 - IR07


#Variável para plotar apenas IR 
IR = np.maximum(IR13,IR07)


ax.imshow(IR,cmap='gray', origin='upper', zorder=2)


print("Second layer...")  # TEMPERATURA DE BRILHO
# SECOND LAYER
plt.pcolormesh(lons_n, lats_n, r, color=color_tuples, zorder=1)


print("Third layer...")
# FIFTH LAYER
ax.imshow(final_image_alpha, origin='upper', zorder=2)




print("Saving Image...")

plt.savefig('geocolor.png', dpi=DPI, facecolor="black", pad_inches=0)


print ('- finished! ') 
 






