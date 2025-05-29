#import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt
plt.switch_backend('agg')

from netCDF4 import Dataset
 # Import the Basemap toolkit
import numpy as np # Import the Numpy package

from matplotlib.colors import LinearSegmentedColormap # Linear interpolation for color maps
#from mpl_toolkits.basemap import Basemap
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
from scipy.interpolate import griddata
import locale
import time as t
#time 
start = t.time()


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



#Abrir o arquivo for the CONUS domain

file_in = ('dados/13/OR_ABI-L2-CMIPF-M6C13_G16_s20221522100206_e20221522109525_c20221522110015.nc')

file_in2 = ('dados/07/OR_ABI-L2-CMIPF-M6C07_G16_s20221522050206_e20221522059525_c20221522059594.nc')

#Canal 13


file_out ='PLOT.png'
# OPEN THE FILE USING THE NETCDF4 LIBRARY
nc = Dataset(file_in)


perspective_point_height = nc.variables['goes_imager_projection'].perspective_point_height

# Extent of data in decimais (2712*0.000056*35786023.0)
#xmin = nc.variables['x'][:].min()
#xmax = nc.variables['x'][:].max()
#ymin = nc.variables['y'][:].min()
#ymax = nc.variables['y'][:].max()

#*perspective_point_height
#extent = [xmin, xmax, ymin, ymax]

extent = [-130.0,-55.0,-25.0,70.0]


# Scan's start time, converted to datetime object
scan_start = datetime.strptime(nc.time_coverage_start, '%Y-%m-%dT%H:%M:%S.%fZ')


resolution = 1.0 
H = nc.variables['goes_imager_projection'].perspective_point_height
x1 = nc.variables['x_image_bounds'][0,0] * H
x2 = nc.variables['x_image_bounds'][0,1] * H
y1 = nc.variables['y_image_bounds'][0,1] * H
y2 = nc.variables['y_image_bounds'][0,0] * H
#print('H= '+ str(H)) 
grid= remap(file_in, extent, resolution, x1,y1,x2,y2)
data = grid.ReadAsArray() 

ori_proj = nc.variables['goes_imager_projection']
# Satellite height
sat_h = ori_proj.perspective_point_height

# Satellite longitude
sat_lon = ori_proj.longitude_of_projection_origin

# Satellite sweep
sat_sweep = ori_proj.sweep_angle_axis

# The projection x and y coordinates equals
# the scanning angle (in radians) multiplied by the satellite height (http://proj4.org/projections/geos.html)
X = nc.variables['x'][:] * sat_h
Y = nc.variables['y'][:] * sat_h

p = Proj(proj='geos', h=sat_h, lon_0=sat_lon, sweep=sat_sweep)


# Convert map points to latitude and longitude with the magic provided by Pyproj
XX, YY = np.meshgrid(X, Y)
lons, lats = p(XX, YY, inverse=True)
 
#print('lon_0= ')
#print(lats)


#----Canal 07-----------------

nc = Dataset(file_in2)

grid07= remap(file_in2,extent,resolution, x1,y1,x2,y2)
data07 = grid07.ReadAsArray() 

#------------RGB--------
#Converter brilho para refletância 
cleanIR07 = data07 

# Normalize the channel between a range.
#       cleanIR = (cleanIR-minimumValue)/(maximumValue-minimumValue)
cleanIR07 = (cleanIR07-90)/(313-90)

# Apply range limits to make sure values are between 0 and 1
cleanIR07 = np.clip(cleanIR07, 0, 1)

# Invert colors so that cold clouds are white
cleanIR07= 1 - cleanIR07

# Lessen the brightness of the coldest clouds so they don't appear so bright
# when we overlay it on the true color image.
cleanIR07 = cleanIR07/1.4

# Yes, we still need 3 channels as RGB values. This will be a grey image.
RGB_cleanIR07 = np.dstack([cleanIR07, cleanIR07, cleanIR07])







# ----- INSERIR NIGHT LIGHTS -------------



#crs = raster.crs
#EPSG: 4326


#bounds = raster.bounds
#boundingBox(left=-180, bottom=-90, right=180,top=90)




raster = gdal.Open('BlackMarble_2016_3km_gray_geo.tif') 
ulx, xres, xskew, uly, yskew, yres = raster.GetGeoTransform()



lrx = ulx + (raster.RasterXSize * xres)
# lrx = 180
lry = uly + (raster.RasterYSize * yres)
# lry = 90
corners = [ulx, lry, lrx, uly]




min_lon = extent[0]; max_lon = extent[2]; min_lat = extent[1]; max_lat = extent[3]


raster = gdal.Translate('teste.tif', raster, projWin = [ min_lon, max_lat, max_lon, min_lat])




array2 = raster.ReadAsArray()
r = array2[0,:,:].astype(float)
g = array2[1,:,:].astype(float)
b = array2[2,:,:].astype(float)

#erro aqui
color_tuples = (np.array([r[:-1,:-1].flatten(), g[:-1,:-1].flatten(), b[:-1,:-1].flatten()]).transpose())/255


geo = raster.GetGeoTransform()
xres = geo[1]
yres = geo[5]
xmin = geo[0]
xmax = geo[0] + (xres * raster.RasterXSize)
ymin = geo[3] + (yres * raster.RasterYSize)
ymax = geo[3]


#coords das luzes 
lons_n = np.arange(xmin,xmax,xres)
lats_n = np.arange(ymax,ymin,yres)



length=len(lats_n)*len(lons_n)
#5062500

lons_n, lats_n = np.meshgrid(lons_n,lats_n)
#print('LATS e LONS do NL', lats_n.shape,lons_n.shape)  #shape(2250,2063)
#print('R:',r.shape)    #type: float64
#print('					')

points = np.empty((length, 2))

print('lons',lons_n.flatten().shape)

#point_end=np.zeros((length,2))
#print(lats_n.flatten().shape)
#points[:, 0] = np.ravel(lats_n)
#points[:, 1] = np.ravel(lons_n)
points[:, 1] = lats_n.flatten() #modifica lat para 1D
points[:, 0] = lons_n.flatten()

#print('points',points.shape)



#print(r.shape)
interpol = griddata(points,r.flatten(),(lons,lats),method='cubic')


np.save('night_light.npy',interpol)
#como puxar> np.load('luzes.npy')







plt.savefig('nightlights.png')





print ('- finished! ') 
 
print ('------finished--------- Time:', t.time() - start, 'seconds')








 
