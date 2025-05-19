"""
Horizontal cross-section plot at 6000 m altitude from WRF output.

This script generates a horizontal cross-section at 6000 meters showing:
- Vertical wind speed (W) as shaded contours
- Ice water mixing ratio (snow + graupel) as solid contour lines
"""


mport datetime
import numpy as np
from numpy import linspace, meshgrid
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cartopy.crs as crs
import cartopy.feature as cfeature
from netCDF4 import Dataset
from wrf import getvar, extract_vars, latlon_coords, destagger

# Load WRF output file
ds = Dataset('../rod_2/wrfout_d02_2016-07-14_22:00:00')

# Extract variables
ht = getvar(ds, "z", timeidx=-1)
vars = extract_vars(ds, timeidx=-1, varnames=("QGRAUP", "QSNOW", "QRAIN", 'U', 'W', 'V'))

qrain = vars["QRAIN"]
qsnow = vars["QSNOW"]
qgraup = vars["QGRAUP"]
w = destagger(vars["W"], 0, meta=True)

# Calculate total ice water mixing ratio (g/kg)
ice_ratio = (qsnow + qgraup) * 1e3

# Select horizontal level near 6000 m (e.g., bottom_top=37)
w_cross = w.sel(bottom_top=37)
ice_cross = ice_ratio.sel(bottom_top=37)

# Get latitude and longitude
lats, lons = latlon_coords(ice_ratio)

# Define plot extent (lon_min, lat_min, lon_max, lat_max)
extent = [-53.75, -27.0, -52.55, -26.40]
step = 0.1
min_lon, max_lon, min_lat, max_lat = extent

# Create figure and axis with map projection
fig = plt.figure(figsize=(10, 10))
ax = plt.axes(projection=crs.PlateCarree())
ax.set_extent(extent, crs=crs.PlateCarree())

# Add background features
ax.add_feature(cfeature.LAND, facecolor='gray', alpha=0.4)
gl = ax.gridlines(draw_labels=True, linewidth=0.25, linestyle='--', color='gray')
gl.top_labels = False
gl.right_labels = False
gl.xlocator = plt.MultipleLocator(step)
gl.ylocator = plt.MultipleLocator(step)

# Vertical wind (W) shaded
interval = 0.2
levels = np.arange(-2, 2, interval)
colors_list = list(plt.cm.PRGn(np.linspace(0, 1, len(levels))))
colors_list[9] = "white"
colors_list[10] = "white"
cmap = colors.ListedColormap(colors_list)
norm = colors.BoundaryNorm(boundaries=levels, ncolors=len(levels))

w_plot = ax.contourf(to_np(lons), to_np(lats), to_np(w_cross),
                     levels=levels, cmap=cmap, norm=norm, extend="both")
plt.colorbar(w_plot, ticks=np.arange(-2, 2, interval), format='%.1f', pad=0.01, fraction=0.025)

# Ice water mixing ratio contours
ice_levels = [0.1, 0.5, 1.0, 1.5, 2.0]
ice_contours = ax.contour(to_np(lons), to_np(lats), to_np(ice_cross),
                          levels=ice_levels, colors='black', linewidths=1.6)
plt.clabel(ice_contours, inline=1, fontsize=12, fmt='%1.1f')

# Titles
ax.set_title("Horizontal cross-section at 6000 m: W (shaded) and ice water mixing ratio (g/kg)", fontsize=13, loc='center', y=1.05)
plt.title("WRF - 3 km resolution", fontsize=10, loc='left')

# Save figure
plt.savefig('wh.png', dpi=300)
