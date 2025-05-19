"""
This script reads WRF model output data and generates a vertical cross-section 
plot of radar reflectivity (dBZ) along a specified latitude/longitude line.
It interpolates reflectivity values vertically between two geographic points, 
fills gaps near the terrain, and overlays the terrain profile on the cross-section. 
The plot uses color shading based on the NWS reflectivity color scale and labels 
the x-axis with latitude and longitude coordinates along the cross-section.
"""



import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from wrf import getvar, vertcross, CoordPair, to_np, interpline, get_cartopy
from matplotlib.colors import from_levels_and_colors

# Open WRF output file
wrf_file = Dataset('run2/wrfout_d01_2016-07-14_22:00:00')

# Define cross-section start and end points (lat, lon)
start_point = CoordPair(lat=-26, lon=-57)
end_point = CoordPair(lat=-26, lon=-52)

# Extract variables at last time step
height = getvar(wrf_file, "z", timeidx=-1)
terrain = getvar(wrf_file, "ter", timeidx=-1)
dbz = getvar(wrf_file, "dbz", timeidx=-1)

# Convert dBZ to linear reflectivity for interpolation
Z = 10 ** (dbz / 10.)

# Calculate vertical cross-section
z_cross = vertcross(Z, height, wrfin=wrf_file,
                    start_point=start_point,
                    end_point=end_point,
                    latlon=True, meta=True)

# Convert back to dBZ
dbz_cross = 10 * np.log10(z_cross)
dbz_cross.attrs.update(z_cross.attrs)
dbz_cross.attrs["description"] = "Radar reflectivity cross section"
dbz_cross.attrs["units"] = "dBZ"

# Fill lower missing values with first non-missing value in each column
dbz_filled = np.ma.copy(to_np(dbz_cross))
for i in range(dbz_filled.shape[-1]):
    column = dbz_filled[:, i]
    first_valid = int(np.transpose((column > -200).nonzero())[0])
    dbz_filled[0:first_valid, i] = dbz_filled[first_valid, i]

# Terrain profile along cross-section
terrain_line = interpline(terrain, wrfin=wrf_file,
                         start_point=start_point,
                         end_point=end_point)

# Setup plot
fig, ax = plt.subplots(figsize=(10, 6))

# Define reflectivity levels and colors (NWS standard)
levels = np.arange(5, 75, 5)
colors_rgb = np.array([[4,233,231],[1,159,244],[3,0,244],[2,253,2],[1,197,1],[0,142,0],
                       [253,248,2],[229,188,0],[253,149,0],[253,0,0],[212,0,0],[188,0,0],
                       [248,0,253],[152,84,198]]) / 255
cmap, norm = from_levels_and_colors(levels, colors_rgb, extend='max')

# Cross-section axes (distance vs height)
x = np.arange(dbz_cross.shape[-1])
y = to_np(dbz_cross.coords["vertical"])

# Plot filled contours of reflectivity
cf = ax.contourf(x, y, to_np(dbz_filled), levels=levels,
                 cmap=cmap, norm=norm, extend='max')

# Add colorbar
cbar = fig.colorbar(cf, ax=ax)
cbar.set_label('Reflectivity (dBZ)', fontsize=10)
cbar.ax.tick_params(labelsize=8)

# Plot terrain
ax.fill_between(x, 0, to_np(terrain_line), color='saddlebrown')

# Format x-axis with lat/lon labels
coords = to_np(dbz_cross.coords["xy_loc"])
num_ticks = 5
step = max(1, int(len(coords) / num_ticks))
ax.set_xticks(x[::step])
ax.set_xticklabels([coord.latlon_str() for coord in coords[::step]], rotation=45, fontsize=8)

# Labels and title
ax.set_xlabel('Longitude', fontsize=10)
ax.set_ylabel('Height (m)', fontsize=10)
ax.set_title('Vertical Cross-Section of Radar Reflectivity (dBZ)', fontsize=14)

# Save figure
plt.tight_layout()
plt.savefig('reflectivity_cross_section.png')
plt.show()
