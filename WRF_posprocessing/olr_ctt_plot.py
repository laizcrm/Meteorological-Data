"""
Script to extract and plot Outgoing Longwave Radiation (OLR) and Cloud Top Temperature (CTT)
from WRF model output files over a specified region.

- Uses WRF-Python package to handle WRF output files.
- Applies smoothing and Gaussian filter to CTT.
- Plots OLR as shaded contours and CTT as contour lines.
- Saves plots for each time step.
  """


import os
from glob import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter
import cartopy.crs as crs
import cartopy.io.shapereader as shpreader
from netCDF4 import Dataset
from wrf import getvar, ALL_TIMES, latlon_coords, smooth2d, to_np, extract_vars

# Load WRF output files
files = sorted(glob("../rod_2/wrfout_d02_*"))
wrfin = [Dataset(f) for f in files]

# Extract OLR variable for all times
wrf_data = extract_vars(wrfin, ALL_TIMES, ("OLR",))

# Extract times from WRF files
times = getvar(wrfin, "times", timeidx=ALL_TIMES, method='cat')

# Define region extent [min_lon, min_lat, max_lon, max_lat]
extent = [-60, -46, -31, -22]
interval = 5  # Gridline interval for map

# Loop over each time step and plot
for idx, ds in enumerate(wrfin):
    olr = wrf_data["OLR"][idx]
    time = times[idx]
    ctt_raw = getvar(ds, "ctt", timeidx=idx)

    # Get lat/lon coordinates
    lats, lons = latlon_coords(olr)

    # Smooth cloud top temperature
    ctt_smooth = smooth2d(ctt_raw, 10, cenweight=4)
    ctt = gaussian_filter(ctt_smooth, sigma=3)

    # Set up plot
    fig = plt.figure(figsize=(10, 10))
    ax = plt.axes(projection=crs.PlateCarree())
    ax.set_extent(extent, crs.PlateCarree())

    # Add shapefile boundaries
    shapefile = list(shpreader.Reader('shapefile/ne_10m_admin_1_states_provinces.shp').geometries())
    ax.add_geometries(shapefile, crs.PlateCarree(), edgecolor='gray', facecolor='none', linewidth=0.3)

    # Add coastlines, borders and gridlines
    ax.coastlines(resolution='10m', color='black', linewidth=0.8)
    ax.add_feature(cartopy.feature.BORDERS, edgecolor='black', linewidth=0.5)

    gl = ax.gridlines(
        crs=crs.PlateCarree(), color='gray', alpha=1.0, linestyle='--', linewidth=0.25,
        xlocs=np.arange(extent[0], extent[2], interval),
        ylocs=np.arange(extent[1], extent[3], interval),
        draw_labels=True
    )

    # Define contour levels
    olr_levels = np.arange(100, 320, 20)
    ctt_levels = np.arange(-80, 0, 10)

    # Plot OLR as filled contours
    olr_cf = ax.contourf(to_np(lons), to_np(lats), olr, levels=olr_levels, cmap='Greys', extend='max')

    # Colorbar for OLR
    cbar = plt.colorbar(olr_cf, ticks=olr_levels, format="%.0f", pad=0.01, fraction=0.03)
    cbar.set_label('W/m²', fontsize=14)

    # Plot CTT as contour lines
    ctt_contours = ax.contour(
        to_np(lons), to_np(lats), ctt, levels=4, colors='yellow', linestyles='solid', alpha=0.6, transform=crs.PlateCarree()
    )
    plt.clabel(ctt_contours, inline=True, fontsize=13, fmt='%1.0f', colors='red')

    # Titles
    plt.title('Outgoing Longwave Radiation (OLR)', loc='center', fontsize=16, y=1.05)
    plt.title('WRF - 9 Km Resolution', loc='left', fontsize=12)
    plt.title(f"{pd.to_datetime(str(time.values)).strftime('%Y-%m-%d %H UTC')}", fontsize=9, loc='right')

    # Save figure
    output_dir = './olr'
    os.makedirs(output_dir, exist_ok=True)
    filename = f"{output_dir}/{pd.to_datetime(str(time.values)).strftime('%Y-%m-%d_%H%M')}_olr9km.png"
    plt.savefig(filename, dpi=300)
    plt.close(fig)

    print(f"Processed and saved plot for time: {time.values}")
