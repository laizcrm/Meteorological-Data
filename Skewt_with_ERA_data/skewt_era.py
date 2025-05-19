##################################################
#      SKEW-T PLOT WITH ERA MODEL DATA
##################################################
# Some regions in Brazil do not have direct observational data
# for vertical wind profiles. Therefore, this script was created
# as an alternative to obtain vertical wind profiles at a specific
# location using ERA5 model data.


# Install required packages
!pip install netCDF4  
!pip install metpy==1.0


import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import pandas as pd
import xarray as xr
from metpy.units import units
import metpy.calc as mpcalc
from metpy.plots import SkewT, Hodograph
import matplotlib.gridspec as gridspec
import os



file = xr.open_dataset("/content/perfilvert1422utc.nc")
print(file)

# Define slice variables for subset domain
lat_slice = -25
lon_slice = -55
level_slice = slice(1000, 10, 100)

# Select data by latitude, longitude, and time
data = file.sel(latitude=lat_slice,
                longitude=lon_slice,
                time='2016-07-14T22:00:00')
#print(data)

# Convert data to DataFrame
df = data.to_dataframe()

# Reset index for level and reverse order
df2 = df.reset_index(level='level')
df3 = df2.iloc[::-1]

# Variables
T = df3['t'][:]
ur = df3['r'][:]
p = df3['level'].values * units.hPa

# Calculate dew point temperature approximation (in Celsius)
df3['Td'] = T - ((100 - ur) / 5) - 273.15
Td = df3['Td']

# Temperature in Celsius
df3['tc'] = T - 273.15
tc = df3['tc']

# Wind components
u = df3['u']
v = df3['v']

# Export to CSV
df3.to_csv('sondagemfoz22as22utc.csv')

# Reload CSV (if necessary)
sond = pd.read_csv('/content/sondagemfoz22as22utc.csv', sep=';', header=0)

T = sond['t'][:]
ur = sond['r'][:]
p = sond['level'].values * units.hPa

# Recalculate dew point and temperature with units
sond['Td'] = T - ((100 - ur) / 5) - 273.15
Td = sond['Td'].values * units.degC
sond['tc'] = (T - 273.15).values * units.degC
tc = sond['tc'].values * units.degC
u = sond['u'].values * units('meters/second')
v = sond['v'].values * units('meters/second')

# Calculate LCL
lcl_pressure, lcl_temperature = mpcalc.lcl(p[0], tc[0], Td[0])

# Calculate parcel profile temperature
parcel_profile = mpcalc.parcel_profile(p, tc[0], Td[0]).to('degC')

print('LCL Pressure:', lcl_pressure)
print('LCL Temperature:', lcl_temperature)

# ---------- PLOT SETTINGS ----------
DPI = 300
output_file = 'skewt_plot.png'

fig = plt.figure(figsize=(12, 15))

gs = gridspec.GridSpec(2, 2, width_ratios=[1, 2], height_ratios=[4, 1])

skew = SkewT(fig, rotation=45, subplot=gs[1])
skew.plot(p, tc, 'g')
skew.plot(p, Td, 'b')

skew.plot_barbs(p[::2], u[::2], v[::2], length=7, fill_empty=False, sizes=dict(emptybarb=0., spacing=0.2, height=0.3))

skew.ax.set_ylim(1000, 100)

plt.title('Skew-T Plot ERA5', loc='center', fontsize=16, y=1.05)
plt.title('Foz do Iguaçu (PR), Brazil', loc='left', fontsize=9)
plt.title('14 July 2016 22:00 UTC', loc='right', fontsize=9)

plt.xlabel('Temperature (°C)', fontsize=12)
plt.ylabel('Pressure (hPa)', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

skew.plot_dry_adiabats()
skew.plot_moist_adiabats()
skew.plot_mixing_lines()

skew.ax.set_xlim(-30, 40)

skew.ax.axvline(0, color='c', linestyle='--', linewidth=2)
skew.ax.axvline(-20, color='c', linestyle='--', linewidth=2)

ax = plt.subplot(gs[0])
h = Hodograph(ax, component_range=60)
h.add_grid(increment=10)
h.plot(u, v, linewidth=1.5)
plt.title('Hodograph', loc='center', fontsize=13)

plt.savefig(output_file, dpi=DPI, bbox_inches='tight', pad_inches=0)
os.system('open -a Preview ' + output_file)
plt.show()
