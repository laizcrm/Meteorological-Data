"""
MERGE Precipitation Product

The MERGE product combines ground-based precipitation observations with satellite estimates
from the GPM-IMERG dataset (version 06). It uses high-resolution (0.1º spatial, 30-min temporal)
satellite data corrected with around 2,500 additional ground stations to reduce bias.

The product provides precipitation estimates in two stages: an initial quick estimate (Early)
and a more accurate delayed estimate (Late). MERGE uses the Late dataset for historical archives.

Data are available hourly since 2010 and daily since 2000, covering latitudes 60ºN to 60ºS.

Variables in the data files include:
- prec: 24-hour accumulated precipitation
- nest: number of observations per grid point

Source and more info: http://ftp.cptec.inpe.br/modelos/tempo/MERGE/rozante_et.al.2010.pdf
"""
from datetime import datetime, timedelta
from google.cloud import storage
import pandas as pd
import geopandas as gpd
import glob
from shapely.geometry import Point
import os

# Initialize Google Cloud Storage client with your service account key JSON file path
SERVICE_ACCOUNT_KEY_PATH = '/path/to/your/service_account_key.json'
storage_client = storage.Client.from_service_account_json(SERVICE_ACCOUNT_KEY_PATH)

# Reference date for processing (change as needed)
reference_date = datetime(2024, 11, 12)
reference_date_str = reference_date.strftime('%Y_%m_%d')


def merge_daily_pentads(storage_client):
    bucket_name = 'your-gcs-bucket-name'  # Replace with your bucket name
    bucket = storage_client.get_bucket(bucket_name)
    local_dir = './merge_daily'  # Local directory to save files

    # Clean local directory before downloading new files
    for file in os.listdir(local_dir):
        file_path = os.path.join(local_dir, file)
        if os.path.isfile(file_path):
            os.remove(file_path)
            print(f"Deleted {file_path}")

    # List of pentad dates (last 5 days before reference_date)
    pentads = [(reference_date - timedelta(days=i)).strftime('%Y_%m_%d') for i in range(1, 6)]

    # Download pentad daily files
    for date_str in pentads:
        gcs_path = f'MERGE/DAILY/{date_str[:4]}/{date_str[5:7]}/rainfall_{date_str}.parquet'
        local_path = os.path.join(local_dir, f'rainfall_{date_str}.parquet')

        blob = bucket.blob(gcs_path)
        blob.download_to_filename(local_path)
        print(f"Downloaded {gcs_path}")

    # Merge all daily pentad files into one DataFrame
    file_list = [os.path.join(local_dir, f'rainfall_{date}.parquet') for date in pentads]

    merged_df = None

    for file in file_list:
        try:
            df = pd.read_parquet(file, engine='pyarrow')
            df['longitude'] = (df['longitude'] + 180) % 360 - 180  # Normalize longitude

            if {'latitude', 'longitude', 'prec'}.issubset(df.columns):
                if merged_df is None:
                    merged_df = df[['latitude', 'longitude', 'prec']].copy()
                else:
                    merged_df = merged_df.merge(df[['latitude', 'longitude', 'prec']],
                                                on=['latitude', 'longitude'],
                                                how='outer',
                                                suffixes=('', '_new'))
                    merged_df['prec'] = merged_df['prec'].fillna(0) + merged_df['prec_new'].fillna(0)
                    merged_df.drop(columns='prec_new', inplace=True)
            else:
                print(f"Missing columns in {file}")
        except FileNotFoundError:
            print(f"File not found: {file}")

    return merged_df


def merge_hourly_pentads(storage_client):
    bucket_name = 'your-gcs-bucket-name'  # Replace with your bucket name
    bucket = storage_client.get_bucket(bucket_name)
    local_dir = './merge_hourly'  # Local directory to save files

    # Clean local directory before downloading new files
    for file in os.listdir(local_dir):
        file_path = os.path.join(local_dir, file)
        if os.path.isfile(file_path):
            os.remove(file_path)
            print(f"Deleted {file_path}")

    # List of pentad dates (last 5 days before reference_date)
    pentads = [(reference_date - timedelta(days=i)).strftime('%Y_%m_%d') for i in range(1, 6)]

    # Download hourly pentad files (for each hour in 24h)
    for date_str in pentads:
        for hour in range(24):
            gcs_path = f'MERGE/HOURLY/{date_str[:4]}/{date_str[5:7]}/{date_str[8:]}/rainfall_{date_str}_{hour:02d}.parquet'
            local_path = os.path.join(local_dir, f'rainfall_{date_str}_{hour:02d}.parquet')

            blob = bucket.blob(gcs_path)
            blob.download_to_filename(local_path)
            print(f"Downloaded {gcs_path}")

    # Merge all hourly pentad files into one DataFrame
    file_list = glob.glob(os.path.join(local_dir, '*.parquet'))

    merged_df = None

    for file in file_list:
        try:
            df = pd.read_parquet(file, engine='pyarrow')
            df['longitude'] = (df['longitude'] + 180) % 360 - 180  # Normalize longitude

            if {'latitude', 'longitude', 'prec'}.issubset(df.columns):
                if merged_df is None:
                    merged_df = df[['latitude', 'longitude', 'prec']].copy()
                else:
                    merged_df = merged_df.merge(df[['latitude', 'longitude', 'prec']],
                                                on=['latitude', 'longitude'],
                                                how='outer',
                                                suffixes=('', '_new'))
                    merged_df['prec'] = merged_df['prec'].fillna(0) + merged_df['prec_new'].fillna(0)
                    merged_df.drop(columns='prec_new', inplace=True)
            else:
                print(f"Missing columns in {file}")
        except FileNotFoundError:
            print(f"File not found: {file}")

    return merged_df


# Execute data merging
hourly_data = merge_hourly_pentads(storage_client)
daily_data = merge_daily_pentads(storage_client)

# Load shapefile (adjust path accordingly, no company info)
shapefile_path = './shapefile/BR_UF_2019.shp'
state_shapes = gpd.read_file(shapefile_path)

# Create GeoDataFrame for hourly data points
hourly_geometry = [Point(xy) for xy in zip(hourly_data['longitude'], hourly_data['latitude'])]
gdf_hourly = gpd.GeoDataFrame(hourly_data, geometry=hourly_geometry)
gdf_hourly.set_crs("EPSG:4326", inplace=True)
gdf_hourly.to_crs(state_shapes.crs, inplace=True)

# Create GeoDataFrame for daily data points
daily_geometry = [Point(xy) for xy in zip(daily_data['longitude'], daily_data['latitude'])]
gdf_daily = gpd.GeoDataFrame(daily_data, geometry=daily_geometry)
gdf_daily.set_crs("EPSG:4326", inplace=True)
gdf_daily.to_crs(state_shapes.crs, inplace=True)

# Spatial join to merge data points with state polygons
gdf_hourly_merged = gpd.sjoin(gdf_hourly, state_shapes, how="inner")
gdf_daily_merged = gpd.sjoin(gdf_daily, state_shapes, how="inner")

 Calculate diff with merge to ensure correct alignment
diff = prec_day.merge(prec_hour, on='SIGLA_UF', suffixes=('_day', '_hour'))
diff['prec'] = diff['prec_day'] - diff['prec_hour']

colors = ["#b4f0f0", "#96d2fa", "#78b9fa", "#3c95f5", "#1e6deb", "#1463d2", "#0fa00f", "#28be28",
            "#50f050", "#72f06e", "#b3faaa", "#fff9aa", "#ffe978", "#ffc13c", "#ffa200", "#ff6200",
            "#ff3300", "#ff1500", "#c00100", "#a50200", "#870000", "#653b32"]
cmap = matplotlib.colors.ListedColormap(colors)

cmap.set_over('#000000')
cmap.set_under('#aaaaaa')
data_min = 20
data_max = 1000
interval = 50
levels = np.arange(data_min, data_max, interval)

norm = mcolors.BoundaryNorm(boundaries=levels, ncolors=cmap.N)

pentadas = [(date_test - timedelta(days=i)).strftime('%Y_%m_%d') for i in range(1, 6)]

# Plot DAILY
extent = [-80.0, -34.0, -45.0, 10.0]
fig, ax = plt.subplots(figsize=(20, 20), subplot_kw=dict(projection=ccrs.PlateCarree()))
ax.set_extent(extent, crs=ccrs.PlateCarree())
ax.add_feature(cfeature.BORDERS, linestyle='-', edgecolor='gray')
ax.coastlines(resolution='10m', color='black', linewidth=0.8)

shape.plot(ax=ax, facecolor='none', edgecolor='gray', linewidth=0.7)

filter_day = gdf_day_points[gdf_day_points['prec'] > 20]
filter_day.plot(
    ax=ax,
    marker='o',
    column='prec',
    cmap=cmap,
    legend=True,
    norm=norm,
    alpha=0.15,
    legend_kwds={'label': "Precipitation (mm)", 'orientation': "vertical"}
)

prec_filtered_day = prec_day[prec_day['prec'] > 20]
for idx, row in prec_filtered_day.iterrows():
    state = shape[shape['SIGLA_UF'] == row['SIGLA_UF']].geometry.iloc[0]
    point = state.representative_point()
    x, y = point.x, point.y
    ax.text(x, y, f'{row["prec"]:.1f}', fontsize=12, ha='center', color='darkred')

gl = ax.gridlines(draw_labels=True, color='gray', linestyle='--', linewidth=0.2)
gl.top_labels = False
gl.right_labels = False
gl.xlocator = mticker.FixedLocator(range(-180, 181, 5))
gl.ylocator = mticker.FixedLocator(range(-90, 91, 5))
gl.xlabel_style = {'size': 12}
gl.ylabel_style = {'size': 12}

plt.suptitle('MERGE DAILY', fontsize=20, fontweight='bold', y=1.00)
plt.title(f'Cumulative precipitation last pentad and average per state', fontsize=17, loc='left', x=0.25, y=1.06, fontweight='bold')

ax.text(0.001, 1.01, f'Pentad: ({pentadas[0]} - {pentadas[4]})', transform=ax.transAxes, fontsize=17, 
        ha='left', va='center', color='black')

ax.set_title('Valid time: ' + datetime.now().strftime('%d/%m/%Y'), fontsize=17, loc='right')

plt.savefig(f'./figures/DAILY_{datetime.now().strftime("%Y%m%d")}.png', dpi=300)

# Plot HOURLY
extent = [-80.0, -34.0, -45.0, 10.0]
fig, ax = plt.subplots(figsize=(20, 20), subplot_kw=dict(projection=ccrs.PlateCarree()))
ax.set_extent(extent, crs=ccrs.PlateCarree())
ax.add_feature(cfeature.BORDERS, linestyle='-', edgecolor='gray')
ax.coastlines(resolution='10m', color='black', linewidth=0.5)

shape.plot(ax=ax, facecolor='none', edgecolor='gray', linewidth=0.7)

gdf_filtered = gdf_hour_points[gdf_hour_points['prec'] > 20]
gdf_filtered.plot(
    ax=ax,
    marker='o',
    column='prec',
    cmap=cmap,
    markersize=20,
    norm=norm,
    alpha=0.15,
    legend=True,
    legend_kwds={'label': "Precipitation (mm)", 'orientation': "vertical"}
)

prec_filtered_hour = prec_hour[prec_hour['prec'] > 20]
for idx, row in prec_filtered_hour.iterrows():
    state = shape[shape['SIGLA_UF'] == row['SIGLA_UF']].geometry.iloc[0]
    point = state.representative_point()
    x, y = point.x, point.y
    ax.text(x, y, f'{row["prec"]:.1f}', fontsize=12, ha='center', color='darkred')

gl = ax.gridlines(draw_labels=True, color='gray', linestyle='--', linewidth=0.2)
gl.top_labels = False
gl.right_labels = False
gl.xlocator = mticker.FixedLocator(range(-180, 181, 5))
gl.ylocator = mticker.FixedLocator(range(-90, 91, 5))
gl.xlabel_style = {'size': 12}
gl.ylabel_style = {'size': 12}

plt.suptitle('MERGE HOURLY', fontsize=20, fontweight='bold',  y=1.00)
plt.title(f'Cumulative precipitation last pentad and average per state', fontsize=17, loc='left', x=0.18, y=1.06, fontweight='bold')

ax.text(0.001, 1.01, f'Pentad: ({pentadas[0]} - {pentadas[4]})', transform=ax.transAxes, fontsize=17, 
        ha='left', va='center', color='black')

ax.set_title('Valid time: ' + datetime.now().strftime('%d/%m/%Y'), fontsize=17, loc='right')

plt.savefig(f'./figures/HOURLY_{datetime.now().strftime("%Y%m%d")}.png', dpi=300)

# BIAS plot
norm_bias = mcolors.TwoSlopeNorm(vmin=-20, vcenter=0, vmax=20)

extent = [-80.0, -34.0, -45.0, 10.0]
fig, ax = plt.subplots(figsize=(20, 20), subplot_kw=dict(projection=ccrs.PlateCarree()))
ax.set_extent(extent, crs=ccrs.PlateCarree())
ax.add_feature(cfeature.BORDERS, linestyle='-', edgecolor='gray')
ax.coastlines(resolution='10m', color='black', linewidth=0.5)

shape.plot(ax=ax, facecolor='none', edgecolor='gray', linewidth=0.5)

for idx, row in diff.iterrows():
    state = shape[shape['SIGLA_UF'] == row['SIGLA_UF']].geometry.iloc[0]
    point = state.representative_point()
    x, y = point.x, point.y
    ax.text(x, y, f'{row["prec"]:.1f}', fontsize=12, ha='center', color='black')
    ax.add_geometries([state], crs=ccrs.PlateCarree(), 
                      facecolor=plt.cm.seismic(norm_bias(row['prec'])), alpha=0.8,
                      edgecolor='gray', linewidth=0.5)

cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm_bias, cmap='seismic'), ax=ax, orientation="vertical")
cbar.set_label("Bias - Precipitation (mm)")

gl = ax.gridlines(draw_labels=True, color='gray', linestyle='--', linewidth=0.2)
gl.top_labels = False
gl.right_labels = False
gl.xlocator = mticker.FixedLocator(range(-180, 181, 5))
gl.ylocator = mticker.FixedLocator(range(-90, 91, 5))
gl.xlabel_style = {'size': 12}
gl.ylabel_style = {'size': 12}

plt.suptitle('BIAS', fontsize=20, fontweight='bold', y=1.00)
plt.title(f'MERGE BIAS - Cumulative precipitation last pentad', fontsize=17, loc='left', x=0.25, y=1.06, fontweight='bold')

ax.text(0.001, 1.01, f'Pentad: ({pentadas[0]} - {pentadas[4]})', transform=ax.transAxes, fontsize=17,
        ha='left', va='center', color='black')

plt.savefig(f'./figures/BIAS_{datetime.now().strftime("%Y%m%d")}.png', dpi=300)


