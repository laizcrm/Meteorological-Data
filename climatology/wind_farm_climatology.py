"""
This script performs the monthly climatology analysis of wind farm generation based on wind speed data.
It compares daily observed wind speeds with the monthly climatological averages and forecasted values,
calculating the anomaly percentages relative to the climatology. The script then generates plots for each
wind farm complex showing the realized and forecasted wind speed anomalies as a percentage of the climatology.

Main steps:
- Load monthly climatology data from BigQuery.
- Load daily observed wind speed data.
- Load ECMWF forecast wind speed data.
- Match forecast grid points to wind farm locations and compute weighted averages by distance.
- Calculate anomalies of realized and forecast wind speed relative to climatology.
- Generate and save time series plots of anomalies for each wind farm complex.
"""


# Custom tools path (adjust or remove if not available)
sys.path.insert(1, '/mnt/c/scripts/tools/')
from bigquery_bucket_tools import *
from plot_tools import *
from tools import *
#-----------------------------------------------------------------------------------------------------------

# Parse the reference date from the command line argument
reference_date = datetime.strptime(sys.argv[1], '%Y-%m-%d')

# Format date string for file naming and queries
date_str = reference_date.date().strftime("%Y-%m-%d")

# Directory to save output plots
output_dir = '/mnt/c/scripts/daily_reports/wind_ne/results/'

# BigQuery parameters for climatology data
project_id = 'modeling-project'
dataset = 'climatologies'
table_climatology = 'era5_wind_monthly_climatology_for_projects'

# Initialize BigQuery client
bq_client = bigquery.Client(project=project_id)

# Define the months range to query climatology (20 days before and 15 days after reference date)
months_range = [
    date_range(reference_date, 'day', '-', 20).date().strftime("%m"),
    date_range(reference_date, 'day', '+', 15).date().strftime("%m")
]
months_range = [str(int(m)) for m in months_range]

# Query monthly climatology for wind speed
climatology_month = query_all_month(bq_client, project_id, dataset, table_climatology, months_range)

# Fix project names if needed (example: replace 'UMR' with 'UMA')
climatology_month['project'] = climatology_month['project'].replace('UMR', 'UMA')

# BigQuery parameters for daily observed wind data
project_id = 'meteorology-series'
dataset = 'era5_data_curated'
table_daily_wind = 'weighted_wind_daily_data_for_projects'

bq_client = bigquery.Client(project=project_id)

# Define date range for observed wind data (21 days before to 1 day before reference date)
start_date = date_range(reference_date, 'day', '-', 21).date().strftime("%Y-%m-%d")
end_date = date_range(reference_date, 'day', '-', 1).date().strftime("%Y-%m-%d")

# Query observed daily wind data
observed_data = query_between_date(bq_client, project_id, dataset, table_daily_wind, start_date, end_date, 'date_ref')

# Standardize project names
observed_data['project'] = observed_data['project'].replace('UMR', 'UMA')

# Extract month number for merging with climatology
observed_data['month'] = observed_data['date_ref'].apply(lambda x: int(pd.to_datetime(x).month))

# Merge observed data with climatology by project and month
merged_observed = observed_data.merge(climatology_month, on=['project', 'month'], how='inner', suffixes=['_observed', '_climatology'])

# BigQuery parameters for ECMWF wind forecast data
project_id = 'short-term-forecast'
dataset = 'processed_input_data'
table_forecast = 'ecmwf_processed_data'

bq_client = bigquery.Client(project=project_id)

# Query forecast data for the reference date
forecast_data = query_equal_date_list_select(
    bq_client,
    project_id,
    dataset,
    table_forecast,
    date_str,
    'model_time',
    ['project', 'model_time', 'timestamp', 'latitude', 'longitude', 'WS100']
)

# Process forecast data: drop model_time, group by lat/lon/project and day, take mean
forecast_data = (
    forecast_data
    .drop('model_time', axis=1)
    .set_index(['timestamp'])
    .groupby(['latitude', 'longitude', 'project', pd.Grouper(freq='1D')])
    .mean()
    .reset_index()
    .sort_values(by=['timestamp', 'project', 'latitude', 'longitude'], ascending=True)
    .rename({'timestamp': 'date_ref'}, axis=1)
    .reset_index(drop=True)
)

# BigQuery parameters for closest ECMWF grid points per project
project_id = 'modeling-project'
dataset = 'metadata'
table_closest_points = 'closest_four_points_of_each_complex_ecmwf_wind_high_resolution'

bq_client = bigquery.Client(project=project_id)

# Query grid points data
closest_points = query_all_date(bq_client, project_id, dataset, table_closest_points)

# Rename columns to match forecast data for merging
closest_points = closest_points.rename({'latitude_ecmwf': 'latitude', 'longitude_ecmwf': 'longitude'}, axis=1)

# Standardize project names
closest_points['project'] = closest_points['project'].replace('UMR', 'UMA')

# Merge forecast data with closest points info
merged_forecast = forecast_data.merge(closest_points, on=['latitude', 'longitude', 'project'], how='outer')

# Drop unnecessary columns
merged_forecast = merged_forecast.drop(merged_forecast.columns[[0, 1, 5, 6]], axis=1)

# Calculate weighted average of wind speed by distance
merged_forecast = weighted_average(merged_forecast, 'WS100', 'dist_km', ['date_ref', 'project']).rename(columns={'project': 'project'})

# Extract month for climatology merge
merged_forecast['month'] = merged_forecast['date_ref'].apply(lambda x: int(pd.to_datetime(x).month))

# Merge forecast data with climatology
merged_forecast_final = merged_forecast.merge(climatology_month, on=['project', 'month'], how='inner', suffixes=['_forecast', '_climatology'])

# Calculate anomaly (%) relative to climatology for observed and forecast data
merged_observed['anomaly'] = (merged_observed['ws100_avg_observed'] / merged_observed['ws100_avg_climatology']) * 100
merged_forecast_final['anomaly'] = (merged_forecast_final['ws_pond'] / merged_forecast_final['ws100_avg']) * 100

# Google Cloud Storage setup for saving plots
storage_client = storage.Client(project='modeling-project')
bucket_name = 'meteorological_reports'

# Plotting parameters
palette = sns.color_palette("tab10")

num_projects = len(merged_observed['project'].unique())
cols = 3
rows = -(-num_projects // cols)  # ceiling division for rows

print('Starting plot generation...')

# Convert date columns to datetime for plotting
merged_observed['date_ref'] = pd.to_datetime(merged_observed['date_ref'])
merged_forecast_final['date_ref'] = pd.to_datetime(merged_forecast_final['date_ref'])

# Loop over each project to generate plots
for idx, project in enumerate(merged_observed['project'].unique()):
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    print(f'Processing project: {project}')
    
    ax.set_title(f'Complex: {project}', y=1.13, fontsize=20, weight='bold')
    
    # Select data for the current project
    df_obs = merged_observed[merged_observed['project'] == project]
    df_fcst = merged_forecast_final[merged_forecast_final['project'] == project]
    
    # Pivot observed anomaly for plotting
    obs_plot = df_obs.pivot(index='date_ref', columns='project', values=['anomaly'])
    obs_plot.plot(ax=ax, linewidth=3, color=palette[idx])
    
    # Pivot forecast anomaly for plotting
    fcst_plot = df_fcst.pivot(index='date_ref', columns='project', values=['anomaly'])
    fcst_plot.plot(ax=ax, linewidth=3,
