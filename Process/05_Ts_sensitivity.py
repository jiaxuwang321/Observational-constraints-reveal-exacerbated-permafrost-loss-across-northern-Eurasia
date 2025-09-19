###
    # This code calculates the average soil temperature sensitivity of observations and various models on seasonal scale in the northern Eurasia.
    # 1. Read soil temperature data includes (1)observations:"../Data/ALLstations.nc" (2)models:"../Data/cmip6_tsl/*.nc"
    # 2. Calculate the seasonal soil temperature sensitivity of each site for each year.
    # 3. Retain sites that more than 16 years of observation (1982-2014 total: 33 years) and within or at the boundaries of permafrost regions.
    # 4. Calculate the annual-average sensitivity of each site.
    # 5. Use bilinear interpolation method to interpolate the simulated soil temperature to station.
    # 6. Use nearest neighbor interpolation method to ensure non missing measurements for stations on land edge.
    # 7. Calculate the spatial average of all sites - grid processing before average - to reduce the uneven impact of sites.
    # 8. Output: "../Data/Ts_seasonal_sensitivity/Ts{depth}_Ts02_all_station_avg_obs_18model.txt", as the x-axis of Figure 3
    # 9. Calculate the simulated Ts Sensitivity of Each Grid ——— tsl_seasonalcycle_avg (lat,lon)
###

import xarray as xr
import glob
import math
import os
import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
import warnings
warnings.filterwarnings("ignore")

## -----------------------------------------------------------------------------------------------------------------------------
    #                                Functions used in this code
## -----------------------------------------------------------------------------------------------------------------------------

# Find the permafrost probability of the grid cell nearest to the station
def find_nearest_extent(lat, lon, per_extent, station_lat, station_lon):
    lat_idx = np.abs(lat - station_lat).argmin()
    lon_idx = np.abs(lon - station_lon).argmin()
    return per_extent[lat_idx, lon_idx]

# Calculate the multi-year average soil temperature sensitivity at stations within or at the boundary of the permafrost region
def process_file(Ts_Ts02, lat, lon, per_extent, boundary_stations):
    
    # Filter stations with more than 16-years observation (1982-2014 total:33 years) for analysis
    filtered_df = Ts_Ts02.groupby('station').filter(lambda x: len(x) > 16)
    
    # Calculate the permafrost probability corresponding to station locations 
    # (using nearest neighbor interpolation; consideration will be given to stations at boundary locations in subsequent steps)
    nearest_extent_values = [find_nearest_extent(lat, lon, per_extent, row['Lat'], row['Lon']) for _, row in filtered_df.iterrows()]
    
    # Select stations within the permafrost extent
    pf_stations = filtered_df[np.array(nearest_extent_values) > 0]
    
    # Add Stations located at the boundary of permafrost area
    boundary_to_add = filtered_df[filtered_df['station'].isin(boundary_stations)]

    # Concat stations within the permafrost extent and boundary stations
    pf_stations = pd.concat([pf_stations, boundary_to_add[~boundary_to_add['station'].isin(pf_stations['station'])]])
    
    return pf_stations

# Searches for non-missing values near a given location
def get_nearest_nonmissing_value(tsl_data, time_idx, depth_idx, lat_idx, lon_idx, max_radius=3):
    """
    Searches for the nearest non-missing value within a maximum radius `max_radius` around (lat_idx, lon_idx) 
    in tsl_data[time_idx, depth_idx, :, :].
    If none is found, returns the original value (which may be NaN).
    """
    # Get the original value
    orig_val = tsl_data[time_idx, depth_idx, lat_idx, lon_idx]
    if not np.isnan(orig_val):
        return orig_val  # The nearest point is already valid, no need to search

    # Initialize the best candidate value and distance
    best_value = np.nan
    best_distance = np.inf
    found = False

    # Get the total number of latitude and longitude points (to ensure the search stays within bounds)
    nlat = tsl_data.shape[2]
    nlon = tsl_data.shape[3]

    # Search the surrounding area sequentially with increasing radius r
    for r in range(1, max_radius + 1):
        # Define the search range, with boundary checks
        lat_min = max(0, lat_idx - r)
        lat_max = min(nlat - 1, lat_idx + r)
        lon_min = max(0, lon_idx - r)
        lon_max = min(nlon - 1, lon_idx + r)
        # Iterate over all grid points within this range
        for i in range(lat_min, lat_max + 1):
            for j in range(lon_min, lon_max + 1):
                candidate = tsl_data[time_idx, depth_idx, i, j]
                if not np.isnan(candidate):
                    # Calculate the Euclidean distance (in grid cell units)
                    dist = np.sqrt((i - lat_idx)**2 + (j - lon_idx)**2)
                    # Select the candidate with the smaller distance
                    if dist < best_distance:
                        best_distance = dist
                        best_value = candidate
                        found = True
        if found:
            # Once a valid value is found within a certain radius, return the closest one found
            return best_value
    # If no valid value is found within the maximum radius, return the original value
    return best_value

# Calculate the grid areas
def area(lat, grid_resolution):    
    rad = np.pi / 180.0
    lat_rad = lat * rad

    # Calculate the latitude and longitude resolution
    re = 6371220.0  # Earth radius (in meters)
    dlon = grid_resolution * rad  # Longitude resolution (in radians)
    dlat = grid_resolution * rad  # Latitude resolution (in radians)

    # Calculate the area of each grid cell
    sin_lat1 = np.sin(lat_rad + dlat / 2)
    sin_lat2 = np.sin(lat_rad - dlat / 2)
    grid_areas = (re ** 2) * dlon * (sin_lat1 - sin_lat2)
    return grid_areas

# Grids the station data and calculates grid areas to facilitate subsequent area-weighted averaging
def grid_data(data, grid_resolution):
    grid_Mean = pd.DataFrame(columns=data.columns[2:])
    grid_Mean['area']=''
    latmin, latmax = data['Lat'].min(), data['Lat'].max()
    m = math.ceil((latmax - latmin) / grid_resolution)
    line = 0
    for j in range(m):
        data1 = data[(data['Lat'] >= latmin + grid_resolution * j) & (data['Lat'] < latmin + grid_resolution * (j + 1))]
        if data1.empty:
            continue
        minlon, maxlon = data1['Lon'].min(), data1['Lon'].max()
        n = math.ceil((maxlon - minlon) / grid_resolution)
        for i in range(n):
            data2 = data1[(data1['Lon'] >= minlon + grid_resolution * i) & (data1['Lon'] < minlon + grid_resolution * (i + 1))]
            grid_Mean.loc[line] = data2.drop(columns=['station', 'year']).mean()
            grid_Mean['area'].loc[line] = area((latmin + grid_resolution * (j+0.5)),grid_resolution)
            line += 1
    return grid_Mean.dropna()

# This function is to check 1.6/2.4/3.2mTs have no missing measurements or have missing measurements in January or June or December 
# (the annual minimum temperature of deep soil usually does not occur in these months), execute subsequent commands
def is_column_valid(col):
    na_mask = col.isna()
    if not na_mask.any():  # No missing values
        return True
    na_indices = na_mask[na_mask].index
    # Check if missing values only appear at the beginning or end
    valid_positions = {0, 6, len(col)-1}
    return all(idx in valid_positions for idx in na_indices)

# This function is to check 0.2m Ts, if all 0.2m Ts are not missing measurements, execute subsequent commands
def is_column_no_missing(col):
    # Check if the column has absolutely no missing values
    return not col.isna().any()

ds = xr.open_dataset('../Data/ALLstations.nc')
longitude, latitude = ds["lon"][:], ds["lat"][:]
tsoil = ds.tsoil.sel(year=slice('1982', '2014'))  # tsoil (station, year, month, depth)

model_order = ["CESM2","CESM2-FV2","CESM2-WACCM","CNRM-CM6-1-HR","CNRM-CM6-1","CNRM-ESM2-1","E3SM-1-1",
               "FGOALS-f3-L","FGOALS-g3","GFDL-CM4","GFDL-ESM4","MIROC6","MIROC-ES2L","MPI-ESM1-2-HR",
               "MPI-ESM1-2-LR","NorESM2-LM","NorESM2-MM","TaiESM1","obs_Ts_Ts02"]

for depth in [320,160,240]:

    ## -----------------------------------------------------------------------------------------------------------------------------
    #                       Calculate the soil temperature sensitivity of each site for each year
    ## -----------------------------------------------------------------------------------------------------------------------------
    
    q = xr.open_dataset('../Data/Probability_pf/probability_lt_threshold_'+str(depth)+'.nc')
    lat, lon = q['lat'][:], q['lon'][:]
    per_extent = q['probability_lt_threshold'][:,:]

    if depth == 160:
        depindex = 5 # The fifth layer of tsoil
        boundary_stations = ['STN0008','STN0009','STN0014','STN0016','STN0019','STN0020','STN0051','STN0052','STN0056','STN0063','STN0068','STN0076','STN0077','STN0232','STN0235','STN0238','STN0243','STN0250','STN0437','STN0440','STN0442','STN0305','STN0309','STN0304','STN0355','STN0356','STN0357','STN0351','STN0346','STN0314','STN0113','STN0261']
    elif depth == 240:
        depindex = 6 # The sixth layer of tsoil
        boundary_stations = ['STN0008','STN0009','STN0014','STN0016','STN0019','STN0020','STN0051','STN0052','STN0056','STN0063','STN0068','STN0076','STN0077','STN0232','STN0235','STN0238','STN0243','STN0250','STN0437','STN0440','STN0442','STN0305','STN0309','STN0304','STN0355','STN0356','STN0357','STN0358','STN0351','STN0346','STN0113','STN0013','STN0014','STN0018','STN0017','STN0061','STN0074','STN0177','STN0233','STN0234','STN0236','STN0261','STN0278','STN0440','STN0436']
    else:
        depindex =7 # The seventh layer of tsoil
        boundary_stations = ['STN0016','STN0019','STN0020','STN0051','STN0052','STN0056','STN0063','STN0068','STN0076','STN0077','STN0192','STN0232','STN0235','STN0238','STN0243','STN0250','STN0292','STN0437','STN0442','STN0305','STN0309','STN0304','STN0355','STN0356','STN0357','STN0358','STN0351','STN0346','STN0113','STN0018','STN0017','STN0074','STN0177','STN0233','STN0234','STN0236','STN0261','STN0278','STN0440','STN0352','STN0390','STN0240','STN0258']
    
    df = tsoil[:,:,:,depindex].to_dataframe(name='tsoil').reset_index() # soil temperature at 1.6m/2.4m/3.2m
    df2 = tsoil[:,:,:,0].to_dataframe(name='tsoil02').reset_index() # soil temperature at 0.2m
    df['tsoil02'] = df2['tsoil02']

    # Prepare a DataFrame to store other stations
    Ts_Ts02 = pd.DataFrame(columns=['station', 'year', 'Lat', 'Lon', 'obs_Ts_Ts02'])
    df_grouped = df.groupby('station')
    s = 0 # Index of each station
    h = 0 # Index of Ts_Ts02's line numbers
    for station, station_df in df_grouped:
        station_df['lat'] = latitude[s].data
        station_df['lon'] = longitude[s].data
        s += 1

        # Group by year and Calculate seasonal Ts sensitivity (relative to Ts02) for each group
        station_df_grouped = station_df.sort_values(by=['year', 'month']).groupby('year')
        # Calculate seasonal Ts sensitivity for each year
        for year, annual_df in station_df_grouped:
            annual_df.index = range(len(annual_df))

            # Check 1.6/2.4/3.2m Ts and 0.2m Ts for each year
            col4 = annual_df.iloc[:, 4] # 1.6/2.4/3.2m Ts
            col5 = annual_df.iloc[:, 5] # 0.2m Ts
            if is_column_valid(col4) and is_column_no_missing(col5):
                # Store station information
                Ts_Ts02.loc[h, 'station'] = annual_df.station.iloc[0]
                Ts_Ts02.loc[h, 'year'] = annual_df.year.mean()
                Ts_Ts02.loc[h, 'Lat'] = annual_df.lat.mean()
                Ts_Ts02.loc[h, 'Lon'] = annual_df.lon.mean()

                # Calculate observed Ts/Ts02 ratio: seasonal cycle amplitude of deep-layer Ts relative to Ts02
                A1 = (max(annual_df.iloc[:, 4]) - min(annual_df.iloc[:, 4])) / (max(annual_df.iloc[:, 5]) - min(annual_df.iloc[:, 5]))
                Ts_Ts02.loc[h, 'obs_Ts_Ts02'] = A1
                h += 1

    ## -----------------------------------------------------------------------------------------------------------------------------
    #     Retain stations that more than 16 years of observation (1982-2014 total: 33 years) and within or at the 
    #     boundaries of permafrost regions, and calculate the multi-year mean of the sensitivity at each station
    ## -----------------------------------------------------------------------------------------------------------------------------
    
    pf_stations = process_file(Ts_Ts02, lat, lon, per_extent, boundary_stations)
    # Calculate the multi-year average soil temperature sensitivity for each station
    PF_stations_mean_pf = pf_stations.groupby('station').mean().reset_index().sort_values(by='station') 
    
    # Index monthly sensitivity observations by station number, to subsequently match monthly sensitivity values from models
    pf_stations_filtered = pf_stations[pf_stations['station'].isin(PF_stations_mean_pf['station'])] 
    pf_stations_filtered.index = range(len(pf_stations_filtered))
    n = len(pf_stations_filtered)

    ## -----------------------------------------------------------------------------------------------------------------------------
    #        Use bilinear interpolation method to interpolate the simulated soil temperature to station
    ## -----------------------------------------------------------------------------------------------------------------------------

    # Match and filter the monthly tsl simulated by the corresponding mode of the backend station
    Model_Soil_df = pd.DataFrame(index=range(n * 12),
                                columns=['station', 'year', 'month'])
    file_list = glob.glob("../Data/cmip6_tsl/*.nc")

    # Prepare a DataFrame to store the tsl of CMIP6 models 
    model_data = {}
    for tsl_file in file_list:
        modelname = tsl_file[18:-22]
        Model_Soil_df[modelname+'_tsl32'] = ''
        Model_Soil_df[modelname+'_tsl02'] = ''
        model_data[modelname] = xr.open_dataset(tsl_file, decode_times=False)

        # Loop through each station and each year, reserving 12 rows of space for monthly soil temperature data for each station per year
    for i in range(n):
        start = i * 12
        end = start + 11 
        # Assign the station and year values from each row of pf_stations_filtered to Model_Soil_df
        Model_Soil_df['station'].loc[start:end] = pf_stations_filtered.loc[i, 'station']
        Model_Soil_df['year'].loc[start:end] = pf_stations_filtered.loc[i, 'year']
        # Assign month values, generating an integer sequence from 1 to 12
        Model_Soil_df['month'].loc[start:end] = np.arange(1, 13)
        for modelname, f in model_data.items():
            tsl = f['tsl'][:, :, :, :].data  # tsl(time, depth, lat, lon)
            model_lat = f['lat'].data  
            model_lon = f['lon'].data    
            # Initialize arrays to store interpolation results
            tsl32_values = np.empty(12)  # Length of time dimension
            tsl02_values = np.empty(12)

            for t in range(12):
                interpolator_tsl32 = RegularGridInterpolator(
                    (model_lat, model_lon),  # Grid points of latitude and longitude
                    tsl[int((pf_stations_filtered.loc[i, 'year']-1982)*12+t), int(depth/10)-1, :, :], # Match 12 months of 3.2/2.4/1.6m soil temperature data for the station's corresponding year
                    bounds_error=False, fill_value=np.nan
                )
                interpolator_tsl02 = RegularGridInterpolator(
                    (model_lat, model_lon),  # Grid points of latitude and longitude
                    tsl[int((pf_stations_filtered.loc[i, 'year']-1982)*12+t), 1, :, :], # Match 12 months of 0.2m soil temperature data for the station's corresponding year
                    bounds_error=False, fill_value=np.nan
                )
                # Perform interpolation based on station coordinates
                tsl32_values[t] = interpolator_tsl32([pf_stations_filtered.loc[i, 'Lat'], pf_stations_filtered.loc[i, 'Lon']])
                tsl02_values[t] = interpolator_tsl02([pf_stations_filtered.loc[i, 'Lat'], pf_stations_filtered.loc[i, 'Lon']])
            # Convert results to Celsius, and place the 12 monthly values of 0.2m and 3.2/2.4/1.6m soil temperature for the year into Model_Soil_df
            Model_Soil_df[modelname + '_tsl32'].loc[start:end] = tsl32_values - 273.15
            Model_Soil_df[modelname + '_tsl02'].loc[start:end] = tsl02_values - 273.15

    ## -----------------------------------------------------------------------------------------------------------------------------
    #       Some stations at land-sea boundaries have missing soil temperature values in the model, 
    #       using nearest neighbor interpolation to obtain soil temperature from corresponding land grids
    ## -----------------------------------------------------------------------------------------------------------------------------
    
    df_missing = Model_Soil_df[Model_Soil_df.isna().any(axis=1)] # Filter rows containing missing values
    df_missing.index = range(len(df_missing))

    # For each model's corresponding two columns (first column is station, second is year, third is month, so starting from column 4, every two columns represent one model's _tsl32 and _tsl02)
    for i in range(3, len(df_missing.columns) - 1, 2):
        col_name_soil32 = df_missing.columns[i]
        col_name_soil02 = df_missing.columns[i + 1]
        # Find the row indices where missing values occur in this column (soil32 and soil02 will have missing values simultaneously, sharing the same row indices)
        missing_indices = df_missing.index[df_missing[col_name_soil32].isna()].tolist()
        if not missing_indices:
            continue
        else:
            print(f"Column {col_name_soil32} missing values at row indices: {missing_indices}")
            # Construct the nc file path based on part of the column name
            tsl_file = glob.glob("../Data/cmip6_tsl/" + col_name_soil32[:-6] + "_*.nc")[0]
            f = xr.open_dataset(tsl_file, decode_times=False)
            tsl = f['tsl'][:, :, :, :].data  # (time, depth, lat, lon)

            for row in missing_indices:
                # Get corresponding latitude and longitude based on station number
                station = df_missing.loc[row, 'station']
                station_info = PF_stations_mean_pf[PF_stations_mean_pf['station'] == station]
                station_lat = station_info['Lat'].values[0]
                station_lon = station_info['Lon'].values[0]

                # Begin nearest neighbor interpolation, first find the nearest grid index
                lat_j = np.argmin(np.abs(f["lat"].data - station_lat))
                lon_j = np.argmin(np.abs(f["lon"].data - station_lon))
                # Find time index based on year and month (starting year is 1982)
                time_idx = int((df_missing.loc[row, 'year'] - 1982) * 12 + df_missing.loc[row, 'month'] - 1)
                # Match tsl at 3.2/2.4/1.6m depth
                raw_val32 = tsl[time_idx, int(depth/10)-1, lat_j, lon_j]

                if np.isnan(raw_val32):
                    # Search for the nearest non-missing value
                    raw_val32 = get_nearest_nonmissing_value(tsl, time_idx, int(depth/10)-1, lat_j, lon_j, max_radius=3)
                df_missing.loc[row, col_name_soil32] = raw_val32 - 273.15

                # Match tsl at 0.2m depth
                raw_val02 = tsl[time_idx, 1, lat_j, lon_j]
                if np.isnan(raw_val02):
                    # Search for the nearest non-missing value
                    raw_val02 = get_nearest_nonmissing_value(tsl, time_idx, 1, lat_j, lon_j, max_radius=3)
                df_missing.loc[row, col_name_soil02] = raw_val02 - 273.15

    # Merge df_missing with Model_Soil_df, replacing the previously missing rows in Model_Soil_df, and store the new Model_Soil_df
    for idx, row in df_missing.iterrows():
        # Construct matching condition: first three columns are identical
        condition = (
            (Model_Soil_df['station'] == row['station']) &
            (Model_Soil_df['year'] == row['year']) &
            (Model_Soil_df['month'] == row['month'])
        )
        # Get the list of indices that satisfy the condition
        matching_indices = Model_Soil_df.index[condition]
        if not matching_indices.empty:
            # Copy row.values to a 2D array with the same number of rows as matching rows
            replacement = pd.DataFrame(
                [row.values] * len(matching_indices),
                columns=Model_Soil_df.columns,
                index=matching_indices
            )
            # Replace all matching rows at once
            Model_Soil_df.loc[matching_indices, :] = replacement
    ## -----------------------------------------------------------------------------------------------------------------------------
    #       Calculate the seasonal Ts sensitivity (relative to Ts02) for each station 
    #       simulated by each model based on soil temperature data
    ## -----------------------------------------------------------------------------------------------------------------------------
    for tsl_file in file_list:
        modelname = tsl_file[18:-22]
        pf_stations_filtered[modelname] = '' # Create a new column for each model to store sensitivity values
    df_grouped = Model_Soil_df.groupby('station')

    h=0
    for station, station_df in df_grouped:
        # Group by year and calculate Ts/Ts02 for each group
        station_df_grouped = station_df.sort_values(by=['year', 'month']).groupby('year')
        for year, annual_df in station_df_grouped:
            # Calculate Ts/Ts02 ratio for each model, one value per year
            for m, tsl_file in enumerate(file_list):
                modelname = tsl_file[18:-22]
                n = m * 2 + 3  # Corresponding column indice of Ts, the column indice of Ts02 is n+1
                A1 = (max(annual_df.iloc[:, n]) - min(annual_df.iloc[:, n])) / (max(annual_df.iloc[:, n + 1]) - min(annual_df.iloc[:, n + 1]))
                pf_stations_filtered.loc[h, modelname] = A1
            h += 1

    pf_stations_filtered.to_csv(f"../Data/Ts_seasonal_sensitivity/Ts_Ts02_obs_NH45_18model_{depth}cm_interpolated.csv")
    
    # Calculate the multi-year average soil temperature sensitivity for each station
    pf_stations_years_ave = pf_stations_filtered.groupby('station').mean().reset_index().sort_values(by='station')

    ## -----------------------------------------------------------------------------------------------------------------------------
    #      Perform gridding processing on soil temperature sensitivity for all stations, 
    #      then calculate spatial average (including one value for observations and each model respectively)
    ## -----------------------------------------------------------------------------------------------------------------------------
    grid_resolution = 2.5
    n_bootstrap = 1000  # Number of sampling iterations
    station_means = []

    # Bootstrap sampling method
    for _ in range(n_bootstrap):
        # Sample len(pf_stations_years_ave) stations with replacement
        sampled_df = pf_stations_years_ave.sample(n=len(pf_stations_years_ave), replace=True)
        
        # Calculate the station average for this sampling iteration (each sampling undergoes gridding processing first, then calculates area-weighted average after gridding)
        grid_station = grid_data(sampled_df, grid_resolution)
        weighted_avg_models = {
            col: (grid_station[col] * grid_station["area"]).sum() / grid_station["area"].sum()
            for col in grid_station.columns[2:-1]
        }
        station_means.append(weighted_avg_models)
    station_means = pd.DataFrame(station_means)

    # Calculate the standard deviation and mean of the 1000 sampling results
    std_dev = station_means.std().obs_Ts_Ts02
    mean_series = station_means.mean().reindex(model_order)
    b = pd.concat([mean_series, pd.Series([std_dev])], ignore_index=True)
    print(b)
    output_file = f"../Data/Ts_seasonal_sensitivity/Ts{depth}cm_Ts02_all_station_avg_obs_18model.txt"
    b.to_csv(output_file, index=False, header=False)

## -----------------------------------------------------------------------------------------------------------------------------
#                        Calculate the simulated Ts Sensitivity of Each Grid 
## -----------------------------------------------------------------------------------------------------------------------------
# Define models and depth levels
models = [
    "CESM2", "CESM2-FV2", "CESM2-WACCM", "CNRM-CM6-1-HR", "CNRM-CM6-1", "CNRM-ESM2-1", 
    "E3SM-1-1", "FGOALS-f3-L", "FGOALS-g3", "GFDL-CM4", "GFDL-ESM4", "MIROC6", "MIROC-ES2L", 
    "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "NorESM2-LM", "NorESM2-MM", "TaiESM1"
]
years = 2014-1982+1

# Main processing loop
for modelname in models:
    # Search for the model file
    alt_file = glob.glob(f"/home/wangjx/Data/cmip6_tsl_05_depth_interpolated/{modelname}_*.nc")[0]

    # Open dataset
    with xr.open_dataset(alt_file, decode_times=False) as f:
        tsl = f['tsl'][:, :, :, :360]  # Limit to 360 longitude
        latitude = f['lat'][:]
        longitude = f['lon'][:360]

        # Reshape tsl data for easy access
        tsl_reshaped = tsl.data.reshape((years, 12, *tsl.shape[1:]))

        # Iterate over each depth level
        for depth in [1.6, 2.4, 3.2]:
            depth_index = int(depth * 10 - 1)

            # Initialize result arrays
            tsl_seasonalcycle_avg = np.full((len(latitude), len(longitude)), np.nan)
            
            # Calculate seasonal cycle sensitivity
            for j in range(len(latitude)):
                for k in range(len(longitude)):
                    Ts_max = np.max(tsl_reshaped[:, :, depth_index, j, k], axis=1)
                    Ts_min = np.min(tsl_reshaped[:, :, depth_index, j, k], axis=1)
                    Ts20_max = np.max(tsl_reshaped[:, :, 1, j, k], axis=1)
                    Ts20_min = np.min(tsl_reshaped[:, :, 1, j, k], axis=1)
                    
                    # Calculate yearly sensitivity and average it
                    tsl_seasonalcycle_yearly = (Ts_max - Ts_min) / (Ts20_max - Ts20_min)
                    avg_sensitivity = np.mean(tsl_seasonalcycle_yearly)
                    
                    # Store the result, applying threshold
                    tsl_seasonalcycle_avg[j, k] = avg_sensitivity if avg_sensitivity <= 3 else np.nan
            
            # Create the xarray dataset
            ds = xr.Dataset(
                {'tsl_seasonalcycle_avg': (['lat', 'lon'], tsl_seasonalcycle_avg)},
                coords={'lat': latitude, 'lon': longitude}
            )
            ds.attrs['description'] = 'Seasonal cycle average of TSL sensitivity values'
            ds.attrs['model'] = modelname

            # Define output directory and save
            output_dir = f'../Data/Ts_seasonal_sensitivity/Simulated_grid_sensitivity/{int(depth*100)}/'
            os.makedirs(output_dir, exist_ok=True)
            output_file = os.path.join(output_dir, f'{modelname}_tsl_seasonalcycle_avg.nc')
            ds.to_netcdf(output_file)
            print(f"Saved: {output_file}")
