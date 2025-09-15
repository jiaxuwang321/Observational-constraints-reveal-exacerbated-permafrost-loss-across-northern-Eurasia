###
    # This code is for calculating the ALT sensitivity on two timesacles
    # (a)seasonal-cycle, and under global warming levels of (b)1.5°C, (c)2°C and (d)3°C relative to pre-industrial period (1850-1900) under SSP5-8.5 scenario.
    # Read tsl data from each model interpolated to 0.5° resolution (since subsequent calculations require averaging sensitivity across grid points, uniform resolution is needed)
    # Calculate the ratio of ALT change to 0.2 m soil temperature change for each grid.
    # Store the results in "../Data/ALT_Sensitivity_lat_lon/" for plotting Figure S4.
###

import pandas as pd
import numpy as np
import xarray as xr
from netCDF4 import Dataset
import glob
import warnings
warnings.filterwarnings("ignore")

models = ["CESM2","CESM2-FV2","CESM2-WACCM","CNRM-CM6-1-HR","CNRM-CM6-1","CNRM-ESM2-1",
          "E3SM-1-1","FGOALS-f3-L","FGOALS-g3","GFDL-CM4","GFDL-ESM4","MIROC6","MIROC-ES2L",
          "MPI-ESM1-2-HR","MPI-ESM1-2-LR","NorESM2-LM","NorESM2-MM","TaiESM1"]
models = ["TaiESM1"]
## -----------------------------------------------------------------------------------------------------------------------------
    #                        Calculate seasonal-cycle of ALT sensitivity 
## -----------------------------------------------------------------------------------------------------------------------------
begin = 1982 
end = 2014
years = end-begin+1
"""
for modelname in models:
    # Open TSL dataset
    tsl_file = glob.glob(f"/home/wangjx/Data/cmip6_tsl_05_depth_interpolated_198010/{modelname}_*.nc")[0]
    f = xr.open_dataset(tsl_file, decode_times=False)
    tsl = f['tsl'][12*(begin-1980):12*(begin-1980)+12*years,:,:,:360] # 1982-2014 monthly lon:0-180°
    latitude = f['lat'][:]
    longitude = f['lon'][:360]

    # Get dimensions
    months, depths, latitudes, longitudes = tsl.shape

    # Reshape TSL data
    tsl_reshaped = tsl.data.reshape((years, 12, depths, latitudes, longitudes))

    # Open ALT dataset
    filename1 = glob.glob("../Data/cmip6_alt_m15_NH45_1850_2100/"+modelname+"_alt*hd.nc")[0]# 1850-2100 monthly
    ds = Dataset(filename1)
    alt = ds["alt"][12*(begin-1850):(12*(begin-1850)+12*years), :, :360]  # 1982-2014  lon:0-180°
    alt[alt == -9999] = np.nan
    alt_reshaped = np.reshape(alt, (years, 12, latitudes, longitudes))
    alt_reshaped = alt_reshaped.filled(np.nan)
    # Initialize arrays for results
    tsl_seasonalcycle_avg = np.full((latitudes, longitudes), np.nan)

    # Vectorized computation of seasonal cycle averages
    for j in range(latitudes):
        for k in range(longitudes):
            alt_max = np.max(alt_reshaped[:, :, j, k], axis=1)
            alt_min = np.min(alt_reshaped[:, :, j, k], axis=1)

            # Avoid division by zero
            Ts20_max = np.max(tsl_reshaped[:, :, 1, j, k], axis=1)
            Ts20_min = np.min(tsl_reshaped[:, :, 1, j, k], axis=1)

            alt_max_avg = np.max(alt_max)
            tsl_seasonalcycle_yearly = (alt_max - alt_min) / (Ts20_max - Ts20_min)
            tsl_seasonalcycle_avg[j, k] = np.mean(tsl_seasonalcycle_yearly)

    # Create and save the dataset
    ds = xr.Dataset(
        {'alt_seasonalcycle_avg': (['lat', 'lon'], tsl_seasonalcycle_avg)},
        coords={'lat': latitude, 'lon': longitude}
    )

    # Add metadata
    ds.attrs['description'] = 'Seasonal cycle average of ALT sensitivity'
    ds.attrs['model'] = modelname

    # Save to NetCDF
    output_file = f"../Data/ALT_Sensitivity_lat_lon/seasonal_scale/{modelname}_alt_seasonalcycle_avg.nc"
    ds.to_netcdf(output_file)
"""
## -----------------------------------------------------------------------------------------------------------------------------
    #                     Calculate ALT sensitivity under three global warming levels
## -----------------------------------------------------------------------------------------------------------------------------

model_year = pd.read_csv("../Data/Tas_data/years_reaching_three_warming_thresholds_ssp585.csv")
model_year = model_year.iloc[:,1:]

degree = [1.5, 2, 3]
degree = [3]
for i in range(len(degree)):
    for modelname in models:
        filename2 = glob.glob("../Data/cmip6_alt_m15_NH45_1850_2100/"+modelname+"_alt*hd.nc")[0]# 1850-2100 monthly
        if degree[i]==1.5:
            begin_year = (int(str(model_year[model_year['model']==modelname][str(degree[i])])[-35:-31])-1850)*12
        else:
            begin_year = (int(str(model_year[model_year['model']==modelname][str(degree[i])])[-33:-29])-1850)*12
                
        # Process future alt data
        f = Dataset(filename2)        
        alt2 = f["alt"][begin_year:begin_year + 240, :, :360]

        # Reshape to yearly data and compute mean max values
        alt_yearly2 = np.reshape(alt2, (-1, 12, alt2.shape[1], alt2.shape[2]))  # Reshape to (years, months, lat, lon)
        alt_yearly2 = alt_yearly2.filled(np.nan)
        max_values2 = np.max(alt_yearly2, axis=1)
        mean_max_values2 = np.mean(max_values2, axis=0)

        # Process historical alt data (1982-2014)
        alt = f["alt"][12*(begin-1850):12*(begin-1850)+12*years, :, :360] 
        alt[alt == -9999] = np.nan
        alt_yearly = np.reshape(alt, (-1, 12, alt.shape[1], alt.shape[2]))
        alt_yearly = alt_yearly.filled(np.nan)
        max_values = np.max(alt_yearly, axis=1)
        max_values[max_values == -9999] = np.nan
        mean_max_values = np.mean(max_values, axis=0)

        # Compute delta for ALT
        delta_alt = mean_max_values2 - mean_max_values

        # Process Ts data for 1980-2100
        ts_file3 = glob.glob("/home/wangjx/Data/cmip6_tsl_05_depth_interpolated_198010/" + modelname+"_*.nc")[0]
        if degree[i]==1.5:
            begin_year = (int(str(model_year[model_year['model']==modelname][str(degree[i])])[-35:-31])-1980)*12
        else:
            if int(str(model_year[model_year['model']==modelname][str(degree[i])])[-33:-29])<1000: # Many models do not reach 4°C of warming
                continue
            else:
                begin_year = (int(str(model_year[model_year['model']==modelname][str(degree[i])])[-33:-29])-1980)*12
        f = xr.open_dataset(ts_file3, decode_times=False)
        tsl3 = f['tsl'][begin_year:begin_year + 240, 1, :, :360] # 20 year slice
        mean_tsl_3 = np.mean(tsl3, axis=0)

        # Process historical 0.2m Ts data
        histsl = f['tsl'][12*(begin-1980):12*(begin-1980)+12*years, 1, :, :360]
        his_mean_tsl = np.mean(histsl, axis=0)

        # Compute delta for TSL
        delta_tsl = mean_tsl_3 - his_mean_tsl
        delta_tsl = np.where(delta_tsl < 0.1, np.nan, delta_tsl)

        # Compute the ALT/TSL sensitivity ratio
        alt_tsl = delta_alt / delta_tsl

        # Create xarray Dataset to store the results
        ds = xr.Dataset(
            {'alt_tsl': (['lat', 'lon'], alt_tsl)},
            coords={'lat': f["lat"][:], 'lon': f["lon"][:360]})

        # Add metadata
        ds.attrs['description'] = 'Climate change of ALT sensitivity'
        ds.attrs['model'] = modelname

        # Store the sensitivity of each model at different degrees of warming as nc files
        output_file = f'../Data/ALT_Sensitivity_lat_lon/{degree[i]}degree/{modelname}_alt_climatechange_avg.nc'
        ds.to_netcdf(output_file)