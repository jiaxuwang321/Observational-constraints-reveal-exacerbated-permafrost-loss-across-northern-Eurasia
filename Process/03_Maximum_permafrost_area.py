### 
    # This code is to calculate the probability of permafrost in each grid, across 18 CMIP6 models.
    # 1. Read the ALT file "../Data/cmip6_alt_m15_NH45_1850_2100/", select 1982-2014.
    # 2. Calculate the probability of permafrost during historical period.
    # 3. Output "../Data/18_models_probability_lt_threshold_{depth}.nc".
###

import glob
import numpy as np
from netCDF4 import Dataset
import warnings
warnings.filterwarnings("ignore")

models = ["CESM2", "CESM2-FV2", "CESM2-WACCM", "CNRM-CM6-1-HR", "CNRM-CM6-1", "CNRM-ESM2-1", 
          "E3SM-1-1", "FGOALS-f3-L", "FGOALS-g3", "GFDL-CM4", "GFDL-ESM4", "MIROC6", "MIROC-ES2L",
          "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "NorESM2-LM", "NorESM2-MM", "TaiESM1"]

depths = [320, 240, 160]
for depth in depths:
    # Initialize the counter array for ALT < threshold
    alt_lt_threshold_count = np.zeros((90, 720), dtype=int)

    # Process each model
    for modelname in models:
        print(f"Processing {modelname}")

        filename = glob.glob("../Data/cmip6_alt_m15_NH45_1850_2100/"+modelname+"_alt*hd.nc")[0] # monthly
        with Dataset(filename) as f:
            alt = f["alt"][(1982-1850)*12:(2015-1850)*12]  # 1982-2014 period
            lat = f["lat"][:]
            lon = f["lon"][:]
            
            # Reshape monthly data into yearly and calculate mean annual maximum ALT
            alt = alt.filled(999)  # Replace NaN with large placeholder
            alt_yearly = alt.reshape((-1, 12, *alt.shape[1:]))  # Reshape to (years, months, lat, lon)
            mean_max_values = np.mean(np.max(alt_yearly, axis=1), axis=0)
            
            # Update counter with cells where ALT < threshold (depth/100)
            alt_lt_threshold_count += (mean_max_values <= depth/100).astype(int)

    # Calculate probability of ALT < threshold
    probability_lt_threshold = alt_lt_threshold_count / len(models)

    # Save the result as a NetCDF file
    output_filename = f'../Data/probability_lt_threshold_{depth}.nc'
    with Dataset(output_filename, 'w', format='NETCDF4') as ncfile:
        # Define dimensions
        ncfile.createDimension('lat', len(lat))
        ncfile.createDimension('lon', len(lon))

        # Define variables
        lats = ncfile.createVariable('lat', 'f4', ('lat',))
        lons = ncfile.createVariable('lon', 'f4', ('lon',))
        prob = ncfile.createVariable('probability_lt_threshold', 'f4', ('lat', 'lon',))

        # Assign data to variables
        lats[:] = lat
        lons[:] = lon
        prob[:, :] = probability_lt_threshold

        # Add attributes
        lats.units = 'degrees_north'
        lons.units = 'degrees_east'
        prob.units = 'probability'
        prob.long_name = f'Probability of ALT < {depth}cm in 18 models'
