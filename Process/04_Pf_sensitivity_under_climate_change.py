### 
    # This code is to calculate the permafrost area sensitivity under Climate Change.
    # 1. Represents permafrost area relative to changes in global temperature increases of 1.5°C, 2°C, and 3°C (area/T; unit: %/°C)
    # 2. Equation: 
    #    (20-year average permafrost area at each temperature thresholds - 1982-2014 average permafrost area) /
    #    (air temperature at each thresholds - 1982-2014 average air temperature).
    # 3. Evaluate the sensitivity under two scenarios (SSP5-8.5 & SSP2-4.5).
    # 4. Step:
    #    Calculate the years when reaching the three warming thresholds for CMIP6 models
    #    Calculate the permafrost area when reaching the three warming thresholds for CMIP6 models
    #    Calculate the permafrost sensitivity
    #    Output: "../Data/Tas_data/years_reaching_three_warming_thresholds_ssp245(585).csv"
    #            "../Data/Permafrost_area_timeseries_ssp245(585)/Permafrost_area_"+depth+"_under_three_warming_levels.csv"
    #            "../Data/Permafrost_sensitivity/"
###

import pandas as pd
import numpy as np
from netCDF4 import Dataset
import glob
import warnings
warnings.filterwarnings("ignore")

# model list under SSP2-4.5
models = ["CESM2","CESM2-WACCM","CNRM-CM6-1-HR","CNRM-CM6-1","CNRM-ESM2-1","FGOALS-f3-L",
           "FGOALS-g3","GFDL-CM4","GFDL-ESM4","MIROC6","MIROC-ES2L","MPI-ESM1-2-HR","MPI-ESM1-2-LR",
           "NorESM2-LM","NorESM2-MM","TaiESM1"]

# model list under SSP5-8.5
# models = ["CESM2","CESM2-FV2","CESM2-WACCM","CNRM-CM6-1-HR","CNRM-CM6-1","CNRM-ESM2-1","E3SM-1-1",
#          "FGOALS-f3-L","FGOALS-g3","GFDL-CM4","GFDL-ESM4","MIROC6","MIROC-ES2L","MPI-ESM1-2-HR",
#          "MPI-ESM1-2-LR","NorESM2-LM","NorESM2-MM","TaiESM1"]

# Initialize DataFrames
model_year = pd.DataFrame({'model': models, '1.5': '/', '2': '/', '3': '/'})
model_tas = pd.DataFrame({'model': models, 'pre': np.nan, 'his': np.nan, 
                         '1.5': '/', '2': '/', '3': '/'})

# Define file path templates
pre_path = '/home/wangjx/cmip6_global_tas_1850_1900_monthly/{}_*.nc'
his_path = '/home/wangjx/cmip6_global_tas_1982_2014_monthly/{}_*.nc'
#future_path = '/home/wangjx/cmip6_global_tas_2015_2100_monthly/{}_*.nc' # SSP5-8.5
future_path = '/home/wangjx/cmip6_global_tas_2015_2100_monthly_ssp245/{}_*.nc' # SSP2-4.5

#-------------------------------------------------------------------------------------------------------------
                                        # Several functions
#-------------------------------------------------------------------------------------------------------------
# load nc data
def load_nc_data(filepath):
    with Dataset(filepath) as file:
        return file["tas"][:]

# Calculate annual average (supports multi-dimensional data)
def calculate_annual_average(data):
    num_years = data.shape[0] // 12
    return np.mean(data.reshape(num_years, 12, -1), axis=1).reshape(num_years, *data.shape[1:])

# Function to find threshold years
def find_threshold_years(annual_tas, pre_mean, window_size=20, start_year=1982):
    thresholds_found = {str(t): False for t in [1.5, 2, 3]}
    result_years = {str(t): '/' for t in [1.5, 2, 3]}
    result_tas = {str(t): '/' for t in [1.5, 2, 3]}
    
    for start_idx in range(len(annual_tas) - window_size + 1):
        avg_20_years = np.mean(annual_tas[start_idx:start_idx + window_size])
        difference = avg_20_years - pre_mean
        current_start_year = start_year + start_idx
        
        for threshold in ['1.5', '2', '3']:
            if not thresholds_found[threshold] and difference > float(threshold):
                result_years[threshold] = f"[{current_start_year},{current_start_year + window_size - 1}]"
                result_tas[threshold] = avg_20_years
                thresholds_found[threshold] = True
        
        if all(thresholds_found.values()):
            break
    
    return result_years, result_tas

#-------------------------------------------------------------------------------------------------------------
                    # Determine the years that reaching the three warming thresholds
#-------------------------------------------------------------------------------------------------------------
# Process each model
window_size = 20
all_annual_tas = []  # Store annual temperature series for all models

for model in models:
    # Load data
    tas_pre = load_nc_data(glob.glob(pre_path.format(model))[0])
    tas_his = load_nc_data(glob.glob(his_path.format(model))[0])
    tas_future = load_nc_data(glob.glob(future_path.format(model))[0])
    
    # Calculate averages
    pre_mean = np.mean(tas_pre)
    his_mean = np.mean(tas_his)
    
    # Update base data
    model_tas.loc[model_tas['model'] == model, 'pre'] = pre_mean
    model_tas.loc[model_tas['model'] == model, 'his'] = his_mean
    
    # Calculate annual averages and concatenate
    annual_his = calculate_annual_average(tas_his).flatten()
    annual_future = calculate_annual_average(tas_future).flatten()
    annual_tas = np.concatenate([annual_his, annual_future])
    all_annual_tas.append(annual_tas)
    
    # Find threshold years
    result_years, result_tas = find_threshold_years(annual_tas, pre_mean, window_size)
    
    # Update results
    for threshold in ['1.5', '2', '3']:
        model_year.loc[model_year['model'] == model, threshold] = result_years[threshold]
        model_tas.loc[model_tas['model'] == model, threshold] = result_tas[threshold]

# Calculate multi-model ensemble mean
if all_annual_tas:
    # Ensure all series have the same length
    min_length = min(len(tas) for tas in all_annual_tas)
    trimmed_tas = [tas[:min_length] for tas in all_annual_tas]
    
    # Calculate ensemble mean
    ensemble_mean = np.mean(trimmed_tas, axis=0)
    reference_value = model_tas['pre'].mean()
    
    # Find threshold years for ensemble mean
    result_years, _ = find_threshold_years(ensemble_mean, reference_value, window_size)
    
    # Create ensemble mean row and add to model_year
    ensemble_row = {'model': 'Ensemble_mean'}
    ensemble_row.update(result_years)
    model_year = pd.concat([model_year, pd.DataFrame([ensemble_row])], ignore_index=True)

model_year.to_csv("../Data/Tas_data/years_reaching_three_warming_thresholds_ssp245.csv") # Table S2 # SSP2-4.5
#model_year.to_csv("../Data/Tas_data/years_reaching_three_warming_thresholds_ssp585.csv") # Table S2 # SSP5-8.5

#-------------------------------------------------------------------------------------------------------------
                    # Calculate the permafrost and its sensitivity of each depth threshold
#-------------------------------------------------------------------------------------------------------------
depths = ['320cm','240cm','160cm']
for depth in depths:
    # Initialize DataFrames for historical (1982–2014) and full period (1982–2100)
    begin = 1982
    end = 2014
    nyears = end-begin+1
    min_pfhis = pd.DataFrame({'year': range(begin, end+1)})
    min_pf = pd.DataFrame({'year': range(begin, 2101)})

    for model in models:
        print(model)
        # Load permafrost area data
        filename1 = glob.glob(f"../Data/Permafrost_area_timeseries_ssp585/cmip6_alt{depth}_NH45_EA0_180_pfarea_1850_2100_monthly/{model}_*.nc")[0]
        filename11 = glob.glob(f"../Data/Permafrost_area_timeseries_ssp585/cmip6_alt{depth}_NH45_EA350_360_pfarea_1850_2100_monthly/{model}_*.nc")[0]
        """
        # For SSP5-8.5
        pfarea_climatology = (Dataset(filename1)["pfarea"][(begin-1850)*12:] / 1e12).filled(np.nan)    # Since 1982
        pfarea_climatology1 = (Dataset(filename11)["pfarea"][(begin-1850)*12:] / 1e12).filled(np.nan)  # Since 1982
        """
        # For SSP2-4.5
        base_path2 = "../Data/Permafrost_area_timeseries_ssp245/cmip6_alt" + depth + "_NH45_"
        file_template2 = "{region}_pfarea_{year_range2}_monthly/"+ model + "_pfarea_NH45_{year_range2}monthly.nc"
        filename2 = base_path2 + file_template2.format(region="EA0_180", year_range2="2015_2100")
        filename22 = base_path2 + file_template2.format(region="EA350_360", year_range2="2015_2100")
        f, ff, f2, ff2 =  Dataset(filename1), Dataset(filename11), Dataset(filename2), Dataset(filename22)
        # For SSP2-4.5, merge 1982-2100
        pfarea_climatology = np.concatenate([f["pfarea"][(begin-1850)*12:(end-1850+1)*12] / 10**12, f2["pfarea"][:] / 10**12]) 
        pfarea_climatology1 = np.concatenate([(ff["pfarea"][(begin-1850)*12:(end-1850+1)*12] / 10**12).filled(np.nan), (ff2["pfarea"][:] / 10**12).filled(np.nan)])
        
        # Combine and fill missing values
        pf_data = np.nansum([pfarea_climatology, pfarea_climatology1], axis=0)

        # Calculate min for historical period and add to DataFrame
        min_pfhis[model] = [np.nanmin(pf_data[12*y:12*(y+1)]) for y in range(0, nyears)]
        min_pf[model] = [np.nanmin(pf_data[12*y:12*(y+1)]) for y in range(0, 2101 - begin)]
        
        # Special handling for specific models in 2015 - because the original tsl data is NaN
        if model in ['NorESM2-LM', 'NorESM2-MM', 'TaiESM1']:
            min_pf.loc[min_pf.year == 2015, model] = np.nan

    # Calculate permafrost area average for each warming threshold
    model_perma = pd.DataFrame({'model': models, 'his': min_pfhis.mean().values[1:], '1.5': '/', '2': '/', '3': '/'})

    for idx, model in enumerate(models):
        for threshold in ['1.5', '2', '3']:
            time_period = model_year.loc[model_year['model'] == model, threshold].values[0]
            if time_period != '/':
                start, end = map(int, (time_period[1:5], time_period[6:10]))
                pf2_mean = min_pf.loc[(min_pf.year >= start) & (min_pf.year <= end), model].mean()
                model_perma.at[idx, threshold] = pf2_mean

    # Add ensemble mean values and standard deviation across models
    mean_data = {'model': 'ensemble_mean'}
    std_data = {'model': 'ensemble_std'}
    # for threshold in ['his','1.5', '2', '3']: # For SSP5-8.5
    for threshold in ['his','1.5', '2']: # For SSP2-4.5, lots of models do not reach 3°C warming in 21st
        mean_data[threshold] = model_perma[threshold].mean()
        std_data[threshold] = model_perma[threshold].iloc[:-1].std()

    model_perma = pd.concat([model_perma, pd.DataFrame([mean_data]), pd.DataFrame([std_data])], ignore_index=True)

    # model_perma.to_csv("../Data/Permafrost_area_timeseries_ssp585/Permafrost_area_"+depth+"_under_three_warming_levels_ssp585.csv") # Table S3 # SSP5-8.5
    model_perma.to_csv("../Data/Permafrost_area_timeseries_ssp245/Permafrost_area_"+depth+"_under_three_warming_levels_ssp245.csv") # Table S3 # SSP2-4.5

    # Initialize DataFrame for permafrost sensitivity under climate change
    model_climsen = pd.DataFrame({'model': models, '1.5': '/', '2': '/', '3': '/'})

    # Calculate sensitivity for each model and threshold
    for idx, modelname in enumerate(model_climsen['model']):
        for threshold in ['1.5', '2', '3']:
            if model_perma.loc[idx, threshold]=='/':
                model_climsen.loc[idx, threshold] = np.nan
            else:
                
                # Calculate percentage decrease in permafrost area and temperature increase for the threshold
                perma_decrease = (float(model_perma.loc[idx, threshold]) - model_perma.loc[idx, 'his']) / model_perma.loc[idx, 'his'] * 100
                tas_decrease = float(model_tas.loc[idx, threshold]) - model_tas.loc[idx, 'his']
                
                # Store the climate sensitivity
                model_climsen.loc[idx, threshold] = perma_decrease / tas_decrease

    # model_climsen.to_csv("../Data/Permafrost_sensitivity/"+depth+"_permafrost_sensitivity_under_climate_change_ssp585.csv")  # SSP5-8.5
    model_climsen.to_csv("../Data/Permafrost_sensitivity/"+depth+"_permafrost_sensitivity_under_climate_change_ssp245.csv")  # SSP2-4.5