###
    # This code interpolates the horizontal resolution of each model to 0.5° resolution 
    # Using bilinear interpolation method
    # Output: "../Data/cmip6_tsl_05_depth_interpolated"
###

#!/bin/bash
for i in TaiESM1 CESM2-FV2 IPSL-CM6A-LR CESM2-WACCM CESM2 CNRM-CM6-1-HR CNRM-CM6-1 CNRM-ESM2-1 FGOALS-g3 MIROC-ES2L MPI-ESM1-2-HR MPI-ESM1-2-LR NorESM2-MM GFDL-CM4 E3SM-1-1 FGOALS-f3-L GFDL-ESM4 MIROC6 NorESM2-LM

do
cdo remapbil,cmip6_Ts_grid05.txt ../Data/cmip6_tsl_interp_depth_198010/"${i}"_*.nc /home/wangjx/Data/cmip6_tsl_05_depth_interpolated/"${i}"_tsl_monthly_0_3_hd.nc
echo "${i}"
done