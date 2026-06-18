###
    # This code is the moisture-driven model, which employed the soil heat transfer model adapted from the Common Land Model (Dai et al. 2003, 2004, 2019a). 
    # 1. Soil physical parameters for the upper 10 layers were derived from Dai et al. (2019b) based on the Global Soil Dataset for Earth System Modeling (Shangguan et al., 2014).
    # 2. This model is driven by surface heat flux and simulate soil temperature evolution via the one-dimensional heat conduction equation, where temperature changes depend on vertical heat flux divergence.
    # 3. In this model, only freeze–thaw processes are considered; precipitation, runoff, and other hydrological processes are excluded so that the total soil water content remains constant.
    # 4. Across different simulations, the only factor varying is the climatological soil water content (liquid+ice) prescribed from different CMIP6 models,so-called mositure-driven model.
    # 5. The supercooled water switch is in the melting step.
    # 6. Output: (1)offsup:"../Data/Moisture_driven/output/{modelname}_moisture_driven_soil_Ts_offsup.nc" (2)onsup:"../Data/Moisture_driven/output/{modelname}_moisture_driven_soil_Ts_onsup.nc". 
###

import xarray as xr
import numpy as np
import pandas as pd
import glob
from datetime import datetime, timedelta
from joblib import Parallel, delayed
import warnings
warnings.filterwarnings('ignore')

def soil_hcap_cond(vf_gravels_s,vf_om_s,vf_sand_s,vf_pores_s,
                   csol,kdry,ksat_u,ksat_f,
                   BA_alpha,BA_beta,
                   temperature,
                   vf_water,vf_ice):
    
    #-----------------------------------------------------------------------
    # DESCRIPTION:
    # Calculate bulk soil heat capacity and soil thermal conductivity
    # -----------------------------------------------------------------------------------------

    # -----------------------------------------------------------------------------------------
    # The heat capacity and thermal conductivity [J(m3 K)]
    # -----------------------------------------------------------------------------------------
    c_water = 4.188e6   # J/(m3 K) = 4188[J/(kg K)]*1000(kg/m3)
    c_ice = 1.94153e6   # J/(m3 K) = 2117.27[J/(kg K)]*917(kg/m3)

    hcap = csol + vf_water*c_water + vf_ice*c_ice 
    
    # -----------------------------------------------------------------------------------------
    # Setting
    # -----------------------------------------------------------------------------------------
    tfrz = 273.16    #freezing temperature [K]
    #tfrz = 271.16    #freezing temperature [K]   # GFDL model

    sr = (vf_water+vf_ice)/vf_pores_s  # wetness or degree of saturation
    sr = min(1.0, sr)
        
    if sr >= 1.0e-10:
     
    # =========================================================================================
    # [4] Balland V. and P. A. Arp, 2005: Modeling soil thermal conductivities over a wide range of conditions. J. Environ. Eng. Sci. 4: 549-558.
    # Be careful in specifying all k affecting fractions as VOLUME FRACTION, whether these fractions are part of the bulk volume, the pore space, or the solid space.
    # =========================================================================================
        if (temperature > tfrz) : # Unfrozen soil
            # Kersten number or normalized thermal conductivity
            ke = sr**(0.5*(1.0+vf_om_s-BA_alpha*vf_sand_s-vf_gravels_s)) * ((1.0/(1.0+np.exp(-BA_beta*sr)))**3-((1.0-sr)/2.0)**3)**(1.0-6)
        else:                     # Fozen or partially frozen soils
            ke = sr**(1.0+vf_om_s)
    
    else:
        ke = 0.0
    
    ke = max(ke, 0.0)
    ke = min(ke, 1.0)
    
    if (temperature > tfrz):      # Unfrozen soil
        thk = (ksat_u-kdry)*ke + kdry
    else:                         # Frozen or partially frozen soils
        thk = (ksat_f-kdry)*ke + kdry
           
    return hcap,thk               # J/(m3 K), W/(m K)

def tridia(n: int, a: np.ndarray, b: np.ndarray, c: np.ndarray, r: np.ndarray) -> np.ndarray:
    """
    Solve tridiagonal linear system using Thomas algorithm (Tridiagonal Matrix Algorithm)
    
    System form:
    b[0]  c[0]   0     ...    0     | u[0]     r[0]
    a[1]  b[1]  c[1]   ...    0     | u[1]     r[1]
     0    a[2]  b[2]   ...    0     | u[2]  =  r[2]
    ...   ...   ...    ...   ...    | ...      ...
     0     0     0   a[n-1] b[n-1]  | u[n-1]   r[n-1]
    
    Parameters
    ----------
    n : int
        Matrix dimension (number of diagonal elements)
    a : np.ndarray
        Sub-diagonal elements (a[0] is unused; valid elements from a[1] to a[n-1])
    b : np.ndarray  
        Main diagonal elements (b[0] to b[n-1])
    c : np.ndarray
        Super-diagonal elements (c[0] to c[n-2]; c[n-1] is unused)
    r : np.ndarray
        Right-hand side vector
    
    Returns
    -------
    u : np.ndarray
        Solution vector
    """
    
    # Initialize arrays
    u = np.zeros(n, dtype=np.float64)
    gam = np.zeros(n, dtype=np.float64)
    
    # Forward elimination
    bet = b[0]
    u[0] = r[0] / bet
    
    for j in range(1, n):
        gam[j] = c[j-1] / bet
        bet = b[j] - a[j] * gam[j]
        u[j] = (r[j] - a[j] * u[j-1]) / bet
    
    # Back substitution
    for j in range(n-2, -1, -1):  # from n-2 down to 0
        u[j] = u[j] - gam[j+1] * u[j+1]
    
    return u

denice     = 917       # density of ice [kg/m3]
denh2o     = 1000      # density of liquid water [kg/m3]
cpliq      = 4188      # specific heat of water [J/kg-K]
cpice      = 2117.27   # specific heat of ice [J/kg-K]
tkair      = 0.023     # thermal conductivity of air [W/m/k]
tkice      = 2.290     # thermal conductivity of ice [W/m/k]
tk_bedrock = 3         # thermal conductivity of bedrock [W/m/K]
c_bedrock  = 2e6       # heat capacity of bedrock [J/m3/K]
scv        = 0         # snow cover, water equivalent [mm, kg/m2]
snowdp     = 0         # snow depth[m]
patchtype  = 0         # land patch type (0=soil,1=urban or built-up,2=wetland,3=land ice, 4=deep lake, 5=shallow lake)
deltim     = 3600      # seconds in a time step [second] = 1h
capr       = 0.34      # tuning factor to turn first layer T into surface T, read from“Case_restart_const_lc2005.nc” capr=0.34

# soil deoth
nl_soil = 15           # upper bound of array
x_array = np.array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])
z_soisno = 0.025*(np.exp(0.5*(x_array-0.5))-1)  # node depth [m]
print('node depth:',[f'{x:.5f}' for x in z_soisno])

dz_soisno = np.zeros(nl_soil)    # soil layer thickiness [m]
dz_soisno[0] = 0.5*(z_soisno[1]+z_soisno[0])
for d in range(1,nl_soil-1):
    dz_soisno[d] = 0.5*(z_soisno[d+1]-z_soisno[d-1])
dz_soisno[nl_soil-1] = z_soisno[nl_soil-1] - z_soisno[nl_soil-2]

zi_soisno = np.zeros(nl_soil+1)  # interface depth [m]
zi_soisno[0]=0.
for d in range(1,nl_soil+1):
    zi_soisno[d]=zi_soisno[d-1]+dz_soisno[d-1]
print('interface depth',[f'{x:.5f}' for x in zi_soisno])

# Constants
days = 365 * 101
seconds_per_day = 86400
seconds_per_year = 86400 * 365.25

# Generate time array (unit: seconds)
t_total = days * seconds_per_day
dt = deltim  # Resolution in hours
t = np.arange(0, t_total, dt)
# Create datetime index
start_date = datetime(1700, 1, 1, 0, 0, 0)
t_dates = [start_date + timedelta(seconds=float(ts)) for ts in t]

def generate_hs_timeseries_corrected(hs_annual_amplitude):
    """
    Generate a realistic soil heat flux (hs) time series.
    Long-term mean is zero, satisfying energy conservation.
    """
    
    # 1. Annual cycle: slightly more heat absorption in summer, more release in winter, zero annual mean
    noise_level = 3
    # Annual cycle signal: positive in summer, negative in winter
    # Using sin(2*pi*t/year) to ensure positive-negative symmetry
    annual_signal = hs_annual_amplitude * np.sin(2 * np.pi * t / seconds_per_year - np.pi/2)
    
    # 2. Diurnal cycle: positive during daytime, negative at night, zero daily mean
    # Diurnal amplitude varies seasonally: large in summer, small in winter
    diurnal_amplitude_summer = 60  # W/m² (summer diurnal amplitude)
    diurnal_amplitude_winter = 20  # W/m² (winter diurnal amplitude)
    
    # Envelope of diurnal amplitude with annual variation
    diurnal_amplitude = (diurnal_amplitude_summer + diurnal_amplitude_winter)/2 + \
                       ((diurnal_amplitude_summer - diurnal_amplitude_winter)/2) * \
                       np.sin(2 * np.pi * t / seconds_per_year - np.pi/2)
    
    # Diurnal cycle signal (peak at 2 PM, zero daily mean)
    phi_daily = 0.58 * 2 * np.pi  # 14/24 * 2pi
    diurnal_signal = diurnal_amplitude * np.sin(2 * np.pi * t / seconds_per_day + phi_daily)
    
    # 3. Combine signals - ensure long-term mean is zero
    hs_smooth = annual_signal + diurnal_signal
    
    # Verify long-term mean
    long_term_mean = np.mean(hs_smooth)
    #print(f"Long-term mean of smooth signal: {long_term_mean:.6f} W/m² (should be close to 0)")
    
    # 4. Add random noise (zero mean)
    noise = np.random.normal(0, noise_level, len(t))
    hs = hs_smooth + noise
    
    return hs

def meltf(patchtype,nl_soil,deltim, csol,
                     fact,brr,hs,
                     hst1,
                     t_soisno_bef,t_soisno,wliq_soisno,wice_soisno,
                     scv,snowdp,porsl,psi0,
                     bsw, # Campbell_SOIL_MODEL
                     dz,z_soisno,zi_soisno):

#-----------------------------------------------------------------------
# DESCRIPTION:
# calculation of the phase change within snow and soil layers:
# (1) check the conditions which the phase change may take place,
#     i.e., the layer temperature is great than the freezing point
#     and the ice mass is not equal to zero (i.e., melting),
#     or layer temperature is less than the freezing point
#     and the liquid water mass is not equal to zero (i.e., freezing);
# (2) assess the rate of phase change from the energy excess (or deficit)
#     after setting the layer temperature to freezing point;
# (3) re-adjust the ice and liquid mass, and the layer temperature
#-----------------------------------------------------------------------

    hfus    = 0.3336e6  # latent heat of fusion for ice [J/kg]
    tfrz    = 273.16    # freezing temperature [K]
    #tfrz    = 271.16    # freezing temperature [K] for GFDL model
    grav    = 9.80616   # gravity constant [m/s2]
    c_water = 4.188e6   # J/(m3 K) = 4188[J/(kg K)]*1000(kg/m3)
    c_ice   = 1.94153e6 # J/(m3 K) = 2117.27[J/(kg K)]*917(kg/m3)
    denice  = 917       # density of ice [kg/m3]
    denh2o  = 1000      # density of liquid water [kg/m3]
    capr    = 0.34      # tuning factor to turn first layer T into surface T, read from “Case_restart_const_lc2005.nc” capr=0.34

    DEF_USE_SUPERCOOL_WATER = False       #！！！！！！the swich！！！！！！#
    imelt = np.zeros(nl_soil)
    hm = np.zeros(nl_soil)
    xm = np.zeros(nl_soil)
    wice0 = np.zeros(nl_soil)
    wliq0 = np.zeros(nl_soil)
    wmass0 = np.zeros(nl_soil)
    supercool = np.zeros(nl_soil)

    vf_water = np.zeros(nl_soil)
    vf_ice = np.zeros(nl_soil)
    cv = np.zeros(nl_soil)
    hcap = np.zeros(nl_soil)
    
    for j in range(nl_soil):
        imelt[j]  = 0  # imelt=1:metling; imelt=2:freezing
        hm[j]     = 0  # energy residual [W/m2]
        xm[j]     = 0  # metling or freezing within a time step [kg/m2]
        wice0[j]  = wice_soisno[j]
        wliq0[j]  = wliq_soisno[j]
        wmass0[j] = wice_soisno[j] + wliq_soisno[j]
        
    # supercooling water
    if (DEF_USE_SUPERCOOL_WATER):
        for j in range(nl_soil):
            supercool[j] = 0  # the maximum liquid water when the soil temperature is below the freezing point [mm3/mm3]
            if (t_soisno[j] < tfrz and patchtype <=2):
                smp = hfus * (t_soisno[j]-tfrz)/(grav*t_soisno[j]) * 1000     # mm
                
                if (porsl[j] > 0):
                    supercool[j] = porsl[j]*(smp/psi0[j])**(-1/bsw[j])
                else:
                    supercool[j] = 0
                supercool[j] = supercool[j]*dz[j]*1000           # mm
    
    for j in range(nl_soil):
        
        # Melting identification
        # if ice exists above melt point, melt some to liquid.
        if (wice_soisno[j] > 0 and t_soisno[j] > tfrz):
            imelt[j] = 1  # flag for melting or freezing [-]
            t_soisno[j] = tfrz
        
        # Freezing identification
        if (DEF_USE_SUPERCOOL_WATER):
            if ((wliq_soisno[j] > supercool[j]) and (t_soisno[j] < tfrz)):
                imelt[j] = 2
                t_soisno[j] = tfrz

        else:
            if (wliq_soisno[j] > 0 and t_soisno[j] < tfrz):
                imelt[j] = 2
                t_soisno[j] = tfrz

        # Calculate the energy surplus and loss for melting and freezing
        if (imelt[j] > 0):
            tinc = t_soisno[j]-t_soisno_bef[j]
            if j > 0:
                hm[j] = brr[j] - tinc/fact[j]        # hm: energy surplus or deficit after temperature adjustment (W/m2)
            else:
                hm[j] = hst1 + brr[j] - tinc/fact[j]
             
        # This error was checked carefully, it results from the the computed error of "Tridiagonal-Matrix" in subroutine "thermal".
        if (imelt[j] == 1 and hm[j] < 0.):
            hm[j] = 0
            imelt[j] = 0
        if (imelt[j] == 2 and hm[j] > 0.):
            hm[j] = 0
            imelt[j] = 0

        # The rate of melting and freezing
        if imelt[j] > 0 and abs(hm[j]) > 0:
            xm[j] = hm[j]*deltim/hfus                 # xm: Mass adjustment due to phase change (kg/m2)
            heatr = 0  # energy residual or loss after melting or freezing
            if (xm[j] > 0):
                wice_soisno[j] = max(0, wice0[j]-xm[j])                        
                wliq_soisno[j] = max(0, wmass0[j]-wice_soisno[j])
                
            else:
                if (DEF_USE_SUPERCOOL_WATER):
                    if (wmass0[j] < supercool[j]):
                        wice_soisno[j] = 0
                    else:
                        wice_soisno[j] = min(wmass0[j]-supercool[j], wice0[j]-xm[j])
                    wliq_soisno[j] = max(0,wmass0[j]-wice_soisno[j])

                else:
                    wice_soisno[j] = min(wmass0[j], wice0[j]-xm[j])
                    wliq_soisno[j] = max(0,wmass0[j]-wice_soisno[j])
            
            # If the moisture adjustment process is insufficient to consume (compensate for) the entire energy surplus (deficit), recalculate the energy residual [W/m2]
            heatr = hm[j] - hfus*(wice0[j]-wice_soisno[j])/deltim 
            if (abs(heatr) > 0):
                
                # Since water and ice content have changed, the heat capacity cv needs to be updated, the updated fact will be used for subsequent calculations
                if j<10:
                    vf_water[j] = wliq_soisno[j]/(dz_soisno[j]*denh2o)  # Volumetric fraction of water in pore space
                    vf_ice[j] = wice_soisno[j]/(dz_soisno[j]*denice) 
                    hcap[j] = csol[j] + vf_water[j]*c_water + vf_ice[j]*c_ice
                    cv[j] = hcap[j]*dz[j]     # heat capacity of soil [J/(m2 K)]
                else:
                    hcap[j] = c_bedrock 
                    cv[j] = hcap[j]*dz[j]     # heat capacity of bedrock [J/(m2 K)]
                if j > 0:
                    fact[j] = deltim/cv[j]
                    t_soisno[j] = t_soisno[j] + fact[j]*heatr

                else:
                    fact[j]= deltim / cv[j] *dz[j]/ (0.5*(z_soisno[j]-zi_soisno[j]+capr*(z_soisno[j+1]-zi_soisno[j])))
                    t_soisno[j] = t_soisno[j] + fact[j]*heatr
                
                if (DEF_USE_SUPERCOOL_WATER):
                    if (j < 0 or patchtype == 3): #snow
                        if (wliq_soisno[j]*wice_soisno[j] > 0):
                            t_soisno[j] = tfrz
                else:
                    if (wliq_soisno[j]*wice_soisno[j] > 0.):
                        t_soisno[j] = tfrz

    return t_soisno[:],wliq_soisno[:],wice_soisno[:]

def monthly_mean(arr):
    df = pd.DataFrame(arr, index=t_dates)
    return df.resample("M").mean().values

def make_blocks(nlat, nlon, block=8):

    blocks = []

    for i in range(0, nlat, block):
        for j in range(0, nlon, block):

            blocks.append((
                i, min(i+block, nlat),
                j, min(j+block, nlon)
            ))

    return blocks

num_times = 24 * 365 * 101  # In hours, for a total of 101 years

models = ["CESM2", "CESM2-FV2", "CESM2-WACCM",
          "CNRM-CM6-1", "CNRM-ESM2-1",
          #"GFDL-ESM4",
          "NorESM2-LM", "NorESM2-MM","TaiESM1"
          ]

for model in models:

    # Soil water content
    mrsll_path = glob.glob(
        f"../Data/Moisture_driven/Soil_water/mrsll_{model}_historical_1982-2014_CoLM_vertical_pf_region.nc"
    )[0]

    mrsfl_path = glob.glob(
        f"../Data/Moisture_driven/Soil_water/mrsfl_{model}_historical_1982-2014_CoLM_vertical_pf_region.nc"
    )[0]

    # Soil temperature
    tsl_path = glob.glob(
        f"../Data/Moisture_driven/Soil_temperature/multi_model_mean_tsl_1982-2014_CoLM_vertical_pf_region.nc"
    )[0]

    # Soil heat flux at the upper boundary
    hs_path = glob.glob(f"../Data/Moisture_driven/hs/multi_model_mean_hs.nc")[0]

    # -----------------------------------------------------
    # Read data
    # -----------------------------------------------------
    ds_mrsll = xr.open_dataset(mrsll_path).sel(lon=slice(0, 180))
    ds_mrsfl = xr.open_dataset(mrsfl_path).sel(lon=slice(0, 180))
    ds_tsl = xr.open_dataset(tsl_path)
    ds_hs = xr.open_dataset(hs_path)

    # Historical mean
    mrsll = ds_mrsll["mrsll"]   # (10, lat, lon)
    mrsfl = ds_mrsfl["mrsfl"]
    tsl = ds_tsl["tsl"]         # (15, lat, lon)
    hs_array = ds_hs["hs"]      # (lat, lon)

    # Lat, lon
    lat = ds_mrsfl['lat']
    lon = ds_mrsfl['lon']
    nlat = len(lat)
    nlon = len(lon)

    # Soil parameter
    ds_param = {}
    # adjusted
    for v in ['theta_s','psi_s','tkdry','tksatf','tksatu']:
        ds = xr.open_dataset(f"../Data/Moisture_driven/soil_parameter/{v}_adjusted_05.nc")
        ds_param[v] = ds[v].transpose("layer","lat","lon").values

    # remapcon
    for v in ['csol','lambda','vf_sand_s','vf_om_s','vf_gravels_s']:
        ds = xr.open_dataset(f"../Data/Moisture_driven/soil_parameter/{v}_remapcon_05.nc")
        ds_param[v] = ds[v].transpose("layer","lat","lon").values

    porsl_10      = ds_param['theta_s']       # porosity [-]
    psi_s         = ds_param['psi_s']         # soil water suction, negative potential [mm]
    tkdry         = ds_param['tkdry']
    tksatu        = ds_param['tksatu']
    tksatf        = ds_param['tksatf']
    csol          = ds_param['csol']
    soil_lambda_l = ds_param['lambda']
    vf_om         = ds_param['vf_om_s']
    vf_sand       = ds_param['vf_sand_s']
    vf_gravels    = ds_param['vf_gravels_s']
    BA_alpha      = np.full(nl_soil, 0.24)
    BA_beta       = np.full(nl_soil, 18.1)
    bsw           = 1./soil_lambda_l          # clapp and hornbereger "b" parameter [-]

    def run_one_pixel(m, n):

        if np.all(np.isnan(mrsll[:, m, n])):
            return None

        hs = generate_hs_timeseries_corrected(hs_array[m,n].data)

        # ========== State variable ==========
        t_soil = np.empty((num_times+1, nl_soil))
        wliq   = np.empty_like(t_soil)
        wice   = np.empty_like(t_soil)
        cv = np.empty((num_times,nl_soil))
        tk = np.empty((num_times,nl_soil))
        t_soil_bef = np.empty((num_times,nl_soil))
        
        # Initial values 
        t_soil[0,:] = tsl[:, m, n]      # soil temperature [K]

        wliq[:] = 0.0
        wice[:] = 0.0

        for i in range(10):  
            wliq[0,i] = mrsll[i, m, n]  # liquid water [kg/m2]
            wice[0,i] = mrsfl[i, m, n]  # ice lens [kg/m2]
        try:
            la, = np.where(wliq[0,:10]+wice[0,:10] < 0.2)
            bedrock_layer = np.min([min(la),10])
        except:
            bedrock_layer = 10

        # ========== Assign parameters ==========
        # The bottom five layers are bedrock: set porosity to 0
        porsl      = np.zeros(15)
        porsl[:bedrock_layer] = porsl_10[:bedrock_layer,m,n]
        psi0  = psi_s[:, m, n] * 10
        bsw0  = bsw[:, m, n]
        
        vf_om0   = vf_om[:, m, n]
        vf_sand0 = vf_sand[:, m, n]
        vf_gra0  = vf_gravels[:, m, n]

        tksatu0 = tksatu[:, m, n]
        tksatf0 = tksatf[:, m, n]
        tkdry0  = tkdry[:, m, n]
        csol0   = csol[:, m, n]

        # ===== Temporary variable =====
        hcap = np.zeros(nl_soil)
        thk  = np.zeros(nl_soil)

        fact = np.zeros(nl_soil)
        fn   = np.zeros(nl_soil)
        fn1  = np.zeros(nl_soil)

        at = np.zeros(nl_soil)
        bt = np.zeros(nl_soil)
        ct = np.zeros(nl_soil)
        rt = np.zeros(nl_soil)

        brr = np.zeros(nl_soil)
        
        # === Time integration ===
        for t in range(num_times-1):

            # ---------- Upper soil layers ----------
            for i in range(bedrock_layer):
        
                vf_water = wliq[t,i]/(dz_soisno[i]*denh2o)           # Volumetric fraction of water in pore space
                vf_ice   = wice[t,i]/(dz_soisno[i]*denice)           # Volumetric fraction of ice in pore space

                hcap[i], thk[i] = soil_hcap_cond(
                    vf_gra0[i], vf_om0[i], vf_sand0[i], porsl[i],
                    csol0[i], tkdry0[i], tksatu0[i], tksatf0[i],
                    BA_alpha[i], BA_beta[i],
                    t_soil[t,i],
                    vf_water, vf_ice
                )

                cv[t,i] = hcap[i] * dz_soisno[i]     # heat capacity of soil [J/(m2 K)]

            # ---------- Bedrock ----------        
            for i in range(bedrock_layer,nl_soil):
                thk[i]  = tk_bedrock                 # thermal conductivity of bedrock [W/(m K)]
                hcap[i] = c_bedrock
                cv[t,i] = c_bedrock*dz_soisno[i]     # heat capacity of bedrock [J/(m2 K)]
            
            # Thermal conductivity at the layer interface [W/(m K)]
            for i in range(nl_soil-1):
            
                tk[t,i] = thk[i]*thk[i+1]*(z_soisno[i+1]-z_soisno[i])/(thk[i]*(z_soisno[i+1]-zi_soisno[i+1])+thk[i+1]*(zi_soisno[i+1]-z_soisno[i]))
            
            tk[t,nl_soil-1] = 0.
        
            t_soil_bef[t,:] = t_soil[t,:]
            
            # The original script directly substituted the fact calculated here into the meltf function to compute the temperature adjustment magnitude (heatr),
            # without promptly updating cv corresponding to the changes in water and ice content.
            # This script incorporates the correction within the meltf function.
            j       = 0
            fact[j] = deltim / cv[t,j] *dz_soisno[j]/ (0.5*(z_soisno[j]-zi_soisno[j]+capr*(z_soisno[j+1]-zi_soisno[j])))
            
            for j in range(1,nl_soil):
                fact[j] = deltim/cv[t,j]
        
            for j in range(nl_soil-1):
                fn[j] = tk[t,j]*(t_soil[t,j+1]-t_soil[t,j])/(z_soisno[j+1]-z_soisno[j])
            
            fn[nl_soil-1] = 0.
            
            # cnfac: Crank Nicholson factor between 0 and 1, read from “Case26_restart_const_lc2005.nc” cnfac=0.5
            cnfac=0.5
        
            # set up vector r and vectors a, b, c that define tridiagonal matrix
            j = 0
            dzp   = z_soisno[j+1]-z_soisno[j]
            at[j] = 0.
            bt[j] = 1+(1.-cnfac)*fact[j]*tk[t,j]/dzp
            ct[j] =  -(1.-cnfac)*fact[j]*tk[t,j]/dzp
            rt[j] = t_soil[t,j] + fact[j]*( hs[t+1] + cnfac*fn[j] )
            
            for j in range(1, nl_soil - 1):
                dzm   = (z_soisno[j]-z_soisno[j-1])
                dzp   = (z_soisno[j+1]-z_soisno[j])
                at[j] =   - (1.-cnfac)*fact[j]* tk[t,j-1]/dzm
                bt[j] = 1.+ (1.-cnfac)*fact[j]*(tk[t,j]/dzp + tk[t,j-1]/dzm)
                ct[j] =   - (1.-cnfac)*fact[j]* tk[t,j]/dzp
                rt[j] = t_soil[t,j] + cnfac*fact[j]*( fn[j] - fn[j-1] )
            
            j     =  nl_soil -1
            dzm   = (z_soisno[j]-z_soisno[j-1])
            at[j] =   - (1.-cnfac)*fact[j]*tk[t,j-1]/dzm
            bt[j] = 1.+ (1.-cnfac)*fact[j]*tk[t,j-1]/dzm
            ct[j] = 0.
            rt[j] = t_soil[t,j] - cnfac*fact[j]*fn[j-1]
            
            # solve for t_soisno[t+1,:]
            i = np.size(at)
            t_soil[t+1,:] =  tridia (i ,at ,bt ,ct ,rt)

            #=======================================================================
            # melting or freezing
            #=======================================================================
            for j in range(nl_soil - 1):
                fn1[j] = tk[t,j]*(t_soil[t+1,j+1]-t_soil[t+1,j])/(z_soisno[j+1]-z_soisno[j])
            fn1[nl_soil-1] = 0.
            
            brr[0] = cnfac*fn[0] + (1.-cnfac)*fn1[0]
            for j in range(1,nl_soil):
                brr[j] = cnfac*(fn[j]-fn[j-1]) + (1.-cnfac)*(fn1[j]-fn1[j-1])
        
            t_soil[t+1,:],wliq[t+1,:],wice[t+1,:] =  meltf (patchtype,nl_soil,deltim,csol0,
                    fact,brr,hs[t],
                    hs[t+1],
                    t_soil_bef[t,:].copy(),t_soil[t+1,:].copy(),wliq[t,:].copy(),wice[t,:].copy(),
                    scv,snowdp,porsl,psi0,
                    bsw0, # Campbell_SOIL_MODEL
                    dz_soisno,z_soisno,zi_soisno )
        
        # calculate monthly mean
        Tmon   = monthly_mean(t_soil[:-1])
        Wliq   = monthly_mean(wliq[:-1])
        Wice   = monthly_mean(wice[:-1])
        CVmon  = monthly_mean(cv)
        TKmon  = monthly_mean(tk)

        return m, n, Tmon, Wliq, Wice, CVmon, TKmon

    def run_block(i1,i2,j1,j2):
        out=[]
        for m in range(i1,i2):
            for n in range(j1,j2):
                r = run_one_pixel(m,n)
                if r is not None:
                    out.append(r)
        return out

    ncore = 60
    block = 8

    blocks = make_blocks(nlat,nlon,block)

    results = Parallel(n_jobs=ncore, backend="loky", verbose=10)(
        delayed(run_block)(*b) for b in blocks
    )

    nmon = 12*101
    Tmon  = np.full((nmon,nl_soil,nlat,nlon),np.nan)
    WliqMon = np.full_like(Tmon,np.nan)
    WiceMon = np.full_like(Tmon,np.nan)
    CVmon = np.full_like(Tmon,np.nan)
    TKmon = np.full_like(Tmon,np.nan)

    for block in results:

        for m,n,T,Wl,Wi,Cv,Tk in block:

            Tmon[:,:,m,n]   = T
            WliqMon[:,:,m,n]= Wl
            WiceMon[:,:,m,n]= Wi
            CVmon[:,:,m,n]  = Cv
            TKmon[:,:,m,n]  = Tk

    ds_out = xr.Dataset(
        data_vars = {
            "tsoil": (("time","layer","lat","lon"), Tmon),
            "wliq":  (("time","layer","lat","lon"), WliqMon),
            "wice":  (("time","layer","lat","lon"), WiceMon),
            "cv":    (("time","layer","lat","lon"), CVmon),
            "tk":    (("time","layer","lat","lon"), TKmon),
        },
        coords = {
            "time": np.arange(nmon),
            "layer": np.arange(nl_soil),
            "lat": lat,
            "lon": lon
        }
    )

    ds_out.to_netcdf("../Data/Moisture_driven/output/"+model+"_moisture_driven_soil_Ts_offsup.nc", format="NETCDF4")
    #ds_out.to_netcdf("../Data/Moisture_driven/output/"+model+"_moisture_driven_soil_Ts_onsup.nc", format="NETCDF4")
