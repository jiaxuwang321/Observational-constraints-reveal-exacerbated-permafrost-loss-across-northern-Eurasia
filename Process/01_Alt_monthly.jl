###-------------------------------------------------------------------------------------------------------------------
   # This code is to calculate thickness of active layer (ALT) in each grid.
   # ALT is detected through the maximum depth at which soil temperature crosses the 0°C threshold. 
   # Greenland is masked out.
   # The output is saved as "../Data/cmip6_alt_m15_NH45/modelname_alt_monthly_hd.nc"
###-------------------------------------------------------------------------------------------------------------------

using NCDatasets
using Interpolations
using DataStructures
using Dates
using DataFrames
using Base.Threads

function pfmask(fname)

    f = Dataset(fname)
    msk = f["RCN_reg"][:,:]
    close(f)
 
   # mask out Greenland
    msk[msk.==5] .= -1
 
    msk2 = copy(msk)
 
   # transformation from 180W-180E to 0-360
    msk2[1:360,:] = msk[361:720,:]
    msk2[361:720,:] = msk[1:360,:]
 
    lon = Array(LinRange(0.25,359.75,720))
    lat = Array(LinRange(-89.75,89.75,360))
 
    return (lon, lat, msk2)
 
 end # pfmask

 function monthly_alt(fname,begyr,endyr,msk,msklon,msklat,regrid::Bool)    

    f = Dataset(fname)

    if f["depth"].attrib["units"] != "m"
    #if f["solth"].attrib["units"] != "m"   # In IPSL-CM6A-LR named solth
       println("units of soil depth error")
       exit()
    end

    depth = f["depth"][:]
    #depth = f["solth"][:]
    ndepth = size(depth)[1]

    if depth[end] < 1.945
       println("shallow model (1.945m)")
       return (missing, missing, missing)
    end

    depth = Float32.(depth)

    lon = f["lon"][:]
    lat = f["lat"][:]

    if !regrid    #regard=True
       N45 = findfirst((lat.-45).>0)
       lat = lat[N45:end]
       nlat = size(lat)[1]
       nlon = size(lon)[1]
    else
       N45 = findfirst((msklat.-45).>0)
       msklat = msklat[N45:end]
       msk = msk[:,N45:end]
       nlat = size(msklat)[1]
       nlon = size(msklon)[1]
    end

    nyear = endyr-begyr+1

    monalt_all = Array{Float32}(undef,nlon,nlat,nyear*12)
    fill!(monalt_all, -9999.)

    znode = Array(LinRange(0.01,15.01,1501))
    #znode = Array(LinRange(0.01,3.21,321))

    tsl = f["tsl"]

    time = f["time"][:]

    yymm = Dates.yearmonth.(time)

    for yr in begyr:endyr
       @show yr

       for mn in 1:12
          t = findfirst(yymm.==[(yr,mn)])

          if t==nothing
             @show "non-exist date", yr, mn
             continue 
          end

          if !regrid
             tsl_new = tsl[:,N45:end,:,t]
          else
             tsl_old = tsl[:,:,:,t]
             tsl_new = Array{Union{Float32,Missing}}(undef,nlon,nlat,ndepth)
             msklon2d = Array{Float32}(undef,nlon,nlat)
             msklat2d = Array{Float32}(undef,nlon,nlat)
             msklon2d .= msklon
             msklat2d .= msklat'

             @threads for k in 1:ndepth
                itp = interpolate((lon,lat), tsl_old[:,:,k], Gridded(Linear()))
                etp = extrapolate(itp, (Periodic(),Flat()))
                tsl_new[:,:,k] .= etp.(msklon2d,msklat2d)
             end

             msk[ismissing.(tsl_new[:,:,1])] .= -1

          end

          monalt = Array{Float32}(undef,nlon,nlat)
          fill!(monalt, -9999.)

          for j in 1:nlat
             @threads for i in 1:nlon
                if !ismissing(tsl_new[i,j,1])
                   if regrid & (msk[i,j] < 0)
                      continue
                   end

                   tsoil = convert.(Float32,tsl_new[i,j,:])
                   itp = interpolate((depth,), tsoil, Gridded(Linear()))
                   etp = extrapolate(itp, Flat())
                   ztemp = etp(znode)

                   frozen = ztemp .< 273.15

                   z1 = findlast(frozen.==false) #deepest unfrozen layer
                   z2 = findfirst(frozen.==true) #shallowest frozen layer
                   z3 = findlast(frozen.==true)  #deepest frozen layer

                   if z1 == nothing              # totally frozen
                      monalt[i,j] = 0.
                   elseif z2 == nothing          # totally thawing
                      monalt[i,j] = -9999.
                   elseif (z1+1) == z2           # top thawing & bottom frozen
                      monalt[i,j] = znode[z1]
                   elseif (z1 > z2) & (z1 < z3)  # top frozen & middle thawing & bottom frozen
                      monalt[i,j] = znode[z1]
                   elseif (z1 > z2) | (z1 > z3)  # top frozen & bottom unfrozen
                      monalt[i,j] = -9999.
                   else                          # unclassified
                      @show "unclassified", t, z1, z2, z3
                   end
                end
             end
          end

          t = (yr-begyr)*12+mn
          monalt_all[:,:,t] = monalt
       end
    end

    close(f)

    if !regrid
       return (lon, lat, monalt_all)
    else
       return (msklon, msklat, monalt_all)
    end

end # monthly_alt

function save_alt(fout,lon,lat,alt,begyr,endyr)

    if isfile(fout)
       rm(fout)
    end
 
    nmonth = (endyr-begyr+1)*12
 
    ds = Dataset(fout,"c")
 
    defDim(ds, "lon", size(lon)[1])
    defDim(ds, "lat", size(lat)[1])
    defDim(ds, "time", nmonth)
 
    x = defVar(ds, "lon", Float32, ("lon",), attrib=OrderedDict("units"=>"degrees_E"))
    y = defVar(ds, "lat", Float32, ("lat",), attrib=OrderedDict("units"=>"degrees_N"))
    t = defVar(ds, "time", Float32, ("time",), attrib=OrderedDict("units"=>"months since "*string(begyr)))
    v = defVar(ds, "alt", Float32, ("lon","lat","time"), attrib=OrderedDict("units"=>"m","_FillValue"=>convert(Float32,-9999.)))
 
    x[:] = lon[:]
    y[:] = lat[:]
    t[:] = Array(1:nmonth)
    v[:,:,:] = alt[:,:,:]
 
    close(ds)
 
 end # save_alt

fmask = "/home/wangjx/Data/PCN_Domain_0.5x0.5.nc"
fpath = "/home/jidy/Data/tas-tsl-merge"  # SSP5-8.5
# TaiESM1 is stored in the following path
# fpath = "/home/wangjx/Data/CMIP6_tas_tsl_1850_2100"  

df = DataFrame(model=String[],ftsl=String[])

push!(df, ("CESM2",            "tsl_CESM2_r1i1p1f1_185001-210012.nc"))
push!(df, ("CESM2-WACCM",      "tsl_CESM2-WACCM_r1i1p1f1_185001-210012.nc"))
push!(df, ("CESM2-FV2",        "tsl_CESM2-FV2_r1i2p2f1_185001-210012.nc"))
push!(df, ("CNRM-CM6-1",       "tsl_CNRM-CM6-1_r1i1p1f2_185001-210012.nc"))
push!(df, ("CNRM-ESM2-1",      "tsl_CNRM-ESM2-1_r1i1p1f2_185001-210012.nc"))
push!(df, ("CNRM-CM6-1-HR",    "tsl_CNRM-CM6-1-HR_r1i1p1f2_185001-210012.nc"))
push!(df, ("GFDL-CM4",         "tsl_GFDL-CM4_r1i1p1f1_185001-210012.nc"))
push!(df, ("GFDL-ESM4",        "tsl_GFDL-ESM4_r1i1p1f1_185001-210012.nc"))
push!(df, ("MPI-ESM1-2-LR",    "tsl_MPI-ESM1-2-LR_r1i1p1f1_185001-210012.nc"))
push!(df, ("MPI-ESM1-2-HR",    "tsl_MPI-ESM1-2-HR_r1i1p1f1_185001-210012.nc"))
push!(df, ("NorESM2-LM",       "tsl_NorESM2-LM_r1i1p1f1_185001-210012.nc"))
push!(df, ("NorESM2-MM",       "tsl_NorESM2-MM_r1i1p1f1_185001-210012.nc"))
push!(df, ("MIROC-ES2L",       "tsl_MIROC-ES2L_r1i1p1f2_185001-210012.nc"))
push!(df, ("MIROC6",           "tsl_MIROC6_r1i1p1f1_185001-210012.nc"))
push!(df, ("FGOALS-g3",        "tsl_FGOALS-g3_r1i1p1f1_185001-210012.nc"))
push!(df, ("FGOALS-f3-L",      "tsl_FGOALS-f3-L_r1i1p1f1_185001-210012.nc"))
push!(df, ("E3SM-1-1",         "tsl_E3SM-1-1_r1i1p1f1_185001-210012.nc"))
push!(df, ("CMCC-CM2-SR5",     "tsl_CMCC-CM2-SR5_r1i1p1f1_185001-210012.nc"))
push!(df, ("CMCC-ESM2",        "tsl_CMCC-ESM2_r1i1p1f1_185001-210012.nc"))
push!(df, ("CAS-ESM2-0",       "tsl_CAS-ESM2-0_r1i1p1f1_185001-210012.nc"))
push!(df, ("CanESM5",          "tsl_CanESM5_r1i1p1f1_185001-210012.nc"))
push!(df, ("CanESM5-CanOE",    "tsl_CanESM5-CanOE_r1i1p2f1_185001-210012.nc"))
push!(df, ("KACE-1-0",         "tsl_KACE-1-0_r1i1p1f1_185001-210012.nc"))

#push!(df, ("IPSL-CM6A-LR",      "tsl_IPSL-CM6A-LR_r1i1p1f1_185001-210012.nc"))
#push!(df, ("TaiESM1",           "tsl_TaiESM1_r1i1p1f1_185001-210012.nc"))

nmodels = nrow(df)

begyr = 1850
endyr = 2100

regrid = true

#############

msklon,msklat,msk = pfmask(fmask)

for i in 1:nmodels
   @show df.model[i]

   fname = df.ftsl[i]
   fname = joinpath(fpath,fname)

   lon,lat,alt = monthly_alt(fname,begyr,endyr,msk,msklon,msklat,regrid)

   if !ismissing(alt)
      fout = joinpath("../Data/cmip6_alt_m15_NH45_1850_2100/",df.model[i]*"_alt_monthly_hd.nc")
      save_alt(fout,lon,lat,alt,begyr,endyr)
   end
end