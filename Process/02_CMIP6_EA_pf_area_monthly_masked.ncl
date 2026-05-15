;;;
    ; This code is to calculate permafrost area based on ALT data (calculated by "../Process/01*.py" script)
    ; Uses ALT less than or equal to three thresholds (3.2, 2.4, 1.6) to determine permafrost existence
    ; Data is stored in the "../Data/permafrost_area_timeseries/" folder
;;;

begin

; Model list for SSP5-8.5
listi = (/"CAS-ESM2-0","CESM2-WACCM","CESM2","CESM2-FV2","CMCC-CM2-SR5","CMCC-ESM2","CNRM-CM6-1-HR","CNRM-CM6-1","CNRM-ESM2-1","FGOALS-g3","MIROC-ES2L","MPI-ESM1-2-HR","MPI-ESM1-2-LR","NorESM2-MM","CanESM5-CanOE","CanESM5","GFDL-CM4","E3SM-1-1","FGOALS-f3-L","GFDL-ESM4","IPSL-CM6A-LR","KACE-1-0","MIROC6","NorESM2-LM","TaiESM1"/)
; Model list for SSP2-4.5
;listi = (/"CESM2-WACCM","CESM2","CNRM-CM6-1-HR","CNRM-CM6-1","CNRM-ESM2-1","FGOALS-g3","MIROC-ES2L","MPI-ESM1-2-HR","MPI-ESM1-2-LR","NorESM2-MM","GFDL-CM4","FGOALS-f3-L","GFDL-ESM4","IPSL-CM6A-LR","MIROC6","NorESM2-LM","TaiESM1"/)
depths = (/3.2, 2.4, 1.6/)
do n = 0,dimsizes(depths)-1,1
  depth = depths(n)
  print(depth)
  do m = 0,dimsizes(listi)-1,1
    modelname = listi(m)
    print(modelname)

    ;Read data
    files = systemfunc("ls ../Data/cmip6_alt_m15_NH45_1850_2100/"+modelname+"_alt_monthly_hd.nc"); 1850-2100 monthly
    ;files = systemfunc("ls ../Data/cmip6_alt_m15_NH45_ssp245/"+modelname+"_alt_monthly_hd.nc"); 2015-2100 monthly
    f=addfile(files,"r")
    lat = f->lat
    lon = f->lon
    Alt = f->alt
    printVarSummary(Alt)

    ;>---------------------------------------------------------------------------------------<
        ;;                                Mask Greenland
    ;>---------------------------------------------------------------------------------------<
    f=addfile("/home/wangjx/Data/PCN_Domain_0.5x0.5.nc","r")
    green_mask = f->RCN_reg(270:359,:)
    
    ;Convert longitude from [-180, 180] to [0, 360]
    dx = 360 / 720  ; Longitude resolution for Alt
    hdx = dx / 2
    lon = fspan(-180 + hdx, 180 - hdx, 720)
    lon!0 = "lon"
    lon&lon = lon
    lon@units = "degrees_east"
    lon = lonFlip(lon) 
    lon = where(lon < 0, lon + 360, lon)

    green_mask := lonFlip(green_mask)
    printVarSummary(green_mask)
    green_mask_3d = conform_dims(dimsizes(Alt), green_mask, (/1, 2/))  ; Convert to time * lat * lon
    Alt := where(green_mask_3d.gt.-1, Alt, Alt@_FillValue)  ; Mask Greenland and oceans
    
    ;>---------------------------------------------------------------------------------------<
    ;                                0.5*0.5 -> 1.0*1.0
    ;>---------------------------------------------------------------------------------------<
    newlat = ispan(455, 895, 10) / 10.
    newlon = ispan(5, 3595, 10) / 10.
    newlat@units = "degrees_north"
    newlon@units = "degrees_east"
    ; Interpolation (area weighting)
    opt = True
    opt@critpc = 50
    alt_grid_1 = area_hi2lores(lon, lat, Alt, True, 1, newlon, newlat, opt)
    printVarSummary(alt_grid_1)
    
    ;>---------------------------------------------------------------------------------------<
    ;                   Calculate permafrost grid area for 0-190 longitude
    ;>---------------------------------------------------------------------------------------<
    alt_grid_1_1 = alt_grid_1(:,:,:190)  ; 0-190 longitude region
    printVarSummary(alt_grid_1_1)
    
    f = addfile("/home/wangjx/Data/sftlf_surf_frac_USGS_grid1.0_gn.nc","r")
    sftlf = f->sftlf({45:90}, {0:190.5})

    pi = 3.1415926
    r = 6372 * 10^3  ; Radius of Earth in meters

    nlat = dimsizes(newlat)
    nlon = dimsizes(sftlf(1,:))
    areacella = new((/nlat, nlon/), float)

    do i = 0, nlat - 1
      do j = 0, nlon - 1
        grid = 1.0  ; Grid resolution
        areacella(i,j) = abs(r * r * 1.0 / 180 * pi * (sin((newlat(i) + grid / 2.0) / 180 * pi) - sin((newlat(i) - grid / 2.0) / 180 * pi)))
      end do
    end do

    areacella_3d = conform_dims(dimsizes(alt_grid_1_1), areacella, (/1, 2/))  ; Convert to time * lat * lon
    ; Depth threshold
    areacella_mask = where(alt_grid_1_1 .le. depth, areacella_3d, 0)  
    sftlf_3d = conform_dims(dimsizes(alt_grid_1_1), sftlf, (/1, 2/))  ; Convert to time * lat * lon
    
    ; Calculate grid area
    alt_area = dim_sum_n(sftlf_3d * areacella_mask, (/1, 2/))

    delete(sftlf)
    delete(areacella)
    delete(areacella_3d)
    delete(sftlf_3d)
    delete(areacella_mask)

    ;>---------------------------------------------------------------------------------------<
    ;                  Calculate permafrost grid area for 350-360 longitude
    ;>---------------------------------------------------------------------------------------<
    alt_grid_1_2 = alt_grid_1(:,:,350:)  ; 350-360 longitude region
    sftlf = f->sftlf({45:90}, {350.5:})
    nlon = dimsizes(sftlf(1,:))
    areacella = new((/nlat, nlon/), float)

    do i = 0, nlat - 1
      do j = 0, nlon - 1
        grid = 1.0  ; Grid resolution
        areacella(i,j) = abs(r * r * 1.0 / 180 * pi * (sin((newlat(i) + grid / 2.0) / 180 * pi) - sin((newlat(i) - grid / 2.0) / 180 * pi)))
      end do
    end do

    areacella_3d = conform_dims(dimsizes(alt_grid_1_2), areacella, (/1, 2/))  ; Convert to time * lat * lon
    ; Depth threshold
    areacella_mask = where(alt_grid_1_2 .le. depth, areacella_3d, 0)
    sftlf_3d = conform_dims(dimsizes(alt_grid_1_2), sftlf, (/1, 2/))  ; Convert to time * lat * lon
    
    ; Calculate grid area
    alt_area2 = dim_sum_n(sftlf_3d * areacella_mask, (/1, 2/))

    ;>---------------------------------------------------------------------------------------<
        ;;                               Output data
    ;>---------------------------------------------------------------------------------------<

    fout = addfile("../Data/Permafrost_area_timeseries_ssp585/cmip6_alt"+toint(depth*100)+"cm_NH45_EA0_180_pfarea_1850_2100_monthly/"+modelname+"_pfarea_NH45_1850_2100monthly.nc","c") ; Monthly area
    fout2 = addfile("../Data/Permafrost_area_timeseries_ssp585/cmip6_alt"+toint(depth*100)+"cm_NH45_EA350_360_pfarea_1850_2100_monthly/"+modelname+"_pfarea_NH45_1850_2100monthly.nc","c") ; Monthly area
    ;fout = addfile("../Data/Permafrost_area_timeseries_ssp245/cmip6_alt"+toint(depth*100)+"cm_NH45_EA0_180_pfarea_2015_2100_monthly/"+modelname+"_pfarea_NH45_2015_2100monthly.nc","c") ; Monthly area
    ;fout2 = addfile("../Data/Permafrost_area_timeseries_ssp245/cmip6_alt"+toint(depth*100)+"cm_NH45_EA350_360_pfarea_2015_2100_monthly/"+modelname+"_pfarea_NH45_2015_2100monthly.nc","c") ; Monthly area    
    fout->pfarea = alt_area
    fout2->pfarea = alt_area2

    delete(f)
    delete(sftlf)
    delete(areacella)
    delete(areacella_3d)
    delete(sftlf_3d)
    delete(areacella_mask)

  end do
end do
end