#!/usr/bin/env python

"""
Use LIS-CABLE half-hour met output to make met file for offline CABLE

Source: generate_CABLE_netcdf_met_imp.py

Includes:

That's all folks.
"""

__author__    = "MU Mengyuan"
__email__     = "mu.mengyuan815@gmail.com"

import os
import sys
import glob
import pandas as pd
import numpy as np
import netCDF4 as nc
import datetime
from scipy.interpolate import interp1d
from scipy.interpolate import griddata

def generate_pixel_met(fname_in, fname_out, loc_lat, loc_lon,
                       CO2_val, css_val, cnsd_val, rho_val, org_val, watr_val, lat_name, lon_name):

    # # convect units
    # DEG_2_KELVIN = 273.15
    # SW_2_PAR     = 2.3
    # PAR_2_SW     = 1.0 / SW_2_PAR
    # HLFHR_2_SEC  = 1.0 / 1800.

    times, times_units, times_long_name, times_calendar = get_lis_times(fname_in)
    nsoil  = 6
    ndim   = 1
    ntimes = len(times)
    
    # create file and write global attributes
    f = nc.Dataset(fname_out, 'w', format='NETCDF4')
    f.description   = 'ERAI met data, created by MU Mengyuan'
    f.source        = fname_in
    f.history       = "Created by: %s" % (os.path.basename(__file__))
    f.creation_date = "%s" % (datetime.datetime.now())

    # set dimensions
    f.createDimension('time', ntimes) # None)
    f.createDimension('z', ndim)
    f.createDimension('y', ndim)
    f.createDimension('x', ndim)
    f.createDimension('soil_depth', nsoil)
    f.Conventions  = "CF-1.0"

    # create variables
    time           = f.createVariable('time', 'f8', ('time',)) 
    time.units     = times_units # "seconds since 2006-01-01 00:00:00"
    time.long_name = times_long_name # "time" ???
    time.calendar  = times_calendar # "standard" ???
    time[:]       = times[:] # ??? data format
    print(time)

    z = f.createVariable('z', 'f8', ('z',))
    z.long_name = "z"
    z.long_name = "z dimension"

    y = f.createVariable('y', 'f8', ('y',))
    y.long_name = "y"
    y.long_name = "y dimension"

    x = f.createVariable('x', 'f8', ('x',))
    x.long_name = "x"
    x.long_name = "x dimension"

    soil_depth = f.createVariable('soil_depth', 'f4', ('soil_depth',))
    soil_depth.long_name = "soil_depth"

    latitude = f.createVariable('latitude', 'f4', ('y', 'x',))
    latitude.units = "degrees_north"
    latitude.missing_value = -9999.
    latitude.long_name = "Latitude"

    longitude = f.createVariable('longitude', 'f4', ('y', 'x',))
    longitude.units = "degrees_east"
    longitude.missing_value = -9999.
    longitude.long_name = "Longitude"

    SWdown = f.createVariable('SWdown', 'f4', ('time', 'y', 'x',))
    SWdown.units = "W/m^2"
    SWdown.missing_value = -9999.
    SWdown.long_name = "Surface incident shortwave radiation"
    SWdown.CF_name = "surface_downwelling_shortwave_flux_in_air"

    Tair = f.createVariable('Tair', 'f4', ('time', 'z', 'y', 'x',))
    Tair.units = "K"
    Tair.missing_value = -9999.
    Tair.long_name = "Near surface air temperature"
    Tair.CF_name = "surface_temperature"

    Rainf = f.createVariable('Rainf', 'f4', ('time', 'y', 'x',))
    Rainf.units = "mm/s"
    Rainf.missing_value = -9999.
    Rainf.long_name = "Rainfall rate"
    Rainf.CF_name = "precipitation_flux"

    Qair = f.createVariable('Qair', 'f4', ('time', 'z', 'y', 'x',))
    Qair.units = "kg/kg"
    Qair.missing_value = -9999.
    Qair.long_name = "Near surface specific humidity"
    Qair.CF_name = "surface_specific_humidity"

    Wind = f.createVariable('Wind', 'f4', ('time', 'z', 'y', 'x',))
    Wind.units = "m/s"
    Wind.missing_value = -9999.
    Wind.long_name = "Scalar windspeed" 
    Wind.CF_name = "wind_speed"

    PSurf = f.createVariable('PSurf', 'f4', ('time', 'y', 'x',))
    PSurf.units = "Pa"
    PSurf.missing_value = -9999.
    PSurf.long_name = "Surface air pressure"
    PSurf.CF_name = "surface_air_pressure"

    LWdown = f.createVariable('LWdown', 'f4', ('time', 'y', 'x',))
    LWdown.units = "W/m^2"
    LWdown.missing_value = -9999.
    LWdown.long_name = "Surface incident longwave radiation"
    LWdown.CF_name = "surface_downwelling_longwave_flux_in_air"

    CO2 = f.createVariable('CO2air', 'f4', ('time', 'z', 'y', 'x',))
    CO2.units = "ppm"
    CO2.missing_value = -9999.
    CO2.long_name = ""
    CO2.CF_name = ""

    LAI = f.createVariable('LAI', 'f4', ('time', 'y', 'x'))
    LAI.setncatts({'long_name': u"Leaf Area Index",})

    elevation = f.createVariable('elevation', 'f4', ('y', 'x',))
    elevation.units = "m" 
    elevation.missing_value = -9999.
    elevation.long_name = "Site elevation above sea level" 

    iveg = f.createVariable('iveg', 'f4', ('y', 'x',))
    iveg.long_name = "vegetation type"
    iveg.units = "-"
    iveg.missing_value = -9999.0

    sand = f.createVariable('sand', 'f4', ('y', 'x',))
    sand.units = "-"
    sand.missing_value = -1.0

    clay = f.createVariable('clay', 'f4', ('y', 'x',))
    clay.units = "-"
    clay.missing_value = -1.0

    silt = f.createVariable('silt', 'f4', ('y', 'x',))
    silt.units = "-"
    silt.missing_value = -1.0

    rhosoil = f.createVariable('rhosoil', 'f4', ('y', 'x',))
    rhosoil.units = "kg m-3"
    rhosoil.long_name = "soil density"
    rhosoil.missing_value = -9999.0

    bch  = f.createVariable('bch', 'f4', ('y', 'x',))
    bch.units = "-"
    bch.long_name = "C and H B"
    bch.missing_value = -9999.0

    hyds = f.createVariable('hyds', 'f4', ('y', 'x',))
    hyds.units = "m s-1"
    hyds.long_name = "hydraulic conductivity at saturation"
    hyds.missing_value = -9999.0

    sucs = f.createVariable('sucs', 'f4', ('y', 'x',))
    sucs.units = "m"
    sucs.long_name = "matric potential at saturation"
    sucs.missing_value = -9999.0

    ssat = f.createVariable('ssat', 'f4', ('y', 'x',))
    ssat.units = "m3 m-3"
    ssat.long_name = "volumetric water content at saturation"
    ssat.missing_value = -9999.0

    swilt= f.createVariable('swilt', 'f4', ('y', 'x',))
    swilt.units = "m3 m-3"
    swilt.long_name = "wilting point"
    swilt.missing_value = -9999.0

    sfc  = f.createVariable('sfc', 'f4', ('y', 'x',))
    sfc.units = "m3 m-3"
    sfc.long_name = "field capcacity"
    sfc.missing_value = -9999.0

    css  = f.createVariable('css', 'f4', ('y', 'x',))
    css.units = "kJ kg-1 K-1"
    css.long_name = "soil specific heat capacity"
    css.missing_value = -9999.0

    cnsd = f.createVariable('cnsd', 'f4', ('y', 'x',))
    cnsd.units = "W m-1 K-1"
    cnsd.long_name = "thermal conductivity of dry soil"
    cnsd.missing_value = -9999.0

    # 3-Dimension variables
    sand_vec = f.createVariable('sand_vec', 'f4', ('soil_depth', 'y', 'x',))
    sand_vec.units = "-"
    sand_vec.missing_value = -1.0

    clay_vec = f.createVariable('clay_vec', 'f4', ('soil_depth', 'y', 'x',))
    clay_vec.units = "-"
    clay_vec.missing_value = -1.0

    silt_vec = f.createVariable('silt_vec', 'f4', ('soil_depth', 'y', 'x',))
    silt_vec.units = "-"
    silt_vec.missing_value = -1.0

    org_vec  = f.createVariable('org_vec', 'f4', ('soil_depth', 'y', 'x',))
    org_vec.units = "-"
    org_vec.missing_value = -1.0

    rhosoil_vec = f.createVariable('rhosoil_vec', 'f4', ('soil_depth', 'y', 'x',))
    rhosoil_vec.units = "kg m-3"
    rhosoil_vec.long_name = "soil density"
    rhosoil_vec.missing_value = -9999.0

    bch_vec  = f.createVariable('bch_vec', 'f4', ('soil_depth', 'y', 'x',))
    bch_vec.units = "-"
    bch_vec.long_name = "C and H B"
    bch_vec.missing_value = -9999.0

    hyds_vec = f.createVariable('hyds_vec', 'f4', ('soil_depth', 'y', 'x',))
    hyds_vec.units = "mm s-1"
    hyds_vec.long_name = "hydraulic conductivity at saturation"
    hyds_vec.missing_value = -9999.0

    sucs_vec = f.createVariable('sucs_vec', 'f4', ('soil_depth', 'y', 'x',))
    sucs_vec.units = "m"
    sucs_vec.long_name = "matric potential at saturation"
    sucs_vec.missing_value = -9999.0

    ssat_vec = f.createVariable('ssat_vec', 'f4', ('soil_depth', 'y', 'x',))
    ssat_vec.units = "m3 m-3"
    ssat_vec.long_name = "volumetric water content at saturation"
    ssat_vec.missing_value = -9999.0

    swilt_vec= f.createVariable('swilt_vec', 'f4', ('soil_depth', 'y', 'x',))
    swilt_vec.units = "m3 m-3"
    swilt_vec.long_name = "wilting point"
    swilt_vec.missing_value = -9999.0

    sfc_vec  = f.createVariable('sfc_vec', 'f4', ('soil_depth', 'y', 'x',))
    sfc_vec.units = "m3 m-3"
    sfc_vec.long_name = "field capcacity"
    sfc_vec.missing_value = -9999.0

    css_vec  = f.createVariable('css_vec', 'f4', ('soil_depth', 'y', 'x',))
    css_vec.units = "kJ kg-1 K-1"
    css_vec.long_name = "soil specific heat capacity"
    css_vec.missing_value = -9999.0

    cnsd_vec = f.createVariable('cnsd_vec', 'f4', ('soil_depth', 'y', 'x',))
    cnsd_vec.units = "W m-1 K-1"
    cnsd_vec.long_name = "thermal conductivity of dry soil"
    cnsd_vec.missing_value = -9999.0

    watr = f.createVariable('watr', 'f4', ('soil_depth', 'y', 'x',))
    watr.units = "m3 m-3"
    watr.long_name = "residual water content of the soil"
    watr.missing_value = -9999.0

    # write data to file
    x[:] = ndim
    y[:] = ndim
    z[:] = ndim
    soil_depth[:] = np.arange(0,nsoil,1)

    # get landinfo values from gridinfo
    css[0,0]       = css_val
    cnsd[0,0]      = cnsd_val
    rhosoil[0,0]   = rho_val
    org_vec[:,0,0] = org_val
    watr[:,0,0]    = watr_val
    CO2[:,0,0]     = CO2_val     

    # get varibale values from lis
    rain,tair,qair,psurf,swdown,lwdown,wind,lai \
                                   = get_met_from_lis(fname_in,loc_lat,loc_lon,lat_name, lon_name)

    latitude[0]   = loc_lat
    longitude[0]  = loc_lon
    SWdown[:,0,0] = swdown
    # print("SWdown = %f" % swdown)
    LWdown[:,0,0] = lwdown
    # print("LWdown = %f" % lwdown)
    PSurf[:,0,0]  = psurf
    # print("PSurf = %f" % psurf)
    Rainf[:,0,0]  = rain
    # print("Rainf = %f" % rain)
    Tair[:,0,0,0] = tair
    # print("Tair = %f" % tair)
    Qair[:,0,0,0] = qair
    # print("Qair = %f" % qair)
    Wind[:,0,0,0] = wind
    # print("Wind = %f" % wind)
    LAI[:,0,0]    = lai
    # print("LAI = %f" % lai)

    # get landinfo values from lis
    mask = mask_by_lat_lon(fname_in[0], loc_lat, loc_lon, lat_name, lon_name)

    para_cable = nc.Dataset(fname_in[0], 'r')
    elevation[0,0] = para_cable.variables['Elevation_inst'][0][mask]
    # print("elevation = %f" % elevation[:])
    iveg[0,0] = para_cable.variables['Landcover_inst'][0][mask]
    # print("iveg = %f" % iveg[:])
    sand[0,0] = para_cable.variables['SandFrac_inst'][0][mask]
    # print("isandveg = %f" % sand[:])
    clay[0,0] = para_cable.variables['ClayFrac_inst'][0][mask]
    # print("clay = %f" % clay[:])
    silt[0,0] = para_cable.variables['SiltFrac_inst'][0][mask]
    # print("silt = %f" % silt[:])
    bch[0,0]  = para_cable.variables['Bch_inst'][0][mask]
    # print("bch = %f" % bch[:])
    hyds[0,0] = para_cable.variables['Hyds_inst'][0][mask]
    # print("hyds = %f" % hyds[:])
    sucs[0,0] = para_cable.variables['Sucs_inst'][0][mask]
    # print("sucs = %f" % sucs[:])
    ssat[0,0] = para_cable.variables['SoilSat_inst'][0][mask]
    # print("ssat = %f" % ssat[:])
    swilt[0,0]= para_cable.variables['SoilWiltPt_inst'][0][mask]
    # print("swilt = %f" % swilt[:])
    sfc[0,0]  = para_cable.variables['SoilFieldCap_inst'][0][mask]
    # print("sfc = %f" % sfc[:])

    # 3-Dimension variables
    sand_vec[:,0,0]    = np.repeat(sand[0,0], 6)
    clay_vec[:,0,0]    = np.repeat(clay[0,0], 6)
    silt_vec[:,0,0]    = np.repeat(silt[0,0], 6)
    rhosoil_vec[:,0,0] = np.repeat(rhosoil[0,0], 6)
    bch_vec[:,0,0]     = np.repeat(bch[0,0], 6)
    hyds_vec[:,0,0]    = np.repeat(hyds[0,0]*1000, 6)
    sucs_vec[:,0,0]    = np.repeat(sucs[0,0], 6)
    ssat_vec[:,0,0]    = np.repeat(ssat[0,0], 6)
    swilt_vec[:,0,0]   = np.repeat(swilt[0,0], 6)
    sfc_vec[:,0,0]     = np.repeat(sfc[0,0], 6)
    css_vec[:,0,0]     = np.repeat(css[0,0], 6)
    cnsd_vec[:,0,0]    = np.repeat(cnsd[0,0], 6)

    f.close()

def get_lis_times(fname_in):

    for i,fn_in in enumerate(fname_in):
        if i == 0:  
            cable           = nc.Dataset(fname_in[i], 'r')
            tmp             = cable.variables['time']
            times           = tmp[:]
            times_units     = tmp.units 
            times_long_name = tmp.standard_name
            times_calendar  = tmp.calendar
            cable.close()
        else:
            cable = nc.Dataset(fn_in, 'r')
            times.append(cable.variables['time'][:])
            cable.close()

    print(times)
    print(times_units)
    print(times_long_name)
    print(times_calendar)      

    return times,times_units,times_long_name,times_calendar

def get_met_from_lis(fname_in,loc_lat,loc_lon, lat_name, lon_name): 

    """
    read met fields from LIS-CABLE output
    """

    print("carry on read_cable_var")

    if len(fname_in) == 1:

        cable  = nc.Dataset(fname_in[0], 'r')

        mask   = mask_by_lat_lon(fname_in[0], loc_lat, loc_lon, lat_name, lon_name)
        ntimes = len(cable.variables['time'])
        print(ntimes)
        mask_multi = [ mask ] * ntimes
        print(np.shape(mask_multi))
        rain   = cable.variables['Rainf_f_inst'][:][mask_multi]#.filled(-9999.)
        print(np.shape(rain))
        tair   = cable.variables['Tair_f_inst'][:][mask_multi]
        print(tair)
        qair   = cable.variables['Qair_f_inst'][:][mask_multi]
        print(qair)
        psurf  = cable.variables['Psurf_f_inst'][:][mask_multi]
        print(psurf)
        swdown = cable.variables['SWdown_f_inst'][:][mask_multi]
        print(swdown)
        lwdown = cable.variables['LWdown_f_inst'][:][mask_multi]
        print(lwdown)
        wind   = cable.variables['Wind_f_inst'][:][mask_multi]
        print(wind)
        lai    = cable.variables['LAI_inst'][:][mask_multi]
        print(lai)
        cable.close()        
    
    else: 
        rain     = []
        tair     = []
        qair     = []
        psurf    = []
        swdown   = []
        lwdown   = []
        wind     = []
        lai      = []

        for i, fn_in in enumerate(fname_in):

            cable = nc.Dataset(fn_in, 'r')
            ntimes = len(cable.variables['times'])

            mask   = mask_by_lat_lon(fn_in, loc_lat, loc_lon, lat_name, lon_name)
            mask_multi = [ mask ] * ntimes

            if i == 0:
                lat      = cable.variables['lat'][mask]
                lon      = cable.variables['lon'][mask]

            tmp = cable.variables['Rainf_f_inst'][mask_multi].filled(-9999.)
            rain.append(tmp)
            tmp = cable.variables['Tair_f_inst'][mask_multi].filled(-9999.)
            tair.append(tmp)
            tmp = cable.variables['Qair_f_inst'][mask_multi].filled(-9999.)
            qair.append(tmp)
            tmp = cable.variables['Psurf_f_inst'][mask_multi].filled(-9999.)
            psurf.append(tmp)
            tmp = cable.variables['SWdown_f_inst'][mask_multi].filled(-9999.)
            swdown.append(tmp)
            tmp = cable.variables['LWdown_f_inst'][mask_multi].filled(-9999.)
            lwdown.append(tmp)
            tmp = cable.variables['Wind_f_inst'][mask_multi].filled(-9999.)
            wind.append(tmp)
            tmp = cable.variables['LAI_inst'][mask_multi].filled(-9999.)
            lai.append(tmp)
            cable.close()        

    print(rain[:])

    return rain,tair,qair,psurf,swdown,lwdown,wind,lai

def mask_by_lat_lon(file_path, loc_lat, loc_lon, lat_name, lon_name):

    file = nc.Dataset(file_path, mode='r')
    lat  = file.variables[lat_name][:]
    lon  = file.variables[lon_name][:]
    print(lat)
    print(lon)
    if len(np.shape(lat)) == 1:
        print("len(np.shape(lat)) == 1")
        lat_spc = lat[1] - lat[0]
        lon_spc = lon[1] - lon[0]
        lons, lats = np.meshgrid(lon, lat)
        mask  = (lats > (loc_lat - lat_spc/2)) & (lats < (loc_lat + lat_spc/2)) & (lons > (loc_lon - lon_spc/2)) & (lons < (loc_lon + lon_spc/2))
    elif len(np.shape(lat)) == 2:
        print("len(np.shape(lat)) == 2")
        ### caution: lat=100, lon=100 is a random pixel, lis run over a small domain may not have such a point
        lat_spc = lat[100,100] - lat[99,100]
        lon_spc = lon[100,100] - lon[100,99]
        print(lat_spc)
        print(lon_spc)
        ### caution: due to irregular space in lis, using lat/lon +lat/lon_spc/2 may includes more than 1 pixel. 
        ### I therefore let the space divied by 2.1 rather than 2
        mask  = (lat > (loc_lat - lat_spc/2.1)) & (lat < (loc_lat + lat_spc/2.1)) & (lon > (loc_lon - lon_spc/2.1)) & (lon < (loc_lon + lon_spc/2.1))
    return mask

def timestep_interpolation(var,dels):

    '''
    Haven't checked whether this spell works
    '''

    tts_05hr = len(var)
    nts      = int(dels/1800)
    tts      = int(tts_05hr/nts)
    print("tts_05hr %s" % tts_05hr)
    print("nts %s" % nts)
    print("tts %s" % tts)
    var_tmp   = np.zeros(tts)
    for j in np.arange(0,tts):
        var_tmp[j]    = np.average(var[j*nts:(j+1)*nts])

    print(var_tmp)
    return(var_tmp)

if __name__ == "__main__":

    ### use lis output files ###
    path_in      = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/ctl_15Sep/LIS_output/"
    fname_in     = [ path_in +"LIS.CABLE.200001-200012.nc" ]
    fname_out    = "./nc_files/ERAI_05hr_pixel_met_from_LIS-CABLE_satfrac_fixed.nc"
    fname_grdinf = "/g/data/w35/mm3972/model/cable/src/CABLE-AUX/offline/gridinfo_AWAP_OpenLandMap_ELEV_DLCM_fix_6layer_uniform.nc"

    # delete output file if it exists.
    os.remove(fname_out)
    
    # lat  lon 
    loc_lat = -34. #-24.255707     
    loc_lon = 145  #135.95001     

    lat_name = "latitude"
    lon_name = "longitude"

    mask     = mask_by_lat_lon(fname_grdinf, loc_lat, loc_lon, lat_name, lon_name)
    grdinf   = nc.Dataset(fname_grdinf, 'r')   
    css_val  = grdinf.variables['css'][:][mask]
    cnsd_val = grdinf.variables['cnsd'][:][mask]
    rho_val  = grdinf.variables['rhosoil'][:][mask] 
    org_val  = grdinf.variables['organic'][:][mask]

    watr_val = 0.01   # from lis screen print
    CO2_val  = 369.55 # 2000


    lat_name = "lat"
    lon_name = "lon"
    generate_pixel_met(fname_in, fname_out, loc_lat, loc_lon, CO2_val, css_val, cnsd_val, rho_val, org_val, watr_val, lat_name, lon_name)

