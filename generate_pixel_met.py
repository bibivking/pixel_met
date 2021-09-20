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

def main(input_fname, out_fname, loc_lat, loc_lon, dels):

    print(input_fname[0])
    print(out_fname)
    DEG_2_KELVIN = 273.15
    SW_2_PAR     = 2.3
    PAR_2_SW     = 1.0 / SW_2_PAR
    HLFHR_2_SEC  = 1.0 / 1800.

    DOY   = 365
    nsoil = 6
    ndim  = 1
    nsoil = nsoil
    n_timesteps = int(DOY*(24*60*60/dels))
    print(n_timesteps)
    times = []
    secs  = 9*3600 + dels
    for i in range(n_timesteps):
        times.append(secs)
        secs += dels

    # create file and write global attributes
    f = nc.Dataset(out_fname, 'w', format='NETCDF4')
    f.description   = 'ERAI met data, created by MU Mengyuan'
    f.source        = input_fname
    f.history       = "Created by: %s" % (os.path.basename(__file__))
    f.creation_date = "%s" % (datetime.datetime.now())

    # set dimensions
    f.createDimension('time', None)
    f.createDimension('z', ndim)
    f.createDimension('y', ndim)
    f.createDimension('x', ndim)
    f.createDimension('soil_depth', nsoil)
    f.Conventions  = "CF-1.0"

    # create variables
    time           = f.createVariable('time', 'f8', ('time',))
    time.units     = "seconds since 2018-01-01 00:00:00"
    time.long_name = "time"
    time.calendar  = "standard"

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
    Wind.long_name = "Scalar windspeed" ;
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
    elevation.units = "m" ;
    elevation.missing_value = -9999.
    elevation.long_name = "Site elevation above sea level" ;

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
    time[:] = times

    lat,lon,rain,tair,qair,psurf,swdown,lwdown,wind,lai = \
                                    get_met_input(input_fname,loc_lat,loc_lon,dels)

    print(lat)
    print(lon)

    latitude[0]  = lat
    longitude[0] = lon
    print(latitude[0])
    print(longitude[0])
    SWdown[:,0,0] = swdown
    LWdown[:,0,0] = lwdown
    PSurf[:,0,0]  = psurf
    Rainf[:,0,0]  = rain
    Tair[:,0,0,0] = tair
    Qair[:,0,0,0] = qair
    Wind[:,0,0,0] = wind
    CO2[:,0,0]    = 335. # 1978
    LAI[:,0,0]    = lai


    para_cable = nc.Dataset(input_fname[0], 'r')
    elevation[0,0] = para_cable.variables['Elevation_inst'][0,loc_lat,loc_lon]
    iveg[0,0] = 5 # veg type 7 in LIS is 5 in CABLE
                # para_cable.variables['Landcover_inst'][0,loc_lat,loc_lon]
    sand[0,0] = para_cable.variables['SandFrac_inst'][0,loc_lat,loc_lon]
    clay[0,0] = para_cable.variables['ClayFrac_inst'][0,loc_lat,loc_lon]
    silt[0,0] = para_cable.variables['SiltFrac_inst'][0,loc_lat,loc_lon]
    rhosoil[0,0] = para_cable.variables['rhosoil_inst'][0,loc_lat,loc_lon]
    bch[0,0]  = para_cable.variables['bch_inst'][0,loc_lat,loc_lon]
    hyds[0,0] = para_cable.variables['hyds_inst'][0,loc_lat,loc_lon]
    sucs[0,0] = para_cable.variables['sucs_inst'][0,loc_lat,loc_lon]
    ssat[0,0] = para_cable.variables['Porosity_inst'][0,loc_lat,loc_lon]
    swilt[0,0]= para_cable.variables['swilt_inst'][0,loc_lat,loc_lon]
    sfc[0,0]  = para_cable.variables['sfc_inst'][0,loc_lat,loc_lon]
    css[0,0]  = para_cable.variables['css_inst'][0,loc_lat,loc_lon]
    cnsd[0,0] = 0.277841448783875
    # 3-Dimension variables
    sand_vec[:,0,0]    = np.repeat(sand[0,0], 6)
    clay_vec[:,0,0]    = np.repeat(clay[0,0], 6)
    silt_vec[:,0,0]    = np.repeat(silt[0,0], 6)
    org_vec[:,0,0]     = [0,0,0,0,0,0]
    rhosoil_vec[:,0,0] = np.repeat(rhosoil[0,0], 6)
    bch_vec[:,0,0]     = np.repeat(bch[0,0], 6)
    hyds_vec[:,0,0]    = np.repeat(hyds[0,0]*1000, 6)
    sucs_vec[:,0,0]    = np.repeat(sucs[0,0], 6)
    ssat_vec[:,0,0]    = np.repeat(ssat[0,0], 6)
    swilt_vec[:,0,0]   = np.repeat(swilt[0,0], 6)
    sfc_vec[:,0,0]     = np.repeat(sfc[0,0], 6)
    css_vec[:,0,0]     = np.repeat(css[0,0], 6)
    cnsd_vec[:,0,0]    = np.repeat(cnsd[0,0], 6)
    watr[:,0,0]        = [0.05, 0.05, 0.05, 0.05, 0.05, 0.05]

    f.close()

def get_met_input(input_fname,loc_lat,loc_lon,dels):

    """
    read met fields from LIS-CABLE output
    """

    print("carry on read_cable_var")

    for month in np.arange(0,12,1):
        print(month)
        cable = nc.Dataset(input_fname[month], 'r')

        if month == 0:
            lat      = cable.variables['lat'][0,loc_lat,loc_lon]
            lon      = cable.variables['lon'][0,loc_lat,loc_lon]
            landmask = cable.variables['Landmask_inst'][0,loc_lat,loc_lon]
            landcover= cable.variables['Landcover_inst'][0,loc_lat,loc_lon]
            print(lat)
            print(lon)
            print(landmask)
            print(landcover)

            rain     = cable.variables['Rainf_f_inst'][:,loc_lat,loc_lon].filled(-9999.)
            tair     = cable.variables['Tair_f_inst'][:,loc_lat,loc_lon].filled(-9999.)
            qair     = cable.variables['Qair_f_inst'][:,loc_lat,loc_lon].filled(-9999.)
            psurf    = cable.variables['Psurf_f_inst'][:,loc_lat,loc_lon].filled(-9999.)
            swdown   = cable.variables['SWdown_f_inst'][:,loc_lat,loc_lon].filled(-9999.)
            lwdown   = cable.variables['LWdown_f_inst'][:,loc_lat,loc_lon].filled(-9999.)
            wind     = cable.variables['Wind_f_inst'][:,loc_lat,loc_lon].filled(-9999.)
            lai      = cable.variables['LAI_inst'][:,loc_lat,loc_lon].filled(-9999.)

        else:
            rain     = np.concatenate((rain,cable.variables['Rainf_f_inst'][:,loc_lat,loc_lon].filled(-9999.)))
            tair     = np.concatenate((tair,cable.variables['Tair_f_inst'][:,loc_lat,loc_lon].filled(-9999.)))
            qair     = np.concatenate((qair,cable.variables['Qair_f_inst'][:,loc_lat,loc_lon].filled(-9999.)))
            psurf    = np.concatenate((psurf,cable.variables['Psurf_f_inst'][:,loc_lat,loc_lon].filled(-9999.)))
            swdown   = np.concatenate((swdown,cable.variables['SWdown_f_inst'][:,loc_lat,loc_lon].filled(-9999.)))
            lwdown   = np.concatenate((lwdown,cable.variables['LWdown_f_inst'][:,loc_lat,loc_lon].filled(-9999.)))
            wind     = np.concatenate((wind,cable.variables['Wind_f_inst'][:,loc_lat,loc_lon].filled(-9999.)))
            lai      = np.concatenate((lai,cable.variables['LAI_inst'][:,loc_lat,loc_lon].filled(-9999.)))

        cable.close()
    print(rain[:])

    if dels != 1800:
        tts_05hr = len(rain)
        nts      = int(dels/1800)
        tts      = int(tts_05hr/nts)
        print("tts_05hr %s" % tts_05hr)
        print("nts %s" % nts)
        print("tts %s" % tts)
        rain_tmp   = np.zeros(tts)
        tair_tmp   = np.zeros(tts)
        qair_tmp   = np.zeros(tts)
        psurf_tmp  = np.zeros(tts)
        swdown_tmp = np.zeros(tts)
        lwdown_tmp = np.zeros(tts)
        wind_tmp   = np.zeros(tts)
        lai_tmp    = np.zeros(tts)
        for j in np.arange(0,tts):
            rain_tmp[j]    = np.average(rain[j*nts:(j+1)*nts])
            tair_tmp[j]    = np.average(tair[j*nts:(j+1)*nts])
            qair_tmp[j]    = np.average(qair[j*nts:(j+1)*nts])
            psurf_tmp[j]   = np.average(psurf[j*nts:(j+1)*nts])
            swdown_tmp[j]  = np.average(swdown[j*nts:(j+1)*nts])
            lwdown_tmp[j]  = np.average(lwdown[j*nts:(j+1)*nts])
            wind_tmp[j]    = np.average(wind[j*nts:(j+1)*nts])
            lai_tmp[j]     = np.average(lai[j*nts:(j+1)*nts])

        rain   = rain_tmp
        tair   = tair_tmp
        qair   = qair_tmp
        psurf  = psurf_tmp
        swdown = swdown_tmp
        lwdown = lwdown_tmp
        wind   = wind_tmp
        lai    = lai_tmp

    print(rain[:])


    return lat,lon,rain,tair,qair,psurf,swdown,lwdown,wind,lai;

if __name__ == "__main__":

    path = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/ERAI_check_one_pixel/LIS_output_backup/"
    input_fname = [path+"LIS.CABLE.2018010100.d01.nc",path+"LIS.CABLE.2018020100.d01.nc",
                   path+"LIS.CABLE.2018030100.d01.nc",path+"LIS.CABLE.2018040100.d01.nc",
                   path+"LIS.CABLE.2018050100.d01.nc",path+"LIS.CABLE.2018060100.d01.nc",
                   path+"LIS.CABLE.2018070100.d01.nc",path+"LIS.CABLE.2018080100.d01.nc",
                   path+"LIS.CABLE.2018090100.d01.nc",path+"LIS.CABLE.2018100100.d01.nc",
                   path+"LIS.CABLE.2018110100.d01.nc",path+"LIS.CABLE.2018120100.d01.nc",
                   ]
    out_fname = "./nc_files/ERAI_05hr_pixel_met_-29_138.nc"

    # lat -24.255707 lon 135.95001
    loc_lat = 20 # 6
    loc_lon = 6 # 20

    #Latitude:   -29.65829, Longitude:    138.3852
    loc_lat = 1 # 6
    loc_lon = 14 # 20

    # # lat = -33.604504 lon = 150.60345 PFT: 2-evergreen broadleaf
    # loc_lat = 47
    # loc_lon = 144

    dels = 3600*0.5
    main(input_fname, out_fname, loc_lat, loc_lon, dels)
