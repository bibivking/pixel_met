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
from scipy.signal import savgol_filter

def main(input_fname, out_fname):

    print(input_fname[0])
    print(out_fname)
    DEG_2_KELVIN = 273.15
    SW_2_PAR     = 2.3
    PAR_2_SW     = 1.0 / SW_2_PAR
    HLFHR_2_SEC  = 1.0 / 1800.

    dels  = 1800
    DOY   = 366
    nsoil = 6
    ndim  = 1
    nsoil = nsoil
    n_timesteps = int(DOY*(24*60*60/dels))
    print(n_timesteps)
    times = []
    secs  = dels
    for i in range(n_timesteps):
        times.append(secs)
        secs += 1800.

    # create file and write global attributes
    f = nc.Dataset(out_fname, 'w', format='NETCDF4')
    f.description   = 'Princeton met data, created by MU Mengyuan'
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
    time.units     = "seconds since 2008-01-01 00:00:00"
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

    # write data to file
    x[:] = ndim
    y[:] = ndim
    z[:] = ndim

    soil_depth[:] = np.arange(0,nsoil,1)
    time[:] = times

    loc_lat = 40
    loc_lon = 140
    lat,lon,rain,tair,qair,psurf,swdown,lwdown,wind,lai = \
                                    get_met_input(input_fname,loc_lat,loc_lon)

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
    CO2[:,0,0]    = 350.
    LAI[:,0,0]    = lai

    f.close()

def get_met_input(input_fname,loc_lat,loc_lon):

    """
    read met fields from LIS-CABLE output
    """

    print("carry on read_cable_var")

    for month in np.arange(0,12,1):
        print(month)
        cable = nc.Dataset(input_fname[month], 'r')
        # rain  =[]
        # tair  =[]
        # qair  =[]
        # psurf =[]
        # swdown=[]
        # lwdown=[]
        # wind  =[]
        # lai   =[]

        if month == 0:
            lat      = cable.variables['lat'][0,loc_lat,loc_lon]
            lon      = cable.variables['lon'][0,loc_lat,loc_lon]
            landmask = cable.variables['Landmask_inst'][0,loc_lat,loc_lon]
            landcover= cable.variables['Landcover_inst'][0,loc_lat,loc_lon]
            print(lat)
            print(lon)
            print(landmask)
            print(landcover)

            rain     = cable.variables['Rainf_f_tavg'][:,loc_lat,loc_lon].filled(-9999.)
            tair     = cable.variables['Tair_f_tavg'][:,loc_lat,loc_lon].filled(-9999.)
            qair     = cable.variables['Qair_f_tavg'][:,loc_lat,loc_lon].filled(-9999.)
            psurf    = cable.variables['Psurf_f_tavg'][:,loc_lat,loc_lon].filled(-9999.)
            swdown   = cable.variables['SWdown_f_tavg'][:,loc_lat,loc_lon].filled(-9999.)
            lwdown   = cable.variables['LWdown_f_tavg'][:,loc_lat,loc_lon].filled(-9999.)
            wind     = cable.variables['Wind_f_tavg'][:,loc_lat,loc_lon].filled(-9999.)
            lai      = cable.variables['LAI_inst'][:,loc_lat,loc_lon].filled(-9999.)

        else:
            rain     = np.concatenate((rain,cable.variables['Rainf_f_tavg'][:,loc_lat,loc_lon].filled(-9999.)))
            tair     = np.concatenate((tair,cable.variables['Tair_f_tavg'][:,loc_lat,loc_lon].filled(-9999.)))
            qair     = np.concatenate((qair,cable.variables['Qair_f_tavg'][:,loc_lat,loc_lon].filled(-9999.)))
            psurf    = np.concatenate((psurf,cable.variables['Psurf_f_tavg'][:,loc_lat,loc_lon].filled(-9999.)))
            swdown   = np.concatenate((swdown,cable.variables['SWdown_f_tavg'][:,loc_lat,loc_lon].filled(-9999.)))
            lwdown   = np.concatenate((lwdown,cable.variables['LWdown_f_tavg'][:,loc_lat,loc_lon].filled(-9999.)))
            wind     = np.concatenate((wind,cable.variables['Wind_f_tavg'][:,loc_lat,loc_lon].filled(-9999.)))
            lai      = np.concatenate((lai,cable.variables['LAI_inst'][:,loc_lat,loc_lon].filled(-9999.)))

        cable.close()

    print(rain)


    return lat,lon,rain,tair,qair,psurf,swdown,lwdown,wind,lai;


if __name__ == "__main__":

    path = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/Princeton_ctl_2008_05hr/LIS_output/"
    input_fname = [path+"LIS.CABLE.2008010100.d01.nc",path+"LIS.CABLE.2008020100.d01.nc",
                   path+"LIS.CABLE.2008030100.d01.nc",path+"LIS.CABLE.2008040100.d01.nc",
                   path+"LIS.CABLE.2008050100.d01.nc",path+"LIS.CABLE.2008060100.d01.nc",
                   path+"LIS.CABLE.2008070100.d01.nc",path+"LIS.CABLE.2008080100.d01.nc",
                   path+"LIS.CABLE.2008090100.d01.nc",path+"LIS.CABLE.2008100100.d01.nc",
                   path+"LIS.CABLE.2008110100.d01.nc",path+"LIS.CABLE.2008120100.d01.nc",
                   ]
    out_fname = "pixel_met.nc"

    main(input_fname, out_fname)
