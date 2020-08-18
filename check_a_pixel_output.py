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
        secs += dels

    # create file and write global attributes
    f = nc.Dataset(out_fname, 'w', format='NETCDF4')
    f.description   = 'check soil temperature and moisture, created by MU Mengyuan'
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
    time.units     = "seconds since 1980-01-01 00:00:00"
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

    SoilMoist = f.createVariable('SoilMoist_tavg', 'f4', ('time','soil_depth', 'y', 'x',))
    SoilMoist.units = "m3/m3"
    SoilMoist.missing_value = -9999.
    SoilMoist.long_name = "Soil moisture"
    SoilMoist.CF_name = "Soil moisture"

    SoilTemp = f.createVariable('SoilTemp_tavg', 'f4', ('time', 'soil_depth','y', 'x',))
    SoilTemp.units = "K"
    SoilTemp.missing_value = -9999.
    SoilTemp.long_name = "Soil temperature"
    SoilTemp.CF_name = "Soil temperature"

    # write data to file
    x[:] = ndim
    y[:] = ndim
    z[:] = ndim

    soil_depth[:] = np.arange(0,nsoil,1)
    time[:] = times

    lat,lon,soilmoist,soiltemp = get_met_input(input_fname,loc_lat,loc_lon,dels)

    print(lat)
    print(lon)

    latitude[0]  = lat
    longitude[0] = lon
    SoilMoist[:,:,0,0] = soilmoist
    SoilTemp[:,:,0,0] = soiltemp

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

            soilmoist = cable.variables['SoilMoist_tavg'][:,:,loc_lat,loc_lon].filled(-9999.)
            soiltemp  = cable.variables['SoilTemp_tavg'][:,:,loc_lat,loc_lon].filled(-9999.)

        else:
            soilmoist = np.concatenate((soilmoist,cable.variables['SoilMoist_tavg'][:,:,loc_lat,loc_lon].filled(-9999.)))
            soiltemp  = np.concatenate((soiltemp,cable.variables['SoilTemp_tavg'][:,:,loc_lat,loc_lon].filled(-9999.)))

        cable.close()

    return lat,lon,soilmoist,soiltemp;

if __name__ == "__main__":

    path = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/ERAI_ctl/LIS_output/"
    input_fname = [path+"LIS.CABLE.1980010100.d01.nc",path+"LIS.CABLE.1980020100.d01.nc",
                   path+"LIS.CABLE.1980030100.d01.nc",path+"LIS.CABLE.1980040100.d01.nc",
                   path+"LIS.CABLE.1980050100.d01.nc",path+"LIS.CABLE.1980060100.d01.nc",
                   path+"LIS.CABLE.1980070100.d01.nc",path+"LIS.CABLE.1980080100.d01.nc",
                   path+"LIS.CABLE.1980090100.d01.nc",path+"LIS.CABLE.1980100100.d01.nc",
                   path+"LIS.CABLE.1980110100.d01.nc",path+"LIS.CABLE.1980120100.d01.nc",
                   ]
    out_fname = "./nc_files/pixel_output_lat_-306_lon_1256.nc"

    loc_lat = 63
    loc_lon = 53

    dels = 3600*24
    main(input_fname, out_fname, loc_lat, loc_lon, dels)
