#!/usr/bin/env python

"""

"""

__author__    = "MU Mengyuan"

import os
import sys
import glob
import pandas as pd
import numpy as np
import netCDF4 as nc
import datetime
from scipy.interpolate import interp1d
from scipy.interpolate import griddata

def main(fname_05hr, fname_05hr_cldstrt, fname_3hr):

    loc_lat = 54
    loc_lon = 149

    f_05hr         = nc.Dataset(fname_05hr, "r")
    f_05hr_cldstrt = nc.Dataset(fname_05hr_cldstrt, "r")
    f_3hr          = nc.Dataset(fname_3hr, 'r')

    print(f_3hr.variables['x'][loc_lon])
    print(f_3hr.variables['y'][loc_lat])
    print(f_05hr.variables['x'][0])
    print(f_05hr.variables['y'][0])

    var_05hr          = pd.DataFrame(f_05hr.variables['Rainf'][:,0,0]*60.*60.*24., columns=['Rainf'])
    var_05hr['Evap']  = f_05hr.variables['Evap'][:,0,0]*60.*60.*24.
    var_05hr['TVeg']  = f_05hr.variables['TVeg'][:,0,0]*60.*60.*24.
    var_05hr['ESoil'] = f_05hr.variables['ESoil'][:,0,0]*60.*60.*24.
    var_05hr['ECanop']= f_05hr.variables['ECanop'][:,0,0]*60.*60.*24.
    var_05hr['Qs']    = f_05hr.variables['Qs'][:,0,0]*60.*60.*24.
    var_05hr['Qsb']   = f_05hr.variables['Qsb'][:,0,0]*60.*60.*24.
    var_05hr['SM1']   = f_05hr.variables['SoilMoist'][:,0,0,0]
    var_05hr['SM2']   = f_05hr.variables['SoilMoist'][:,1,0,0]
    var_05hr['SM3']   = f_05hr.variables['SoilMoist'][:,2,0,0]
    var_05hr['SM4']   = f_05hr.variables['SoilMoist'][:,3,0,0]
    var_05hr['SM5']   = f_05hr.variables['SoilMoist'][:,4,0,0]
    var_05hr['SM6']   = f_05hr.variables['SoilMoist'][:,5,0,0]
    var_05hr['GWMoist'] = f_05hr.variables['GWMoist'][:,0,0]


    var_05hr_cldstrt          = pd.DataFrame(f_05hr_cldstrt.variables['Rainf'][:,0,0]*60.*60.*24., columns=['Rainf'])
    var_05hr_cldstrt['Evap']  = f_05hr_cldstrt.variables['Evap'][:,0,0]*60.*60.*24.
    var_05hr_cldstrt['TVeg']  = f_05hr_cldstrt.variables['TVeg'][:,0,0]*60.*60.*24.
    var_05hr_cldstrt['ESoil'] = f_05hr_cldstrt.variables['ESoil'][:,0,0]*60.*60.*24.
    var_05hr_cldstrt['ECanop']= f_05hr_cldstrt.variables['ECanop'][:,0,0]*60.*60.*24.
    var_05hr_cldstrt['Qs']    = f_05hr_cldstrt.variables['Qs'][:,0,0]*60.*60.*24.
    var_05hr_cldstrt['Qsb']   = f_05hr_cldstrt.variables['Qsb'][:,0,0]*60.*60.*24.
    var_05hr_cldstrt['SM1']   = f_05hr_cldstrt.variables['SoilMoist'][:,0,0,0]
    var_05hr_cldstrt['SM2']   = f_05hr_cldstrt.variables['SoilMoist'][:,1,0,0]
    var_05hr_cldstrt['SM3']   = f_05hr_cldstrt.variables['SoilMoist'][:,2,0,0]
    var_05hr_cldstrt['SM4']   = f_05hr_cldstrt.variables['SoilMoist'][:,3,0,0]
    var_05hr_cldstrt['SM5']   = f_05hr_cldstrt.variables['SoilMoist'][:,4,0,0]
    var_05hr_cldstrt['SM6']   = f_05hr_cldstrt.variables['SoilMoist'][:,5,0,0]
    var_05hr_cldstrt['GWMoist'] = f_05hr_cldstrt.variables['GWMoist'][:,0,0]

    var_3hr          = pd.DataFrame(f_3hr.variables['Rainf'][:,loc_lat,loc_lon]*60.*60.*24., columns=['Rainf'])
    var_3hr['Evap']  = f_3hr.variables['Evap'][:,loc_lat,loc_lon]*60.*60.*24.
    var_3hr['TVeg']  = f_3hr.variables['TVeg'][:,loc_lat,loc_lon]*60.*60.*24.
    var_3hr['ESoil'] = f_3hr.variables['ESoil'][:,loc_lat,loc_lon]*60.*60.*24.
    var_3hr['ECanop']= f_3hr.variables['ECanop'][:,loc_lat,loc_lon]*60.*60.*24.
    var_3hr['Qs']    = f_3hr.variables['Qs'][:,loc_lat,loc_lon]*60.*60.*24.
    var_3hr['Qsb']   = f_3hr.variables['Qsb'][:,loc_lat,loc_lon]*60.*60.*24.
    var_3hr['SM1']   = f_3hr.variables['SoilMoist'][:,0,loc_lat,loc_lon]
    var_3hr['SM2']   = f_3hr.variables['SoilMoist'][:,1,loc_lat,loc_lon]
    var_3hr['SM3']   = f_3hr.variables['SoilMoist'][:,2,loc_lat,loc_lon]
    var_3hr['SM4']   = f_3hr.variables['SoilMoist'][:,3,loc_lat,loc_lon]
    var_3hr['SM5']   = f_3hr.variables['SoilMoist'][:,4,loc_lat,loc_lon]
    var_3hr['SM6']   = f_3hr.variables['SoilMoist'][:,5,loc_lat,loc_lon]
    var_3hr['GWMoist'] = f_3hr.variables['GWMoist'][:,loc_lat,loc_lon]

    f_05hr.close()
    f_05hr_cldstrt.close()
    f_3hr.close()

    var_05hr.to_csv("./csv/water_balance_05hr.csv")
    var_05hr_cldstrt.to_csv("./csv/water_balance_05hr_cldstrt.csv")
    var_3hr.to_csv("./csv/water_balance_3hr.csv")

if __name__ == "__main__":

    fname_05hr        = "/g/data/w35/mm3972/model/cable/runs/my_version/run_Princeton_single_pixel/ctl/outputs/cable_out_2008.nc"
    fname_05hr_cldstrt= "/g/data/w35/mm3972/model/cable/runs/my_version/run_Princeton_single_pixel/ctl/outputs/cable_out_2008_coldstart.nc"
    fname_3hr         = "/g/data/w35/mm3972/model/cable/runs/my_version/run_Princeton/ctl_daily/outputs/cable_out_2008.nc"
    main(fname_05hr, fname_05hr_cldstrt, fname_3hr)
