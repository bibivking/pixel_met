#!/usr/bin/env python

"""
Run add gamma to restart file in order to run Hvrd's plant water stress function
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
from scipy.signal import savgol_filter

def main(sm_frm, restart_fname,hist_fname,loc_lat,loc_lon):

    # create file and write global attributes
    f     = nc.Dataset(restart_fname, "r+", format="NETCDF4")
    cable = nc.Dataset(hist_fname, 'r')


    if sm_frm == "HESS":
        print(cable.variables['y'][loc_lat])
        print(cable.variables['x'][loc_lon])
        f.variables['wb'][:,0]   = cable.variables['SoilMoist'][364,:,loc_lat,loc_lon].filled(-9999.)
        f.variables['GWwb'][0] = cable.variables['GWMoist'][364,loc_lat,loc_lon].filled(-9999.)
        print(cable.variables['SoilMoist'][364,:,loc_lat,loc_lon])
        print(cable.variables['GWMoist'][364,loc_lat,loc_lon])
        print(f.variables['wb'][:,0])
        print(f.variables['GWwb'][0])
    if sm_frm == "LIS":
        print(f.variables['wb'][:,0])
        print(f.variables['GWwb'][0])
        f.variables['wb'][:,0] = cable.variables['SoilMoist_tavg'][30,:,loc_lat,loc_lon].filled(-9999.)
        f.variables['GWwb'][0] = cable.variables['GWwb_tavg'][30,loc_lat,loc_lon].filled(-9999.)
        print(cable.variables['SoilMoist_tavg'][30,:,loc_lat,loc_lon])
        print(cable.variables['GWwb_tavg'][30,loc_lat,loc_lon])
        print(f.variables['wb'][:,0])
        print(f.variables['GWwb'][0])


    f.close()
    cable.close()

if __name__ == "__main__":

    sm_frm = "LIS"
    case_name = ''
    restart_fname = "/g/data/w35/mm3972/model/cable/runs/my_version/run_Princeton_single_pixel/ctl/outputs/restart_2007.nc"
    #restart_fname = "/g/data/w35/mm3972/model/cable/runs/my_version/run_Princeton_single_pixel/ctl/outputs/%s/restart_2007.nc" % case_name

    if sm_frm == "HESS":
        hist_fname    = "/g/data/w35/mm3972/model/cable/runs/my_version/run_Princeton/ctl_daily/outputs/cable_out_2007.nc"
        # # lat_-355_lon_1495
        # loc_lat = 54  # -35.5
        # loc_lon = 149 # 149
        # lat = -33.5 lon = 150.5 PFT: 2-evergreen broadleaf
        loc_lat = 56
        loc_lon = 150
    if sm_frm == "LIS":
        hist_fname    = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/Princeton_ctl/LIS_output/LIS.CABLE.2007120100.d01.nc"
        # # lat_-355_lon_1495
        # loc_lat = 40  # -35.5
        # loc_lon = 140 # 149
        # lat = -33.604504 lon = 150.60345 PFT: 2-evergreen broadleaf
        loc_lat = 47
        loc_lon = 144
    main(sm_frm, restart_fname,hist_fname,loc_lat,loc_lon)
