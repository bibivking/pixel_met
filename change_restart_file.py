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

def main(restart_fname,hist_fname,loc_lat,loc_lon):


    # create file and write global attributes
    f     = nc.Dataset(restart_fname, "r+", format="NETCDF4")
    cable = nc.Dataset(hist_fname, 'r')

    print(cable.variables['y'][loc_lat])
    print(cable.variables['x'][loc_lon])
    f.variables['wb'][:,0]   = cable.variables['SoilMoist'][364,:,loc_lat,loc_lon].filled(-9999.)
    f.variables['GWwb'][0] = cable.variables['GWMoist'][364,loc_lat,loc_lon].filled(-9999.)
    print(cable.variables['SoilMoist'][364,:,loc_lat,loc_lon])
    print(cable.variables['GWMoist'][364,loc_lat,loc_lon])
    print(f.variables['wb'][:,0])
    print(f.variables['GWwb'][0])
    f.close()
    cable.close()

if __name__ == "__main__":

    # lat_-355_lon_1495
    loc_lat = 54  # -35.5
    loc_lon = 149 # 149

    # # lat = -33.5 lon = 150.5 PFT: 2-evergreen broadleaf
    # loc_lat = 56
    # loc_lon = 150


    restart_fname = "/g/data/w35/mm3972/model/cable/runs/my_version/run_Princeton_single_pixel/ctl/outputs/restart_2007.nc"
    hist_fname    = "/g/data/w35/mm3972/model/cable/runs/my_version/run_Princeton/ctl_daily/outputs/cable_out_2007.nc"
    main(restart_fname,hist_fname,loc_lat,loc_lon)
