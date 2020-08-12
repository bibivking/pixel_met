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

def main(fname_lis,fname_hess,loc_lat_lis,loc_lon_lis):

    lis  = nc.Dataset(fname_lis, "r")
    hess = nc.Dataset(fname_hess, "r")


    para = pd.DataFrame([lis.variables['Soiltype_inst'][0,loc_lat_lis,loc_lon_lis],
                        hess.variables['isoil'][0,0]],
                        columns=['Soiltype'])

    para['sand'] = [lis.variables['SandFrac_inst'][0,loc_lat_lis,loc_lon_lis],
                    hess.variables['sand'][0,0]]
    para['silt'] = [lis.variables['SiltFrac_inst'][0,loc_lat_lis,loc_lon_lis],
                    hess.variables['silt'][0,0]]
    para['clay'] = [lis.variables['ClayFrac_inst'][0,loc_lat_lis,loc_lon_lis],
                    hess.variables['clay'][0,0]]
    para['ssat'] = [lis.variables['Porosity_inst'][0,loc_lat_lis,loc_lon_lis],
                    hess.variables['ssat'][0,0]]
    para['bch']  = [lis.variables['bch_inst'][0,loc_lat_lis,loc_lon_lis],
                    hess.variables['bch'][0,0]]
    para['sfc']  = [lis.variables['sfc_inst'][0,loc_lat_lis,loc_lon_lis],
                    hess.variables['sfc'][0,0]]
    para['swilt']= [lis.variables['swilt_inst'][0,loc_lat_lis,loc_lon_lis],
                    hess.variables['swilt'][0,0]]
    para['hyds'] = [lis.variables['hyds_inst'][0,loc_lat_lis,loc_lon_lis],
                    hess.variables['bch'][0,0]]
    para['sucs'] = [lis.variables['sucs_inst'][0,loc_lat_lis,loc_lon_lis],
                    hess.variables['sucs'][0,0]]
    para['css']  = [lis.variables['css_inst'][0,loc_lat_lis,loc_lon_lis],
                    hess.variables['css'][0,0]]
    para['rhosoil'] = [lis.variables['rhosoil_inst'][0,loc_lat_lis,loc_lon_lis],
                      hess.variables['rhosoil'][0,0]]
    para['elev']    = [lis.variables['Elevation_inst'][0,loc_lat_lis,loc_lon_lis],
                      hess.variables['elev'][0,0]]

    hess.close()
    lis.close()

    para.to_csv("./csv/soil_parameters.csv")


if __name__ == "__main__":


    # lat_-355_lon_1495
    loc_lat = 54  # -35.5
    loc_lon = 149 # 149
    loc_lat_lis = 40
    loc_lon_lis = 140

    fname_hess = "/g/data/w35/mm3972/model/cable/runs/my_version/run_Princeton_single_pixel/ctl/outputs/lat_-355_lon_1495_grass/cable_out_2008.nc"
    fname_lis  = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/Princeton_ctl_para/LIS_output/LIS.CABLE.2008020100.d01.nc"

    main(fname_lis,fname_hess,loc_lat_lis,loc_lon_lis)
