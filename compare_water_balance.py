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

def main( fname, fname_cldstrt, fname_3hr_hess, fname_3hr_LIS,
          loc_lat, loc_lon, loc_lat_lis, loc_lon_lis, case_name, timestep):

    f         = nc.Dataset(fname, "r")
    # f_cldstrt = nc.Dataset(fname_cldstrt, "r")
    f_3hr_hess= nc.Dataset(fname_3hr_hess, 'r')

    print(f_3hr_hess.variables['x'][loc_lon])
    print(f_3hr_hess.variables['y'][loc_lat])
    print(f.variables['x'][0])
    print(f.variables['y'][0])

    var          = pd.DataFrame(f.variables['Rainf'][:,0,0]*60.*60.*24., columns=['Rainf'])
    var['Evap']  = f.variables['Evap'][:,0,0]*60.*60.*24.
    var['TVeg']  = f.variables['TVeg'][:,0,0]*60.*60.*24.
    var['ESoil'] = f.variables['ESoil'][:,0,0]*60.*60.*24.
    var['ECanop']= f.variables['ECanop'][:,0,0]*60.*60.*24.
    var['Qs']    = f.variables['Qs'][:,0,0]*60.*60.*24.
    var['Qsb']   = f.variables['Qsb'][:,0,0]*60.*60.*24.
    var['SM1']   = f.variables['SoilMoist'][:,0,0,0]
    var['SM2']   = f.variables['SoilMoist'][:,1,0,0]
    var['SM3']   = f.variables['SoilMoist'][:,2,0,0]
    var['SM4']   = f.variables['SoilMoist'][:,3,0,0]
    var['SM5']   = f.variables['SoilMoist'][:,4,0,0]
    var['SM6']   = f.variables['SoilMoist'][:,5,0,0]
    var['GWMoist'] = f.variables['GWMoist'][:,0,0]


    # var_cldstrt          = pd.DataFrame(f_cldstrt.variables['Rainf'][:,0,0]*60.*60.*24., columns=['Rainf'])
    # var_cldstrt['Evap']  = f_cldstrt.variables['Evap'][:,0,0]*60.*60.*24.
    # var_cldstrt['TVeg']  = f_cldstrt.variables['TVeg'][:,0,0]*60.*60.*24.
    # var_cldstrt['ESoil'] = f_cldstrt.variables['ESoil'][:,0,0]*60.*60.*24.
    # var_cldstrt['ECanop']= f_cldstrt.variables['ECanop'][:,0,0]*60.*60.*24.
    # var_cldstrt['Qs']    = f_cldstrt.variables['Qs'][:,0,0]*60.*60.*24.
    # var_cldstrt['Qsb']   = f_cldstrt.variables['Qsb'][:,0,0]*60.*60.*24.
    # var_cldstrt['SM1']   = f_cldstrt.variables['SoilMoist'][:,0,0,0]
    # var_cldstrt['SM2']   = f_cldstrt.variables['SoilMoist'][:,1,0,0]
    # var_cldstrt['SM3']   = f_cldstrt.variables['SoilMoist'][:,2,0,0]
    # var_cldstrt['SM4']   = f_cldstrt.variables['SoilMoist'][:,3,0,0]
    # var_cldstrt['SM5']   = f_cldstrt.variables['SoilMoist'][:,4,0,0]
    # var_cldstrt['SM6']   = f_cldstrt.variables['SoilMoist'][:,5,0,0]
    # var_cldstrt['GWMoist'] = f_cldstrt.variables['GWMoist'][:,0,0]

    var_3hr_hess          = pd.DataFrame(f_3hr_hess.variables['Rainf'][:,loc_lat,loc_lon]*60.*60.*24., columns=['Rainf'])
    var_3hr_hess['Evap']  = f_3hr_hess.variables['Evap'][:,loc_lat,loc_lon]*60.*60.*24.
    var_3hr_hess['TVeg']  = f_3hr_hess.variables['TVeg'][:,loc_lat,loc_lon]*60.*60.*24.
    var_3hr_hess['ESoil'] = f_3hr_hess.variables['ESoil'][:,loc_lat,loc_lon]*60.*60.*24.
    var_3hr_hess['ECanop']= f_3hr_hess.variables['ECanop'][:,loc_lat,loc_lon]*60.*60.*24.
    var_3hr_hess['Qs']    = f_3hr_hess.variables['Qs'][:,loc_lat,loc_lon]*60.*60.*24.
    var_3hr_hess['Qsb']   = f_3hr_hess.variables['Qsb'][:,loc_lat,loc_lon]*60.*60.*24.
    var_3hr_hess['SM1']   = f_3hr_hess.variables['SoilMoist'][:,0,loc_lat,loc_lon]
    var_3hr_hess['SM2']   = f_3hr_hess.variables['SoilMoist'][:,1,loc_lat,loc_lon]
    var_3hr_hess['SM3']   = f_3hr_hess.variables['SoilMoist'][:,2,loc_lat,loc_lon]
    var_3hr_hess['SM4']   = f_3hr_hess.variables['SoilMoist'][:,3,loc_lat,loc_lon]
    var_3hr_hess['SM5']   = f_3hr_hess.variables['SoilMoist'][:,4,loc_lat,loc_lon]
    var_3hr_hess['SM6']   = f_3hr_hess.variables['SoilMoist'][:,5,loc_lat,loc_lon]
    var_3hr_hess['GWMoist'] = f_3hr_hess.variables['GWMoist'][:,loc_lat,loc_lon]

    f.close()
    # f_cldstrt.close()
    f_3hr_hess.close()

    var_3hr_LIS = get_daily_LIS_output(fname_3hr_LIS,loc_lat_lis,loc_lon_lis)


    var.to_csv("./csv/water_balance_%s_%s_%s_%s.csv" % (timestep, loc_lat, loc_lon,case_name))
    # var_cldstrt.to_csv("./csv/water_balance_%s_cldstrt_%s_%s_%s.csv" % (timestep, loc_lat, loc_lon,case_name))
    var_3hr_hess.to_csv("./csv/water_balance_3hr_hess_%s_%s.csv" % (loc_lat, loc_lon))
    var_3hr_LIS.to_csv("./csv/water_balance_3hr_LIS_%s_%s.csv" % (loc_lat, loc_lon))



def get_daily_LIS_output(input_fname,loc_lat_lis,loc_lon_lis):

    """
    read met fields from LIS-CABLE output
    """

    print("carry on read_cable_var")

    for month in np.arange(0,12,1):
        print(month)
        cable = nc.Dataset(input_fname[month], 'r')

        if month == 0:
            rain  = cable.variables['Rainf_f_tavg'][:,loc_lat_lis,loc_lon_lis].filled(-9999.)
            evap  = cable.variables['Evap_tavg'][:,loc_lat_lis,loc_lon_lis].filled(-9999.)
            tveg  = cable.variables['TVeg_tavg'][:,loc_lat_lis,loc_lon_lis].filled(-9999.)
            esoil = cable.variables['ESoil_tavg'][:,loc_lat_lis,loc_lon_lis].filled(-9999.)
            ecanop= cable.variables['ECanop_tavg'][:,loc_lat_lis,loc_lon_lis].filled(-9999.)
            qs    = cable.variables['Qs_tavg'][:,loc_lat_lis,loc_lon_lis].filled(-9999.)
            qsb   = cable.variables['Qsb_tavg'][:,loc_lat_lis,loc_lon_lis].filled(-9999.)
            sm1   = cable.variables['SoilMoist_tavg'][:,0,loc_lat_lis,loc_lon_lis].filled(-9999.)
            sm2   = cable.variables['SoilMoist_tavg'][:,1,loc_lat_lis,loc_lon_lis].filled(-9999.)
            sm3   = cable.variables['SoilMoist_tavg'][:,2,loc_lat_lis,loc_lon_lis].filled(-9999.)
            sm4   = cable.variables['SoilMoist_tavg'][:,3,loc_lat_lis,loc_lon_lis].filled(-9999.)
            sm5   = cable.variables['SoilMoist_tavg'][:,4,loc_lat_lis,loc_lon_lis].filled(-9999.)
            sm6   = cable.variables['SoilMoist_tavg'][:,5,loc_lat_lis,loc_lon_lis].filled(-9999.)
            gwwb  = cable.variables['GWwb_tavg'][:,loc_lat_lis,loc_lon_lis].filled(-9999.)

        else:
            rain  = np.concatenate((rain,cable.variables['Rainf_f_tavg'][:,loc_lat_lis,loc_lon_lis].filled(-9999.)))
            evap  = np.concatenate((evap,cable.variables['Evap_tavg'][:,loc_lat_lis,loc_lon_lis].filled(-9999.)))
            tveg  = np.concatenate((tveg,cable.variables['TVeg_tavg'][:,loc_lat_lis,loc_lon_lis].filled(-9999.)))
            esoil = np.concatenate((esoil,cable.variables['ESoil_tavg'][:,loc_lat_lis,loc_lon_lis].filled(-9999.)))
            ecanop= np.concatenate((ecanop,cable.variables['ECanop_tavg'][:,loc_lat_lis,loc_lon_lis].filled(-9999.)))
            qs    = np.concatenate((qs,cable.variables['Qs_tavg'][:,loc_lat_lis,loc_lon_lis].filled(-9999.)))
            qsb   = np.concatenate((qsb,cable.variables['Qsb_tavg'][:,loc_lat_lis,loc_lon_lis].filled(-9999.)))
            sm1   = np.concatenate((sm1,cable.variables['SoilMoist_tavg'][:,0,loc_lat_lis,loc_lon_lis].filled(-9999.)))
            sm2   = np.concatenate((sm2,cable.variables['SoilMoist_tavg'][:,1,loc_lat_lis,loc_lon_lis].filled(-9999.)))
            sm3   = np.concatenate((sm3,cable.variables['SoilMoist_tavg'][:,2,loc_lat_lis,loc_lon_lis].filled(-9999.)))
            sm4   = np.concatenate((sm4,cable.variables['SoilMoist_tavg'][:,3,loc_lat_lis,loc_lon_lis].filled(-9999.)))
            sm5   = np.concatenate((sm5,cable.variables['SoilMoist_tavg'][:,4,loc_lat_lis,loc_lon_lis].filled(-9999.)))
            sm6   = np.concatenate((sm6,cable.variables['SoilMoist_tavg'][:,5,loc_lat_lis,loc_lon_lis].filled(-9999.)))
            gwwb  = np.concatenate((gwwb,cable.variables['GWwb_tavg'][:,loc_lat_lis,loc_lon_lis].filled(-9999.)))

        cable.close()

    Var          = pd.DataFrame(rain*60.*60.*24., columns=['Rainf'])
    Var['Evap']  = evap*60.*60.*24.
    Var['TVeg']  = tveg*60.*60.*24.
    Var['ESoil'] = esoil*60.*60.*24.
    Var['ECanop']= ecanop*60.*60.*24.
    Var['Qs']    = qs*60.*60.*24.
    Var['Qsb']   = qsb*60.*60.*24.
    Var['SM1']   = sm1
    Var['SM2']   = sm2
    Var['SM3']   = sm3
    Var['SM4']   = sm4
    Var['SM5']   = sm5
    Var['SM6']   = sm6
    Var['GWMoist'] = gwwb

    return Var


if __name__ == "__main__":


    # # lat_-355_lon_1495
    # loc_lat = 54  # -35.5
    # loc_lon = 149 # 149
    # loc_lat_lis = 40
    # loc_lon_lis = 140
    # case_name = "lat_-355_lon_1495_grass_LIS_root"

    # lat = -33.5 lon = 150.5
    loc_lat = 56
    loc_lon = 150
    loc_lat_lis = 47
    loc_lon_lis = 144
    case_name = "lat_-335_lon_1505_broadleaf_forest"
    # case_name = "lat_-335_lon_1505_forest_LIS_all"
    #"lat_-355_lon_1495_grass_LIS_g1_vcmax"

    timestep = "05hr"
    fname        = "/g/data/w35/mm3972/model/cable/runs/my_version/run_Princeton_single_pixel/ctl/outputs/%s/cable_out_2008.nc" % case_name
    #fname        = "/g/data/w35/mm3972/model/cable/runs/my_version/run_Princeton_single_pixel/ctl/outputs/cable_out_2008.nc"
    fname_cldstrt= "/g/data/w35/mm3972/model/cable/runs/my_version/run_Princeton_single_pixel/ctl/outputs/%s/cable_out_2008_cldstrt.nc" % case_name
    fname_3hr_hess = "/g/data/w35/mm3972/model/cable/runs/my_version/run_Princeton/ctl_daily/outputs/cable_out_2008.nc"


    path_LIS = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/Princeton_ctl/LIS_output/"
    fname_3hr_LIS = [path_LIS+"LIS.CABLE.2008010100.d01.nc",path_LIS+"LIS.CABLE.2008020100.d01.nc",
                     path_LIS+"LIS.CABLE.2008030100.d01.nc",path_LIS+"LIS.CABLE.2008040100.d01.nc",
                     path_LIS+"LIS.CABLE.2008050100.d01.nc",path_LIS+"LIS.CABLE.2008060100.d01.nc",
                     path_LIS+"LIS.CABLE.2008070100.d01.nc",path_LIS+"LIS.CABLE.2008080100.d01.nc",
                     path_LIS+"LIS.CABLE.2008090100.d01.nc",path_LIS+"LIS.CABLE.2008100100.d01.nc",
                     path_LIS+"LIS.CABLE.2008110100.d01.nc",path_LIS+"LIS.CABLE.2008120100.d01.nc",
                    ]


    main(fname, fname_cldstrt, fname_3hr_hess, fname_3hr_LIS,
         loc_lat, loc_lon, loc_lat_lis, loc_lon_lis, case_name, timestep)
