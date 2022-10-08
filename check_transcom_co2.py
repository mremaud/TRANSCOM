import pandas as pd
import numpy as np
import os
import xarray as xr
from netCDF4 import Dataset
import datetime
from useful import *
import matplotlib.pyplot as plt
import copy

begy=2010
endy=2020

dir_result="/home/users/mremaud/CIF/LAUNCH/TRANSCOM/"
dir_flux="/home/surface1/mremaud/COS/INPUT_sflx/TRANSCOM_FLUX/CO2_Flux_1x1/"

components=os.listdir(dir_result)
components=[ff for ff in components if "TIER1"  in ff]
components=[ff for ff in components if "CO2"  in ff]


nn=1
for dd in components:
 if '._' in dd: continue
 nn+=1

f,ax=plt.subplots(1,1)
row=0
col=0
for dd in ["SIB4","ORC"]:
  gpp_m=xr.open_dataset(dir_result+"CO2_CO2-GPP-"+dd+"-Month_LSCE_LMDz_TIER1_20200207.nc")
  ter_m=xr.open_dataset(dir_result+"CO2_CO2-TER-"+dd+"-Month_LSCE_LMDz_TIER1_20200207.nc")
  ter_d=xr.open_dataset(dir_result+"CO2_CO2-TER-"+dd+"-Diurn_LSCE_LMDz_TIER1_20200207.nc")
  gpp_d=xr.open_dataset(dir_result+"CO2_CO2-GPP-"+dd+"-Diurn_LSCE_LMDz_TIER1_20200207.nc")
  gpp_d.co2.values=gpp_d.co2.values+ter_d.co2.values
  gpp_m.co2.values=gpp_m.co2.values+ter_m.co2.values
  gpp_d=gpp_d.to_dataframe()
  gpp_m=gpp_m.to_dataframe()
  gpp_m=gpp_m[gpp_m.station=="MLO"]
  gpp_d=gpp_d[gpp_d.station=="MLO"]
  gpp_m.set_index("time",inplace=True)
  gpp_d.set_index("time",inplace=True)

  ax.plot(gpp_m.index,gpp_m.co2*10**(6),label=dd+"-Month")
  ax.plot(gpp_d.index,gpp_d.co2*10**(6),label=dd+"-Diurn")

  ax.set_ylabel("CO2 [ppm]") 
  ax.set_title("TER+GPP")
plt.suptitle("CO2")
