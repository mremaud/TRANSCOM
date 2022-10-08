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
dir_flux="/home/surface1/mremaud/COS/INPUT_sflx/TRANSCOM_FLUX/COS_Flux_1x1/"

components=os.listdir(dir_result)
components=[ff for ff in components if "TIER1"  in ff]

components=os.listdir(dir_flux)

nn=1
for dd in components:
 if '._' in dd: continue
 nn+=1

f,ax=plt.subplots(3,2)
row=0
col=0
for dd in components:
 if '._' in dd: continue
 list_flux=os.listdir(dir_flux+dd)
 list_flux=[ff[:-12] for ff in list_flux if "phy"  in ff]
 list_flux=np.unique(list_flux)
 compound=copy.copy(dd[:3])
 for flux in  list_flux:
  if '._' in flux: continue
  #if "Diurn" in flux: continue
  data=xr.open_dataset(dir_result+compound+"_"+flux+"_LSCE_LMDz_TIER1_20200207.nc")
  data=data.to_dataframe()
  data=data[data.station=="MLO"]
  data.set_index("time",inplace=True)
  ax[row,col].plot(data.index,data.cos*10**(12),label=flux)
  ax[row,col].set_ylabel("COS [ppt]") 
  ax[row,col].set_title(flux[4:8])
 ax[row,col].legend(fontsize=8)
 if  row<2: ax[row,col].set_xlabel(" ")
 if  col==1: ax[row,col].set_ylabel(" ")
  
 row+=1
 if row==3: 
   row=0
   col=1
plt.subplots_adjust(wspace=0.6,hspace=0.6)
plt.suptitle("COS")

#######DIURNE EFFECT ON GPP ORC######
plt.figure()
gpp_d=xr.open_dataset(dir_result+"COS_COS-VEG-SIB4-Diurn_LSCE_LMDz_TIER1_20200207.nc")
gpp_d=gpp_d.to_dataframe()
gpp_d=gpp_d[gpp_d.station=="MLO"]
gpp_d=gpp_d[gpp_d.time.dt.year==2016]
gpp_d.set_index("time",inplace=True)
gpp_m=xr.open_dataset(dir_result+"COS_COS-VEG-SIB4-Month_LSCE_LMDz_TIER1_20200207.nc")
gpp_m=gpp_m.to_dataframe()
gpp_m=gpp_m[gpp_m.station=="MLO"]
gpp_m=gpp_m[gpp_m.time.dt.year==2016]
gpp_m.set_index("time",inplace=True)
plt.plot((gpp_d.cos-gpp_m.cos)*10**(12))
plt.ylabel("COS [ppt]")
plt.title("COS-VEG-SIB4-Diurn-COS-VEG-SIB4-Month at MLO")
