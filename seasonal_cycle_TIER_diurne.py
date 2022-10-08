import pandas as pd
import numpy as np
import os
import xarray as xr
from netCDF4 import Dataset
import datetime
from useful import *
import matplotlib.pyplot as plt
import copy
import collections

begy=2010
endy=2020
models=["TOMCAT_ULeic","LMDZ_LSCE","TM5_UU","TM3_MPI","MIROC4_JAMSTEC","NICAM5_TM","NICAM6_TM"]
#models=["LMDZ_LSCE"]
#,"MIROC4_JAMSTEC","TM3_MPI","NICAM5_TM","NICAM6_TM","TOMCAT_ULeic"]
nmodels=["TOMCAT","LMDZ","TM5","TM3","MIROC4","NICAM5","NICAM6"]
colors=["darkgreen","red","blue","brown","violet","orange","saddlebrown"]
lines=["solid","dotted","dashdot","dashed","dashed"]

dir_result="/home/surface1/mremaud/COS/TRANSCOM/OUTPUTS/"
dir_flux="/home/surface1/mremaud/COS/INPUT_sflx/TRANSCOM_FLUX/COS_Flux_1x1/"


#############################DIURNAL###################################
compound="COS"
flux=["SIB4","ORC"]
colors=["k","k"]
list_station=["SPO","PSA","CGO","SMO","KUM","MLO","WIS","NWR","THD","HFM","LEF","GIF","MHD","BRW","SUM","ALT"]
LEGEND_HEIGHT = 0.1
#f,ax=plt.subplots(4,4,figsize=(5,4))
for iff,ff in enumerate(flux):
 f,ax=plt.subplots(4,4,figsize=(8,6))
 row=0; col=0
 for stat in list_station:
   CycleS=pd.DataFrame()
   for im,mm in enumerate(models):
     if not os.path.exists(dir_result+"/"+mm+"/TIER1/"+compound+"_COS-VEG-"+flux[iff]+"-Diurn_"+mm+"_TIER1.nc"): continue
     for itt,tt in enumerate(["SOIL","VEG"]):
#     for itt,tt in enumerate(["VEG"]):

      data_D2=xr.open_dataset(dir_result+"/"+mm+"/TIER1/"+compound+"_COS-"+tt+"-"+flux[iff]+"-Diurn_"+mm+"_TIER1.nc").to_dataframe()
      data_M2=xr.open_dataset(dir_result+"/"+mm+"/TIER1/"+compound+"_COS-"+tt+"-"+flux[iff]+"-Month_"+mm+"_TIER1.nc").to_dataframe()
      data_M2=data_M2[(data_M2.station==stat)&(data_M2.time.dt.year==2016)].copy(deep=True).set_index("time")
      data_D2=data_D2[(data_D2.station==stat)&(data_D2.time.dt.year==2016)].copy(deep=True).set_index("time")
      if itt==0:
       data_D=data_D2.copy(deep=True); data_M=data_M2.copy(deep=True)
      else:
       data_D.cos+=data_D2.cos; data_M.cos+=data_M2.cos
     if data_D.empty: continue
     data_D=data_D.resample("M").mean().dropna()
     data_M=data_M.resample("M").mean().dropna()
     data_D["cos"]=data_D["cos"]-data_M["cos"]
     data_D["model"]=nmodels[im]
     CycleS=CycleS.append(data_D.reset_index())
   CycleS2=CycleS.groupby(["time"]).mean()
   CycleS2["std"]=CycleS.groupby(["time"]).std().cos
   ax[row,col].plot(np.arange(1,13),CycleS2.cos,color=colors[iff])
   ax[row,col].fill_between(np.arange(1,13),CycleS2.cos-CycleS2["std"],CycleS2.cos+CycleS2["std"],color=colors[iff],alpha=0.5)
   ax[row,col].set_xticks(np.arange(1,13,2))
   ax[row,col].set_yticks(np.arange(-15,20,5))
   ax[row,col].set_yticklabels(['-15',' ',' ','0',' ',' ','15'])
   ax[row,col].set_xticklabels(['Jan','Mar','May','Jul','Sep','Nov'])
   ax[row,col].set_xlim(1,12)
   ax[row,col].tick_params(rotation=70,labelsize=12,axis="x",direction='in')
   ax[row,col].tick_params(labelsize=12,axis="y",direction='in')

   if col==0 : ax[row,col].set_ylabel("COS [ppt]",fontsize=13)
   if col != 0: ax[row,col].set_yticks([])
   if row !=3: ax[row,col].set_xticks([])
   ax[row,col].set_title(stat,pad=0.03,fontsize=12)
   ax[row,col].set_ylim(-19,19)
   col+=1
   if col==4: col=0; row+=1
 f.subplots_adjust(wspace=0.2, hspace=0.3)
 plt.savefig('/home/users/mremaud/PYTHON/COS/TRANSCOM/FIG/Diurnal'+ff+'.png',format='png', bbox_inches='tight', dpi = 500)

