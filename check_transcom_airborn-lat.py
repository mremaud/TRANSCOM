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
remove_bias=1
models=["LMDZ_LSCE","TM5_UU","TM3_MPI","MIROC4_JAMSTEC","NICAM5_TM","NICAM6_TM"]
nmodels=["LMDZ","TM5","TM3","MIROC","NICAM5","NICAM6"]
colors=["red","blue","orange","violet","peru","saddlebrown"]
lines=["solid","dotted","dashdot","dashed","dashed"]

process={"Name":["Opt"],"process":["OPT"]}
process=pd.DataFrame(process)

latref=np.arange(-90,90,10)
latref2=np.arange(-90,90,30)
y2=np.arange(-10,70,20)

dir_result="/home/surface1/mremaud/COS/TRANSCOM/OUTPUTS/"

####FIND THE FIXED STATION
list_stations=[]
obs=pd.read_csv("/home/surface1/mremaud/COS/TRANSCOM/TRANSCOM_OBS/DATA_TIER2/TIER2.csv")
obs["date (UTC)"]=pd.to_datetime(obs["date (UTC)"])
obs=obs[obs["date (UTC)"].dt.year>=2010]
obs=obs.groupby(["lat","flight"]).mean().reset_index()
obs=obs[obs.lat>0]
obs=obs.sort_values("lat")

list_stations=["ATOM","HIPPO"]

f,ax=plt.subplots(2,1)
row=0 ; col=0

which=["JM","LSCE"]

for iss,ss in enumerate(list_stations):
 ##Open observations###############################################################
 obs_stat=pd.read_pickle("/home/satellites16/mremaud/OBS/"+ss.lower()+".pkl")
 if ss== "ATOM":
  obs_stat["lat"]=obs_stat.lat.round(1)
  obs_stat["lon"]=obs_stat.lon.round(1)
  obs_stat["alt"]=obs_stat.lon.round(1)

 print(obs_stat)
 obs_stat=obs_stat.rename(columns={"COS":"cos","date":"time"})
 obs_stat2=obs_stat.sort_values(by='time')
 obs_stat["lat"]=obs_stat2.apply(lambda row: latref[np.abs(row.lat-latref).argmin()],axis=1)
 obs_stat=obs_stat.groupby(["time","lat"]).mean().reset_index()
 obs_stat=obs_stat.set_index("time").groupby("lat").resample("M").mean().dropna()
 del obs_stat["lat"];obs_stat.reset_index(inplace=True)
 obs_stat=obs_stat.groupby(["lat"]).mean().reset_index()
 mean_obs=obs_stat.cos.mean()
 ax[row].plot(obs_stat.lat,obs_stat.cos,label="obs",color="k")
 ##################################################################################

 for ww,wh in enumerate(which):
  for im,mm in enumerate(models):
   name_file=dir_result+"/"+mm+"/TIER2/COS_COS-OPT-"+wh+"-Month_"+mm+"_TIER2.nc"
   if not os.path.exists(name_file): continue
   data=xr.open_dataset(name_file)
   data=data.to_dataframe()
   data["station"]=data.apply(lambda row: "HIPPO" if "HIPPO" in row.station else row.station,axis=1)
   data=data[data.station==ss]
   if data.empty: continue
   data=data.sort_values(by='time')
   data["lat"]=data.apply(lambda row: latref[np.abs(row.lat-latref).argmin()],axis=1)
   data=data.groupby(["time","lat"]).mean().reset_index()
   data=data.set_index("time").groupby("lat").resample("M").mean()
   del data["lat"];data.reset_index(inplace=True)
   data=data.groupby(["lat"]).mean()
   if remove_bias:
    bias=data.cos.mean()-mean_obs
    print(ss,"obs",mean_obs,mm,wh,data.cos.mean())
    bias=data.cos.mean()-mean_obs
   ax[row].plot(data.index,data.cos-bias,label=mm+" "+wh,color=colors[im],linestyle=lines[ww])
 ax[row].set_xlabel("Latitude") 
 #ax[row].set_ylim(-15-bias,30-bias)
 ax[row].set_xlim(-90,90)

 #ax[row].set_yticks(y2)
 ax[row].tick_params(labelsize=10,axis="x") 
 ax[row].set_title(ss,fontsize=10)
 if (row==1): ax[row].legend(fontsize=8,bbox_to_anchor=(2.5,1.5), loc="upper right")
 ax[row].set_ylabel("COS [ppt]")

 if  row==0: ax[row].set_xlabel(" ")
 # plt.show()  
 row+=1

plt.subplots_adjust(wspace=0.6,hspace=0.6)
plt.suptitle("Optimized fluxes")


