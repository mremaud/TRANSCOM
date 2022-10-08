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
import seaborn as sn

####Mean obs####
list_station=['MLO','ALT',"KUM",'BRW','SUM','MHD','LEF','MLO','CGO','HFM','THD','GIF','NWR','SMO','PSA','SPO']
mean_obs=0
for station in list_station:
 data=pd.read_pickle("/home/surface1/mremaud/COS/SURFACE/STATION/"+station+".pkl")
 data=data[(data.date.dt.year==2010)&(data.date.dt.month==1)]
 if data.empty: print(station);continue
 mean_obs+=data.obs.mean()
mean_obs=mean_obs/len(list_station)*10**6
print(mean_obs)
begy=2010
endy=2018

models=["OBS","TOMCAT_ULeic","LMDZ_LSCE","TM5_UU","TM3_MPI","MIROC4_JAMSTEC","NICAM5_TM","NICAM6_TM"]
nmodels=["OBS","TOMCAT","LMDZ","TM5","TM3","MIROC4","NICAM5","NICAM6"]
colors=["k","darkgreen","red","blue","skyblue","violet","orange","saddlebrown"]
dir_result="/home/surface1/mremaud/COS/TRANSCOM/OUTPUTS/"
compound="COS"

ts=pd.DataFrame()

list_flux=["BB-Stin","SOIL-SIB4",'VEG-SIB4','OCEAN-Lennartz','ANTHR-Zumkehr']
#list_flux=["OPT-JM"]
list_station=['BRW'] #,'MLO','ALT',"KUM",'BRW','SUM','MHD','LEF','MLO','CGO','HFM','THD','GIF','NWR','SMO','PSA','SPO']
f,ax=plt.subplots(1,1,figsize=(8,2))
row=0
for station in list_station:
 for im,mm in enumerate(models):
  if mm=="OBS":
    dataf=pd.read_pickle("/home/surface1/mremaud/COS/SURFACE/STATION/"+station+".pkl")
    dataf = dataf.rename({'date': 'time', 'obs': 'cos'}, axis='columns')
    dataf.cos*=10**6
    dataf=dataf[(dataf.time.dt.year>=begy)&(dataf.time.dt.year<=endy)]
    #mean_obs=dataf.iloc[0].cos
    ax.plot(dataf.time,dataf.cos,color=colors[im],label=nmodels[im])
    TS=dataf.time.values
  dataf=pd.DataFrame()
  for iff,flux in enumerate(list_flux):
   if mm =="OBS": continue
   if not os.path.exists(dir_result+"/"+mm+"/TIER1/"+compound+"_"+compound+"-"+flux+"-Month_"+mm+"_TIER1.nc"): continue
   data=xr.open_dataset(dir_result+"/"+mm+"/TIER1/"+compound+"_"+compound+"-"+flux+"-Month_"+mm+"_TIER1.nc")
   data=data.to_dataframe()
   data=data[data.station==station]
   data=data.sort_values(by='time')
   data.set_index("time",inplace=True)
   if iff==0:
    dataf=data.copy(deep=True)
   else:
    dataf.cos+=data.cos
  dataf.reset_index(inplace=True)
  print(mm,len(dataf))
  if dataf.empty: continue
  dataf=dataf[(dataf.time.dt.year>=begy)&(dataf.time.dt.year<=endy)]
  dataf["time"]=TS
  dataf.cos=dataf.cos+mean_obs
  dataf["model"]=mm
  ts=ts.append(dataf)

ts2=ts.groupby(["time"]).mean().reset_index()
ts2["std"]=ts.groupby(["time"]).std().reset_index().cos
ax.plot(ts2.time,ts2.cos,color="orange",label="ref")
ax.fill_between(ts2.time,ts2.cos-ts2["std"],ts2.cos+ts2["std"],alpha=0.1,color="orange")
ax.tick_params(direction='in', length=3, color='k',grid_color="grey",colors="grey",
                grid_alpha=0.5,labelsize=13)
ax.set_ylabel("COS [ppt]",color="grey",fontsize=15)

plt.legend()
plt.savefig("FIG/time_series_ref.png",format='png', bbox_inches='tight', dpi = 400)
