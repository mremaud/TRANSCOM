"""
Marine Remaud, february 2022
Aim: Make plots of the mean vertical gradients averaged over all NOAA airborne stations
"""

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
import seaborn as sns
begy=2011
endy=2018
which= ["COS-OPT-JM-Month"]  #["COS-ANTHR-Zumkehr-Month","COS-BB-Stin-Month","COS-SOIL-SIB4-Month","COS-VEG-SIB4-Month","COS-OCEAN-LennartzPICSES-Month"]  #["COS-OPT-JM-Month"]
name_plot="opt"
models=["LMDZ_LSCE","TM5_UU","TM3_MPI","MIROC4_JAMSTEC","NICAM6_TM","NICAM5_TM"]

nmodels=["LMDZ","TM5","TM3","MIROC4","NICAM6","NICAM5"]
colors=["red","skyblue","blue","orange","violet","peru","saddlebrown"]
lines=["solid","dotted","dashdot","dashed","dashed"]
process={"Name":["Opt"],"process":["OPT"]}
process=pd.DataFrame(process)

zref=np.arange(0,12000,1000)
dir_result="/home/surface1/mremaud/COS/TRANSCOM/OUTPUTS/"
dir_obs_transcom="/home/surface1/mremaud/COS/TRANSCOM/TRANSCOM_OBS/DATA_TIER2/TIER2.csv"
dir_obs="/home/satellites16/mremaud/OBS/"

####FIND THE FIXED STATION
list_stations=[]
obs=pd.read_csv(dir_obs_transcom)
obs["date (UTC)"]=pd.to_datetime(obs["date (UTC)"])
obs=obs[obs["date (UTC)"].dt.year>=2010]
obs=obs.groupby(["lat","flight"]).mean().reset_index()
obs=obs[obs.lat>0]
obs=obs.sort_values("lat")
for stat in obs.flight.unique():
 if not np.abs(obs[obs.flight==stat].lat.min()-obs[obs.flight==stat].lat.max())>1:
  list_stations.append(stat)
print(list_stations)

GradV={"model":[],"cos":[],"station":[],"season":[],"lat":[]}
for iss,ss in enumerate(list_stations):
 ##Open observations###############################################################
 obs_stat=pd.read_pickle(dir_obs+ss.lower()+".pkl")
 print(obs_stat)
 obs_stat.rename(columns={"date":"time","COS":"cos"},inplace=True)
 obs_stat["zref"]=obs_stat.apply(lambda row: zref[np.abs(row.alt-zref).argmin()],axis=1)
 obs_stat.zref/=1000
 obs_stat=obs_stat.sort_values(by='time')
 obs_stat=obs_stat[(obs_stat.time.dt.year>=begy)&(obs_stat.lon>-900)&(obs_stat.time.dt.year<=endy)]
 obs_stat=obs_stat.groupby(["time","zref"]).mean().reset_index()
 obs_stat=obs_stat.set_index("time").groupby("zref").resample("M").mean()
 del obs_stat["zref"];obs_stat.reset_index(inplace=True)
 obs_stat["season"]="DJF"
 obs_stat.loc[(obs_stat.time.dt.month==1)|(obs_stat.time.dt.month==2)|(obs_stat.time.dt.month==12),"season"]="DJF"
 obs_stat.loc[(obs_stat.time.dt.month==3)|(obs_stat.time.dt.month==4)|(obs_stat.time.dt.month==5),"season"]="MAM"
 obs_stat.loc[(obs_stat.time.dt.month==6)|(obs_stat.time.dt.month==7)|(obs_stat.time.dt.month==8),"season"]="JJA"
 obs_stat.loc[(obs_stat.time.dt.month==9)|(obs_stat.time.dt.month==10)|(obs_stat.time.dt.month==11),"season"]="SON"
 obs_stat=obs_stat.groupby(["zref","season"]).mean().reset_index()
 lat_stat=obs_stat.lat.values[0]
 for sea in obs_stat.season.unique():
  obs_sea=obs_stat[(obs_stat.season==sea)].copy(deep=True)
  if obs_sea.empty: continue
  GradV["station"].append(ss)
  GradV["model"].append("obs")
  GradV["season"].append(sea)
  GradV["lat"].append(lat_stat)
  if not obs_stat[obs_stat.zref==1].empty:
   gradv=-obs_sea[obs_sea.zref==4].cos.values[0]+obs_sea[(obs_sea.zref==1)].cos.values[0]
  else:
   gradv=-obs_stat[(obs_stat.zref==4)&(obs_stat.season==sea)].cos.values[0]+obs_stat[(obs_stat.zref==2)&(obs_stat.season==sea)].cos.values[0]
  GradV["cos"].append(gradv)
 ##################################################################################

 for im,mm in enumerate(models):
  for ww,wh in enumerate(which):
   name_file=dir_result+"/"+mm+"/TIER2/COS_"+wh+"_"+mm+"_TIER2.nc"
   print(name_file)
   if not os.path.exists(name_file): continue
   data=xr.open_dataset(name_file,decode_times=True)
   print(name_file,data) 
   data=data.to_dataframe()
   print(data.station.unique())
   data=data[(data.time.dt.year>=begy)&(data.time.dt.year<=endy)]
   print(data.station.unique())
   data=data[data.lon>-999]
   data=data[data.station==ss]
   print(data) 
   if data.empty: continue
   data["zref"]=data.apply(lambda row: zref[np.abs(row.alt-zref).argmin()],axis=1)
   data.zref/=1000
   data=data.sort_values(by='time')
   print(data.time,data)
   data=data.groupby(["time","zref"]).mean().reset_index()
   if ww==0:
    tmp=data.copy(deep=True)
   else:
    print(wh,mm)
    tmp.cos=tmp.cos.values+data.cos.values
  tmp["season"]="DJF"
  tmp=tmp.set_index("time").groupby("zref").resample("M").mean()
  del tmp["zref"];tmp.reset_index(inplace=True)
  tmp.loc[(tmp.time.dt.month==1)|(tmp.time.dt.month==2)|(tmp.time.dt.month==12),"season"]="DJF"
  tmp.loc[(tmp.time.dt.month==3)|(tmp.time.dt.month==4)|(tmp.time.dt.month==5),"season"]="MAM"
  tmp.loc[(tmp.time.dt.month==6)|(tmp.time.dt.month==7)|(tmp.time.dt.month==8),"season"]="JJA"
  tmp.loc[(tmp.time.dt.month==9)|(tmp.time.dt.month==10)|(tmp.time.dt.month==11),"season"]="SON"
  tmp=tmp.groupby(["zref","season"]).mean().reset_index()
  for sea in ["DJF","MAM","JJA","SON"]:
    if tmp[(tmp.season==sea)].empty: continue
    GradV["station"].append(ss)
    GradV["model"].append(nmodels[im])
    if not tmp[tmp.zref==1].empty:
     gradv=-tmp[(tmp.zref==4)&(tmp.season==sea)].cos.values[0]+tmp[(tmp.zref==1)&(tmp.season==sea)].cos.values[0]
    else:
     gradv=-tmp[(tmp.zref==4)&(tmp.season==sea)].cos.values[0]+tmp[(tmp.zref==2)&(tmp.season)].cos.values[0]
    GradV["cos"].append(gradv)
    GradV["season"].append(sea)
    GradV["lat"].append(lat_stat)
GradV=pd.DataFrame(GradV)



#1) Gradients par station
markers = {"obs": "k","TM5":"deepskyblue" ,"TOMCAT": "darkgreen","LMDZ":"red","TM3":"skyblue","MIROC4":"violet","NICAM5":"orange","NICAM6":"saddlebrown"}
fig, axes = plt.subplots(1, 1,figsize=(10,3))
GradV_St=GradV.groupby(["station","model"]).mean().reset_index()
GradV_St=GradV_St.sort_values("lat")
sns.swarmplot(x='station', y='cos', hue="model",s=6,data=GradV_St,ax=axes,palette=markers)
axes.tick_params(direction='in', length=3, color='k',grid_color="grey",colors="grey",
                grid_alpha=0.5,labelsize=13)
box = axes.get_position()
axes.set_position([box.x0, box.y0, box.width * 0.9, box.height])
axes.legend(loc='center left', bbox_to_anchor=(1., 0.8),frameon=False)
axes.set_xlabel("")
axes.set_ylabel("COS [ppt]",fontsize=13)
axes.set_title("Annual vertical gradient",fontsize=13)
plt.savefig('/home/users/mremaud/PYTHON/COS/TRANSCOM/FIG/GradV_stat_'+name_plot+'.png',format='png', bbox_inches='tight', dpi = 200)

#Seasonal variation
name_sort = {'DJF':0,'MAM':1,'JJA':2,"SON":3}
GradV['name_sort'] = GradV.season.map(name_sort)
GradV=GradV.sort_values('name_sort')

GradV2=GradV.groupby(["season","model"]).mean().reset_index()
GradV2["std"]=GradV.groupby(["season","model"]).std().reset_index().cos.values

fig, axes = plt.subplots(1, 1,figsize=(6,4))
sns.pointplot(x='season', y='cos', hue="model",data=GradV,ax=axes,palette=markers)
x = axes.get_position()
axes.set_position([box.x0, box.y0, box.width * 0.9, box.height])
axes.legend(loc='center left', bbox_to_anchor=(1., 0.8),frameon=False)
axes.set_xlabel("")
axes.set_ylabel("COS [ppt]",fontsize=13)
axes.set_title("Seasonal vertical gradient",fontsize=13)
axes.tick_params(direction='in', length=3, color='k',grid_color="grey",colors="grey",
                grid_alpha=0.5,labelsize=13)


axes.grid(alpha=0.2)
plt.savefig('/home/users/mremaud/PYTHON/COS/TRANSCOM/FIG/GradV_sea_'+name_plot+'.png',format='png', bbox_inches='tight', dpi = 200)

