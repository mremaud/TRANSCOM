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
models=["LMDZ_LSCE","TM5_UU","TM3_MPI","MIROC4_JAMSTEC","NICAM6_TM","NICAM5_TM"]
#models=["TM5_UU"]
nmodels=["LMDZ","TM5","TM3","MIROC4","NICAM6","NICAM5"]
colors=["red","blue","orange","violet","peru","saddlebrown"]
lines=["solid","dotted","dashdot","dashed","dashed"]
process={"Name":["Opt"],"process":["OPT"]}
process=pd.DataFrame(process)

zref=np.arange(500,11500,1000)
zref2=np.arange(0,10,2)
dir_result="/home/surface1/mremaud/COS/TRANSCOM/OUTPUTS/"

####FIND THE FIXED STATION
list_stations=[]
obs=pd.read_csv("/home/surface1/mremaud/COS/TRANSCOM/TRANSCOM_OBS/DATA_TIER2/TIER2.csv")
obs["date (UTC)"]=pd.to_datetime(obs["date (UTC)"])
obs=obs[obs["date (UTC)"].dt.year>=2010]
obs=obs.groupby(["lat","flight"]).mean().reset_index()
obs=obs[obs.lat>0]
obs=obs.sort_values("lat")
for stat in obs.flight.unique():
 if not np.abs(obs[obs.flight==stat].lat.min()-obs[obs.flight==stat].lat.max())>1:
  list_stations.append(stat)
print(list_stations)
f,ax=plt.subplots(4,3,figsize=(20,8))
row=0 ; col=0

which=["JM","LSCE"]
for iss,ss in enumerate(list_stations):

 ##Open observations###############################################################
 obs_stat=pd.read_pickle("/home/satellites16/mremaud/OBS/"+ss.lower()+".pkl")
 obs_stat.rename(columns={"date":"time","COS":"cos"},inplace=True)
 obs_stat["zref"]=obs_stat.apply(lambda row: zref[np.abs(row.alt-zref).argmin()],axis=1)
 obs_stat.zref/=1000
 obs_stat2=obs_stat.sort_values(by='time')
 obs_stat2=obs_stat2[obs_stat2.time.dt.year>2009]
 obs_stat2=obs_stat2[obs_stat2.lon>-900]
 obs_stat=obs_stat2.groupby(["time","zref"]).mean().reset_index()
 obs_stat=obs_stat.set_index("time").groupby("zref").resample("M").mean()
 del obs_stat["zref"];obs_stat.reset_index(inplace=True)
 obs_stat=obs_stat.groupby(["zref"]).mean().reset_index()
 mean_obs=np.copy(obs_stat[obs_stat.zref==3.5].cos.iloc[0])
 mean_obs=np.copy(obs_stat.mean().cos)
 ax[row,col].plot(obs_stat.cos,obs_stat.zref,label="obs",color="k")
 ##################################################################################

 for ww,wh in enumerate(which):
  for im,mm in enumerate(models):
   name_file=dir_result+"/"+mm+"/TIER2/COS_COS-OPT-"+wh+"-Month_"+mm+"_TIER2.nc"
   print(name_file)
   if not os.path.exists(name_file): continue
   data=xr.open_dataset(name_file)
   data=data.to_dataframe()
   data=data[data.time.dt.year>2009]
   data=data[data.lon>-999]
   print(data)
   data=data[data.station==ss]
   if data.empty: continue
   data["zref"]=data.apply(lambda row: zref[np.abs(row.alt-zref).argmin()],axis=1)
   data.zref/=1000
   data2=data.sort_values(by='time')
   data=data2.groupby(["time","zref"]).mean().reset_index()
   data=data.set_index("time").groupby("zref").resample("M").mean()
   del data["zref"];data.reset_index(inplace=True)
   data=data.groupby(["zref"]).mean()
   mean_mod=data[data.index==3.5].cos.iloc[0]
   mean_mod=data.mean().cos

   data.cos=data.cos-mean_mod+mean_obs
   print(data[data.index==4.5].cos,mean_obs)
   ax[row,col].plot(data.cos,data.index,label=mm+" "+wh,color=colors[im],linestyle=lines[ww])
 ax[row,col].set_ylabel("Altitude [km]") 
 ax[row,col].set_ylim(0.5,8.5)
 ax[row,col].set_xlim(440,535)

 ax[row,col].set_yticks(zref2)
 ax[row,col].tick_params(rotation=30,labelsize=8,axis="x") 
 ax[row,col].set_title(ss,fontsize=10)
 if (row==3) &(col==1): ax[row,col].legend(fontsize=8,bbox_to_anchor=(2.4,-3,0.,0.5),ncol=2, loc="upper right")
 if (row!=3)&(col!=2): ax[row,col].get_xaxis().set_ticks([])
 if (row!=2)&(col==2): ax[row,col].get_xaxis().set_ticks([])
 if row==3: ax[row,col].set_xlabel("COS [ppt]")

 if  row<3: ax[row,col].set_xlabel(" ")
 if  col>=1: ax[row,col].set_ylabel(" ")
 # plt.show()  
 row+=1
 if (row==4)&(col<2): 
   row=0
   col+=1
 if (row==4)&(col==2):
   row=0
   col=0
   f,ax=plt.subplots(3,3)

plt.subplots_adjust(wspace=0.1,hspace=0.2)
f.delaxes(ax[3][2])
plt.suptitle("Optimized fluxes")

plt.savefig("NOAA-prof.pdf",format='pdf')
