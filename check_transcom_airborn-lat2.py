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
obs.reset_index(inplace=True)
obs["nobs"]=obs.index
obs=obs.groupby(["lat","flight"]).mean().reset_index()
obs=obs[obs.lat>0]
obs=obs.sort_values("lat")

list_stations=["ATOM","HIPPO"]

f,ax=plt.subplots(2,1)
row=0 ; col=0

which=["JM","LSCE"]

for iss,ss in enumerate(list_stations):
 ##Open observations###############################################################
 obs=pd.read_csv("/home/surface1/mremaud/COS/TRANSCOM/TRANSCOM_OBS/DATA_TIER2/TIER2.csv")
 obs["flight"]=obs.apply(lambda row: "HIPPO" if "HIPPO" in row.flight else row.flight,axis=1)
 obs_stat=obs[obs.flight==ss].copy(deep=True)
 obs_stat=obs_stat.rename(columns={"date (UTC)":"time","Unnamed: 0":"nobs"})
 obs_stat2=obs_stat.sort_values(by='time')
 obs_stat2["time"]=pd.to_datetime(obs_stat2["time"])
 
 if ss == "ATOM":
  obs_stat2.loc[(obs_stat2.time.dt.year==2016),"flight"]="ATOM1"
  obs_stat2.loc[(obs_stat2.time.dt.month==1)|(obs_stat2.time.dt.month==2)&(obs_stat2.time.dt.year==2017),"flight"]="ATOM2"
  obs_stat2.loc[(obs_stat2.time.dt.month==9)|(obs_stat2.time.dt.month==10)&(obs_stat2.time.dt.year==2017),"flight"]="ATOM3"
  obs_stat2.loc[(obs_stat2.time.dt.month==4)|(obs_stat2.time.dt.month==5)&(obs_stat2.time.dt.year==2018),"flight"]="ATOM4"
  #obs_stat2=obs_stat2[obs_stat2.flight=="ATOM1"]
 # obs_stat2=obs_stat2[obs_stat2.lon<-52]
 obs_stat["lat"]=obs_stat2.apply(lambda row: latref[np.abs(row.lat-latref).argmin()],axis=1)
 obs_stat=obs_stat.groupby(["time","lat"]).mean().reset_index()
 obs_stat["time"]=pd.to_datetime(obs_stat["time"])
 obs_stat=obs_stat.set_index("time").groupby("lat").resample("M").mean().dropna()
 del obs_stat["lat"];obs_stat.reset_index(inplace=True)
 obs_stat=obs_stat.groupby(["lat"]).mean().reset_index()
 mean_obs=obs_stat.cos.mean()
 ax[row].plot(obs_stat.lat,obs_stat.cos,label="obs",color="k")
 ##################################################################################

 for ww,wh in enumerate(which):
  for im,mm in enumerate(models):
   name_file=dir_result+"/"+mm+"/TIER2/COS_COS-OPT-"+wh+"-Month_"+mm+"_TIER2.nc"
   print(name_file)
   if not os.path.exists(name_file): continue
   data=xr.open_dataset(name_file)
   data=data.to_dataframe()
   data=data[data.cos>-900]
   data["station"]=data.apply(lambda row: "HIPPO" if "HIPPO" in row.station else row.station,axis=1)
   data=data[data.station==ss]
   if data.empty: continue
   if ss == "ATOM":
    data.loc[(data.time.dt.year==2016),"station"]="ATOM1"
    data.loc[(data.time.dt.month==1)|(data.time.dt.month==2)&(data.time.dt.year==2017),"station"]="ATOM2"
    data.loc[(data.time.dt.month==9)|(data.time.dt.month==10)&(data.time.dt.year==2017),"station"]="ATOM3"
    data.loc[(data.time.dt.month==4)|(data.time.dt.month==5)&(data.time.dt.year==2018),"station"]="ATOM4"
    data=data[data.station=="ATOM4"]
   # data=data[data.lon<-52]
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
   ax[row].plot(data.index,data.cos-bias,label=mm+" (OPT: "+wh+")",color=colors[im],linestyle=lines[ww])
 ax[row].set_xlabel("Latitude") 
 ax[row].set_ylim(420,540)
 ax[row].set_xlim(-90,90)

 #ax[row].set_yticks(y2)
 ax[row].tick_params(labelsize=10,axis="x") 
 ax[row].set_title(ss,fontsize=10)
 if (row==1): ax[row].legend(fontsize=6,bbox_to_anchor=(0.5,-0.8), loc="lower center",ncol=4)
 ax[row].set_ylabel("COS [ppt]")

 if  row==0: ax[row].set_xlabel(" ")
 # plt.show()  
 row+=1

plt.subplots_adjust(wspace=0.1,hspace=0.05)
plt.suptitle("Optimized fluxes")
#Legend outside the plot
p0 = ax[0].get_position()
p1 = ax[1].get_position()
p12 = [p1.x0, p1.y0+0.3*p1.height, p1.width*1, p1.height*0.7]
p02 = [p0.x0, p0.y0+0.3*p0.height, p0.width*1, p0.height*0.7]
ax[0].set_position(p02)
ax[1].set_position(p12)
plt.savefig("G-IH-ATOM4.pdf",format="pdf")


