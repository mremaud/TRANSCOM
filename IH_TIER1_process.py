"""
Author: Marine Remaud
Aim   : Contribution of the COS budget component to the latitudinal gradient of COS at the NOAA network
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
import seaborn as sn

dir_result="/home/surface1/mremaud/COS/TRANSCOM/OUTPUTS/"
compound="COS"


models=["TOMCAT_ULeic","LMDZ_LSCE","TM5_UU","TM3_MPI","MIROC4_JAMSTEC","NICAM5_TM","NICAM6_TM"]
nmodels=["TOMCAT","LMDZ","TM5","TM3","MIROC4","NICAM5","NICAM6"]
colors=["darkgreen","red","blue","skyblue","violet","orange","saddlebrown"]

process={"Name":["Ocean","Soil","Opt","Anthro","Bioburn","Vegetation"],"process":["OCEA","SOIL","OPT","ANT","BB","VEG"]}
process=pd.DataFrame(process)


list_flux=["BB-UU","SOIL-SIB4",'VEG-SIB4','OCEAN-Lennartz','ANTHR-Zumkehr']
Gd={"cos":[],'lat':[],"station":[],'process':[],"model":[]}
name_flux=["Burn","Soil","Plant","Ocean","Anth"]
color_flux=["red","darkgoldenrod","darkgreen","skyblue","orange"]

for station in ['ALT',"KUM",'BRW','SUM','MHD','LEF','MLO','CGO','HFM','THD','NWR','SMO','PSA','SPO']:
 lat_station=pd.read_pickle("/home/surface1/mremaud/COS/SURFACE/STATION/"+station+".pkl").lat.iloc[0]
 for model in models:
  for iff,flux in enumerate(list_flux):
   if not os.path.exists(dir_result+"/"+model+"/TIER1/"+compound+"_"+compound+"-"+flux+"-Month_"+model+"_TIER1.nc"): continue
   data=xr.open_dataset(dir_result+"/"+model+"/TIER1/"+compound+"_"+compound+"-"+flux+"-Month_"+model+"_TIER1.nc")
   data=data.to_dataframe()
   data=data[data.station==station]
   if data.empty: continue
   data=data.sort_values(by='time')
   data2=ccgv(data,0)
   data2=data2[(data2.date.dt.year<=2018)&(data2.date.dt.year>=2011)]
   data2.set_index("date",inplace=True)
   data2=data2.resample("M").mean().dropna().reset_index()
   Gd["lat"].append(lat_station)
   Gd["station"].append(station)
   Gd["cos"].append(data2.meas.mean())
   Gd["process"].append(flux)
   Gd["model"].append(model)
   data.set_index("time",inplace=True)
   if iff==0:
    dataf=data.copy(deep=True)
   else:
    dataf.cos+=data.cos
  if "level_0" not in dataf.columns :
   dataf.reset_index(inplace=True)
   data2=ccgv(dataf,0)
   data2=data2[(data2.date.dt.year<=2018)&(data2.date.dt.year>=2011)]
   data2.set_index("date",inplace=True)
   data2=data2.resample("M").mean().dropna().reset_index()
   Gd["lat"].append(lat_station)
   Gd["station"].append(station)
   Gd["cos"].append(data2.meas.mean())
   Gd["process"].append("Net")
   Gd["model"].append(model)
Gd=pd.DataFrame(Gd)
Gd=Gd.sort_values(by="lat",ascending=True)
Gd2=Gd.groupby(["lat","process"]).mean()
Gd2["STD"]=Gd.groupby(["lat","process"]).std().cos
Gd2.reset_index(inplace=True)

fig, ax = plt.subplots()
for iff,flux in enumerate(list_flux):
  Gd2.loc[Gd2.process==flux,"cos"]-=Gd2[Gd2.process==flux].cos.mean()
  ax.plot(Gd2[Gd2.process==flux].lat, Gd2[Gd2.process==flux].cos, color=color_flux[iff],linewidth=2,ms=5, label=name_flux[iff])
Gd2.loc[Gd2.process=="Net","cos"]-=Gd2[Gd2.process=="Net"].cos.mean()
ax.plot(Gd2[Gd2.process=="Net"].lat, Gd2[Gd2.process=="Net"].cos, color="dimgrey",linewidth=3,ms=5, label="Net")
plt.fill_between(Gd2[Gd2.process=="Net"].lat,Gd2[Gd2.process=="Net"].cos-Gd2[Gd2.process=="Net"].STD.values,Gd2[Gd2.process=="Net"].cos+Gd2[Gd2.process=="Net"].STD.values,axes=ax,color="k",alpha=0.1)
ax.set_xlabel("Latitude",color="grey",fontsize=15)
ax.set_ylabel("COS [ppt]",fontsize=15,color="grey")
ax.tick_params(direction='in', length=6, width=2, color='k',grid_color="grey",colors="grey",
                grid_alpha=0.5,labelsize=13)
ax.set_xticks(np.arange(-80,100,20))
ax.legend(frameon=False,fontsize=12)
ax.set_title("a) Latitudinal distribution",fontsize=18)
plt.savefig('/home/users/mremaud/PYTHON/COS/TRANSCOM/FIG/TIER1_IH_process.png',format='png', bbox_inches='tight',pad_inches = 0, dpi = 200)



