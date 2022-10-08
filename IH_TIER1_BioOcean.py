"""
Author: Marine Remaud
Aim: Impact of the ocean (DMS) and biosphere fluxes on the interhemispheric gradient. The modeled COS are sampled at the NOAA stations.
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


models=["OBS","TOMCAT_ULeic","LMDZ_LSCE","TM5_UU","TM3_MPI","MIROC4_JAMSTEC","NICAM5_TM","NICAM6_TM"]
nmodels=["OBS","TOMCAT","LMDZ","TM5","TM3","MIROC4","NICAM5","NICAM6"]
colors=["k","darkgreen","red","blue","skyblue","violet","orange","saddlebrown"]
lines=["solid","dotted","dashdot","dashed","dashed"]
process={"Name":["Ocean","Soil","Opt","Anthro","Bioburn","Vegetation"],"process":["OCEA","SOIL","OPT","ANT","BB","VEG"]}
process=pd.DataFrame(process)

dir_result="/home/surface1/mremaud/COS/TRANSCOM/OUTPUTS/"

compound="COS"


####Gradient Inter hemispherique
models=["OBS","TOMCAT_ULeic","LMDZ_LSCE","TM5_UU","TM3_MPI","MIROC4_JAMSTEC","NICAM5_TM","NICAM6_TM"]
list_flux=["BB-Stin","SOIL-SIB4",'VEG-SIB4','OCEAN-Lennartz','ANTHR-Zumkehr']
dico_flux={"Ctl":["BB-Stin","SOIL-SIB4",'VEG-SIB4','OCEAN-Lennartz','ANTHR-Zumkehr'],"Bio 2":["BB-Stin","SOIL-SIB4",'VEG-ORC','OCEAN-Lennartz','ANTHR-Zumkehr'],"Ocean 2":["BB-Stin","SOIL-SIB4",'VEG-SIB4','OCEAN-LennartzPICSES','ANTHR-Zumkehr']}
nmodels=["OBS","TOMCAT","LMDZ","TM5","TM3","MIROC4","NICAM5","NICAM6"]
colors=["k","darkgreen","red","blue","skyblue","violet","orange","saddlebrown"]
GdIH={"station":[],"lat":[],"cos":[],"model":[],"flux":[],'season':[],"color":[]}
for station in ['ALT',"KUM",'BRW','SUM','MHD','LEF','WIS','MLO','CGO','HFM','THD','GIF','NWR','SMO','PSA','SPO']:
 for im,mm in enumerate(models):
  dataf=pd.DataFrame()
  if mm=="OBS":
    dataf=pd.read_pickle("/home/surface1/mremaud/COS/SURFACE/STATION/"+station+".pkl")
    dataf = dataf.rename({'date': 'time', 'obs': 'cos'}, axis='columns')
    dataf.cos*=10**6
    dataf.set_index("time",inplace=True)
    lat_station=dataf.lat.iloc[0]
    #dataf.reset_index(inplace=True)
  for kflux in ["Ctl","Bio 2","Ocean 2"]:
   list_flux=dico_flux[kflux]
   for iff,flux in enumerate(list_flux):
     if mm == "OBS": continue
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
   if (mm =="OBS")&(kflux!="Ctl"):continue 
   dataf.reset_index(inplace=True)
   if dataf.empty: continue
   ####CCGV#################################
   data2=ccgv(dataf,0)
   data2=data2[(data2.date.dt.year<=2018)&(data2.date.dt.year>=2011)]
   data2.set_index("date",inplace=True)
   data2=data2.resample("M").mean().dropna().reset_index()
   GdIH["lat"].append(lat_station)
   GdIH["cos"].append(data2.mean().meas)
   GdIH["model"].append(nmodels[im])
   GdIH["station"].append(station)
   GdIH["flux"].append(kflux)
   GdIH["season"].append("Annual")
   GdIH["color"].append(colors[im])

   GdIH["lat"].append(lat_station)
   GdIH["cos"].append(data2[data2.date.dt.month==8].mean().meas)
   GdIH["model"].append(nmodels[im])
   GdIH["station"].append(station)
   GdIH["flux"].append(kflux)
   GdIH["season"].append("August")
   GdIH["color"].append(colors[im])

   GdIH["lat"].append(lat_station)
   GdIH["cos"].append(data2[data2.date.dt.month==2].mean().meas)
   GdIH["model"].append(nmodels[im])
   GdIH["station"].append(station)
   GdIH["flux"].append(kflux)
   GdIH["season"].append("February")
   GdIH["color"].append(colors[im])
GdIH=pd.DataFrame(GdIH)
GdIH=GdIH.sort_values(by='lat',ascending=True)

#Suppress the biais
GdIH=GdIH[GdIH.station!="KUM"]


for sea in GdIH.season.unique():
 GdIH3=GdIH[GdIH.season==sea]
 for im,mm in enumerate(GdIH3.model.unique()):
  for ff in GdIH3.flux.unique():
   bias=GdIH[(GdIH.model=="OBS")&(GdIH.station!="GIF")&(GdIH.flux=="Ctl")].mean().cos-GdIH[(GdIH.model==mm)&(GdIH.station!="GIF")&(GdIH.flux==ff)].mean().cos
   GdIH3.loc[(GdIH3.model==mm)&(GdIH3.flux==ff),"cos"]+=np.copy(bias)

 GdIH3.loc[(GdIH3.flux=="Bio 2"),"lat"]+=1
 GdIH3.loc[(GdIH3.flux=="Ocean 2"),"lat"]-=1


 # Plot
 fig, ax = plt.subplots()
 ax.plot(GdIH3[GdIH3.model=="OBS"].lat, GdIH3[GdIH.model=="OBS"].cos, linestyle='-', ms=3, label="obs",color="k",marker="*")
 GdIH3=GdIH3[GdIH3.model!="OBS"]
 ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling

 GdIH2=GdIH3.groupby(["flux","station"]).mean().reset_index()
 GdIH2["std"]=GdIH3.groupby(["flux","station"]).std().reset_index().cos.values

 groups = GdIH2.groupby('flux')
 palette={"Ctl": "grey", "Bio 2": "green","Ocean 2":"deepskyblue"} 
 for name, group in groups:
    print(name,group)
    ax.errorbar(group.lat, group.cos,yerr=group["std"] ,linestyle='',fmt='o', ms=2, label=name,color=palette[name])
    

 #ax.legend()
 ax.set_xlabel("latitude",color="darkgrey",fontsize=15)
 ax.set_ylabel("COS [ppt]",color="darkgrey",fontsize=15)
 ax.tick_params(direction='in', length=6, width=1, color='k',grid_color="grey",colors="grey",
                grid_alpha=0.5,labelsize=14)
 ax.set_xticks(np.arange(-80,100,20))
 ax.set_ylim(350,650)
 ax.set_xlim((-92,90))

 #if sea=="February":ax.set_ylim(350,680)
 if sea != "February": plt.legend(frameon=False)

 ax2 = ax.twiny()
 ax2.tick_params(axis=u'x',length=3,labeltop="on",direction="in")
 ax2.set_xticks(GdIH3[GdIH3.flux=="Ctl"]['lat'].unique())
 station2=['SPO', 'PSA', 'CGO', 'SMO', 'MLO','WIS',' ', ' ',
       'HFM', ' ',' ','MHD', 'BRW', ' ', 'ALT']
 ax2.set_xticklabels(station2,rotation="vertical",fontsize=12,color="grey")
 ax2.set_xlim((-92,90))

 plt.title(sea,fontsize=15)


 plt.show()
 plt.savefig('/home/users/mremaud/PYTHON/COS/TRANSCOM/FIG/TIER1_IH'+sea+'.png',format='png', bbox_inches='tight', dpi = 200)
