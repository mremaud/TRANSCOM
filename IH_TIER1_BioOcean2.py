"""

Author: Marine Remaud
Aim: Impact of the ocean (DMS) and biosphere fluxes on the IH gradient.
     Transcom paper

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

###Suppress MLO
for station in GdIH.station.unique():
 if station=="MLO": continue
 for flux in GdIH.flux.unique():
  for model in GdIH.model.unique():
   for sea in GdIH.season.unique():
    mask=(GdIH.model==model)&(GdIH.flux==flux)&(GdIH.station==station)&(GdIH.season==sea)
    if GdIH[mask].empty: continue
    GdIH.loc[mask,"cos"]=GdIH.loc[mask,"cos"]-GdIH[(GdIH.season==sea)&(GdIH.model==model)&(GdIH.flux==flux)&(GdIH.station=="MLO")].cos.values[0]
    if model=="OBS": print(GdIH[(GdIH.season==sea)&(GdIH.model==model)&(GdIH.flux==flux)&(GdIH.station=="MLO")].cos.values[0])
GdIH.loc[GdIH.station=="MLO","cos"]=1
for sea in ["August","February","Annual"]:
 fig, ax = plt.subplots(figsize=(13,4))
 GdIH3=GdIH[GdIH.season==sea].copy(deep=True)
 sn.boxplot(x='station',y='cos',hue="flux",
            data=GdIH3[(GdIH.model!="OBS")],
            showfliers=False,
            linewidth=0.75,
            ax=ax,boxprops=dict(alpha=.3),palette={"Ctl": "grey", "Bio 2": "green","Ocean 2":"deepskyblue"})
 sn.regplot(x=np.arange(len(GdIH3[GdIH3.model=="OBS"])), y=np.array(GdIH3[GdIH3.model=="OBS"].cos.values), scatter=True, fit_reg=False, marker='*',color="k",scatter_kws={"s": 100},ax=ax,label="obs")

 ax.plot(np.arange(len(GdIH3[GdIH3.model=="OBS"])), GdIH3[GdIH3.model=="OBS"].cos, '--k*', lw=1)
 ax.legend()
 #ax.get_legend().remove()
 ax.set_xlabel("")
 ax.set_ylabel("COS [ppt]",fontsize=14,labelpad=0.02)

#ax.set_ylim(0,250)

 ax.tick_params(direction='in', length=3, color='k',grid_color="grey",colors="grey",
                grid_alpha=0.5,labelsize=13)
 plt.savefig('/home/users/mremaud/PYTHON/COS/TRANSCOM/FIG/TIER1_GdH_flux_'+sea+'.png',format='png', bbox_inches='tight', dpi = 400)

