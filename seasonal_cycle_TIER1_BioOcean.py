"""

Author: Marine Remaud
Aim: Impact of the ocean (DMS) and biosphere fluxes on the seasonal cycle of COS at the NOAA surface stations

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
list_flux=["BB-Stin","SOIL-SIB4",'VEG-SIB4','OCEAN-Lennartz','ANTHR-Zumkehr']


CycleS={"flux":[],"lat":[],"model":[],"Amplitude":[],"MonthMin":[],"station":[]}
for station in ['ALT',"KUM",'BRW','SUM','MHD','LEF','WIS','MLO','CGO','HFM','THD','GIF','NWR','SMO','PSA','SPO']:
 for im,mm in enumerate(models):
  dataf=pd.DataFrame()
  if mm=="OBS":
    dataf=pd.read_pickle("/home/surface1/mremaud/COS/SURFACE/STATION/"+station+".pkl")
    dataf = dataf.rename({'date': 'time', 'obs': 'cos'}, axis='columns')
    dataf.cos*=10**6
    dataf.set_index("time",inplace=True)
    lat_station=dataf.lat.iloc[0]
  for iff,flux in enumerate(list_flux):
   if flux == "OBS": continue
   print(dir_result+"/"+mm+"/TIER1/"+compound+"_"+compound+"-"+flux+"-Month_"+mm+"_TIER1.nc")
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
  if dataf.empty: continue
   ####CCGV#################################
  dataf['frac']=dataf.apply(lambda row: fraction_an(row),axis=1)    
  file_out=station+'.txt'
  os.system('rm -f '+file_out)  
  np.savetxt(file_out,dataf[['frac','cos']],fmt='%4.8f %3.3f')
  file_fit=station+'2.txt'
  #Model: without cal and ori
  os.system('./ccgcrv -cal -interv 7 -short 80 -equal -smooth -npoly 2 -nharm 4  -func -trend  -f  '+file_fit+' '+file_out)
  cgv=np.loadtxt(file_fit)
  data2=pd.DataFrame()
  data2['meas']=np.copy(cgv[:,4] -cgv[:,5])
  data2['date']=[datetime.datetime(int(cgv[ii,0]),int(cgv[ii,1]),int(cgv[ii,2])) for ii in range(len(cgv))]
  data2=data2[(data2.date.dt.year<=2018)&(data2.date.dt.year>=2011)]
  data2=pd.DataFrame(data2)
  data2.set_index("date",inplace=True)

  data2=data2.resample("M").mean().dropna().reset_index()
  data2=data2.groupby(data2.date.dt.month).mean()
  CycleS["Amplitude"].append(data2.meas.max()-data2.meas.min())
  CycleS["station"].append(station)
  CycleS["model"].append(nmodels[im])
  CycleS["MonthMin"].append(data2[data2.meas==data2.meas.min()].index[0])
  CycleS["lat"].append(lat_station)
  CycleS["flux"].append("Ctl")

list_flux=["BB-Stin","SOIL-SIB4",'VEG-ORC','OCEAN-Lennartz','ANTHR-Zumkehr']
for station in ['ALT',"KUM",'BRW','SUM','MHD','LEF','WIS','MLO','CGO','HFM','THD','GIF','NWR','SMO','PSA','SPO']:
 for im,mm in enumerate(models):
  dataf=pd.DataFrame()
  if mm=="OBS":
    lat_station=pd.read_pickle("/home/surface1/mremaud/COS/SURFACE/STATION/"+station+".pkl").lat.iloc[0]
  for iff,flux in enumerate(list_flux):
   if flux == "OBS": continue
   print(dir_result+"/"+mm+"/TIER1/"+compound+"_"+compound+"-"+flux+"-Month_"+mm+"_TIER1.nc")
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
  if dataf.empty: continue
   ####CCGV#################################
  dataf['frac']=dataf.apply(lambda row: fraction_an(row),axis=1)
  file_out=station+'.txt'
  os.system('rm -f '+file_out)
  np.savetxt(file_out,dataf[['frac','cos']],fmt='%4.8f %3.3f')
  file_fit=station+'2.txt'
  #Model: without cal and ori
  os.system('./ccgcrv -cal -interv 7 -short 80 -equal -smooth -npoly 2 -nharm 4  -func -trend  -f  '+file_fit+' '+file_out)
  cgv=np.loadtxt(file_fit)
  data2=pd.DataFrame()
  data2['meas']=np.copy(cgv[:,4] -cgv[:,5])
  data2['date']=[datetime.datetime(int(cgv[ii,0]),int(cgv[ii,1]),int(cgv[ii,2])) for ii in range(len(cgv))]
  data2=data2[(data2.date.dt.year<=2018)&(data2.date.dt.year>=2011)]
  data2=pd.DataFrame(data2)
  data2.set_index("date",inplace=True)

  data2=data2.resample("M").mean().dropna().reset_index()
  data2=data2.groupby(data2.date.dt.month).mean()
  CycleS["Amplitude"].append(data2.meas.max()-data2.meas.min())
  CycleS["station"].append(station)
  CycleS["model"].append(nmodels[im])
  CycleS["MonthMin"].append(data2[data2.meas==data2.meas.min()].index[0])
  CycleS["lat"].append(lat_station)
  CycleS["flux"].append("Bio 2")


list_flux=["BB-Stin","SOIL-SIB4",'VEG-SIB4','OCEAN-LennartzPICSES','ANTHR-Zumkehr']
for station in ['ALT',"KUM",'BRW','SUM','MHD','LEF','WIS','MLO','CGO','HFM','THD','GIF','NWR','SMO','PSA','SPO']:
 for im,mm in enumerate(models):
  dataf=pd.DataFrame()
  if mm=="OBS":
    lat_station=pd.read_pickle("/home/surface1/mremaud/COS/SURFACE/STATION/"+station+".pkl").lat.iloc[0]
  for iff,flux in enumerate(list_flux):
   if flux == "OBS": continue
   print(dir_result+"/"+mm+"/TIER1/"+compound+"_"+compound+"-"+flux+"-Month_"+mm+"_TIER1.nc")
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
  if dataf.empty: continue
   ####CCGV#################################
  dataf['frac']=dataf.apply(lambda row: fraction_an(row),axis=1)
  file_out=station+'.txt'
  os.system('rm -f '+file_out)
  np.savetxt(file_out,dataf[['frac','cos']],fmt='%4.8f %3.3f')
  file_fit=station+'2.txt'
  #Model: without cal and ori
  os.system('./ccgcrv -cal -interv 7 -short 80 -equal -smooth -npoly 2 -nharm 4  -func -trend  -f  '+file_fit+' '+file_out)
  cgv=np.loadtxt(file_fit)
  data2=pd.DataFrame()
  data2['meas']=np.copy(cgv[:,4] -cgv[:,5])
  data2['date']=[datetime.datetime(int(cgv[ii,0]),int(cgv[ii,1]),int(cgv[ii,2])) for ii in range(len(cgv))]
  data2=data2[(data2.date.dt.year<=2018)&(data2.date.dt.year>=2011)]
  data2=pd.DataFrame(data2)
  data2.set_index("date",inplace=True)

  data2=data2.resample("M").mean().dropna().reset_index()
  data2=data2.groupby(data2.date.dt.month).mean()
  CycleS["Amplitude"].append(data2.meas.max()-data2.meas.min())
  CycleS["station"].append(station)
  CycleS["model"].append(nmodels[im])
  CycleS["MonthMin"].append(data2[data2.meas==data2.meas.min()].index[0])
  CycleS["lat"].append(lat_station)
  CycleS["flux"].append("Ocean 2")

CycleS=pd.DataFrame(CycleS)
CycleS=CycleS.sort_values(by="lat",ascending=True)
CycleS.loc[(CycleS.model=="OBS"),"flux"]="OBS"
CycleS=CycleS.sort_values(by="lat",ascending=True)

nmodels=["OBS","TOMCAT","LMDZ","TM5","TM3","MIROC4","NICAM5","NICAM6"]
fig, axes = plt.subplots(2, 1, figsize=(13, 4), sharex=True)
markers = {"OBS": "k", "TOMCAT": "darkgreen","LMDZ":"red","TM5":"blue","TM3":"skyblue","MIROC4":"violet","NICAM5":"orange","NICAM6":"saddlebrown"}
CycleS2=CycleS[(CycleS.model!="OBS")]
CycleS2=CycleS2.sort_values(by="lat",ascending=True)

PROPS = {
    'boxprops':{'facecolor':'none', 'edgecolor':'k'},
    'medianprops':{'color':'k'},
    'whiskerprops':{'color':'k'},
    'capprops':{'color':'k'}
}

sn.boxplot(x='station',y='Amplitude',hue="flux",
            data=CycleS2,
            showfliers=False,
            linewidth=0.75, 
            ax=axes[0],boxprops=dict(alpha=.3),palette={"Ctl": "grey", "Bio 2": "green","Ocean 2":"deepskyblue"})
sn.regplot(x=np.arange(len(CycleS[CycleS.model=="OBS"])), y=np.array(CycleS[CycleS.model=="OBS"].Amplitude.values), scatter=True, fit_reg=False, marker='*',color="k",scatter_kws={"s": 100},ax=axes[0],label="obs") 
axes[0].plot(np.arange(len(CycleS[CycleS.model=="OBS"])), CycleS[CycleS.model=="OBS"].Amplitude, '--k*', lw=1)

axes[0].set_ylabel('Amplitude [ppt]',fontsize=14)

axes[1].set_ylabel('Month of the COS minimum',fontsize=14)
#sn.boxplot(x='station',y='MonthMin',
#            data=CycleS2,
#            showfliers=False,
#            linewidth=0.75,
#            **PROPS,ax=axes[1])
#CycleS2=CycleS2[CycleS2.model=="TM5"]
CycleS2=CycleS2.groupby(["station","flux","MonthMin"]).count().reset_index()
for stat in CycleS2.station.unique():
 for flux in CycleS2.flux.unique():
  max_lat=CycleS2[(CycleS2.station==stat)&(CycleS2.flux==flux)].lat.max()
  CycleS2=CycleS2.drop( CycleS2[(CycleS2.station==stat)&(CycleS2.flux==flux)&(CycleS2.lat!=max_lat)].index)
  if len(CycleS2[(CycleS2.station==stat)&(CycleS2.flux==flux)])>1:
   CycleS2=CycleS2.drop( CycleS2[(CycleS2.station==stat)&(CycleS2.flux==flux)].index[0])
 CycleS2.loc[CycleS2.station==stat,"lat"]=CycleS[CycleS.station==stat].lat.iloc[0]
CycleS2=CycleS2.sort_values(by="lat",ascending=True)

sn.swarmplot(data=CycleS2, x="station", y="MonthMin",hue="flux",palette={"Ctl": "grey", "Bio 2": "green","Ocean 2":"deepskyblue"},ax=axes[1])
sn.regplot(x=np.arange(len(CycleS[CycleS.model=="OBS"])), y=np.array(CycleS[CycleS.model=="OBS"].MonthMin.values), scatter=True, fit_reg=False, marker='*',color="k",scatter_kws={"s": 100},ax=axes[1],label="obs") 
axes[1].plot(np.arange(len(CycleS[CycleS.model=="OBS"])), CycleS[CycleS.model=="OBS"].MonthMin, '--k*', lw=1)

axes[1].get_legend().remove()
axes[1].set_yticks(np.arange(1,13,2))

fig.subplots_adjust(hspace=0)
axes[1].set_xlabel("")
axes[0].set_ylim(0,250)

box = axes[1].get_position()
axes[1].set_position([box.x0, box.y0, box.width * 0.8, box.height])
box = axes[0].get_position()
axes[0].set_position([box.x0, box.y0, box.width * 0.8, box.height])
axes[0].legend(loc='center left', bbox_to_anchor=(1., 0.5))
box = axes[1].get_position()

axes[0].tick_params(direction='in', length=3, color='k',grid_color="grey",colors="grey",
                grid_alpha=0.5,labelsize=13)
axes[0].set_yticks(np.arange(0,300,50))
axes[1].tick_params(direction='in', length=3, color='k',grid_color="grey",colors="grey",
                grid_alpha=0.5,labelsize=13)
#axes[0].yaxis.set_label_position("right")
#axes[0].yaxis.tick_right()

axes[1].grid(linewidth=0.3)
axes[0].grid(linewidth=0.3)
plt.savefig('/home/users/mremaud/PYTHON/COS/TRANSCOM/FIG/TIER1_CycleS_flux.png',format='png', bbox_inches='tight', dpi = 200)

