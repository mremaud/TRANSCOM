import pandas as pd
import numpy as np
import os
import xarray as xr
from netCDF4 import Dataset
import datetime
"""
Marine Remaud
Courbes du papier dans Atmosphere
Cycle saisonnier moyen de chaque composante du budget atmospherique du COS avec LMDz
Avec la chimie (Puit OH+ Puit strato)
"""

from useful import *
import matplotlib.pyplot as plt
import copy
import collections
import seaborn as sn
begy=2011
endy=2015

colors=["k","darkgreen","red","blue","skyblue","violet","orange","saddlebrown"]
lines=["solid","dotted","dashdot","dashed","dashed"]
dir_result="/home/surface1/mremaud/COS/TRANSCOM/OUTPUTS/"
compound="COS"
list_flux=["BB-Stin","SOIL2-ORC",'VEG-ORC','OCEAN-LennartzPICSES','ANTHR-Zumkehr']
#list_flux=["BB-Stin","SOIL2-ORC",'VEG-ORC','OCEAN-OPT','ANTHR-Zumkehr']
name_flux=["Burn","Soil","Plant","Ocean","Anth"]
color_flux=["red","darkgoldenrod","darkgreen","skyblue","orange"]
models=["LMDZ_LSCE"]
list_station=['BRW'] #,'GIF','MHD'] #'BRW','SUM','MHD','LEF','MLO','CGO','HFM','THD','GIF','NWR','SMO','PSA','SPO']
for station in list_station: 
  fig, ax = plt.subplots()
  ax.set_title("b) Seasonal Cycle",fontsize=18)
  obs=pd.read_pickle("/home/surface1/mremaud/COS/SURFACE/STATION/"+station+".pkl")
  obs = obs.rename({'date': 'time', 'obs': 'cos'}, axis='columns')
  obs.cos*=10**6
  #obs.set_index("time",inplace=True)
  obs['frac']=obs.apply(lambda row: fraction_an(row),axis=1)
  file_out=station+'.txt'
  os.system('rm -f '+file_out)
  np.savetxt(file_out,obs[['frac','cos']],fmt='%4.8f %3.3f')
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
  data2.plot(y="meas",ax=ax,label="Obs",color="darkgrey",linewidth=2,linestyle="--")
  for iff,flux in enumerate(list_flux):
   if flux == "OBS": continue
   data3=pd.DataFrame()
   for model in models:
    print(dir_result+"/"+model+"/TIER11/"+compound+"_"+compound+"-"+flux+"-Month2_"+model+"_TIER1.nc")
    if not os.path.exists(dir_result+"/"+model+"/TIER11/"+compound+"_"+compound+"-"+flux+"-Month2_"+model+"_TIER1.nc"): continue
    data=xr.open_dataset(dir_result+"/"+model+"/TIER11/"+compound+"_"+compound+"-"+flux+"-Month2_"+model+"_TIER1.nc")
    data=data.to_dataframe()
    data=data[data.station==station]
    data=data.sort_values(by='time')
    data.set_index("time",inplace=True)
    data.reset_index(inplace=True)
    if data.empty: continue
    ####CCGV#################################
    data['frac']=data.apply(lambda row: fraction_an(row),axis=1)
    file_out=station+'.txt'
    os.system('rm -f '+file_out)
    np.savetxt(file_out,data[['frac','cos']],fmt='%4.8f %3.3f')
    file_fit=station+'2.txt'
    #Model: without cal and ori
    os.system('./ccgcrv -cal -interv 7 -short 80 -equal -smooth -npoly 2 -nharm 4  -func -trend  -f  '+file_fit+' '+file_out)
    cgv=np.loadtxt(file_fit)
    data2=pd.DataFrame()
    data2['meas']=np.copy(cgv[:,4] -cgv[:,5])
    data2['date']=[datetime.datetime(int(cgv[ii,0]),int(cgv[ii,1]),int(cgv[ii,2])) for ii in range(len(cgv))]
    data2=data2[(data2.date.dt.year<=endy)&(data2.date.dt.year>=begy)]
    data2=pd.DataFrame(data2)
    data2.set_index("date",inplace=True)
    data2=data2.resample("M").mean().dropna().reset_index()
    data2=data2.groupby(data2.date.dt.month).mean().reset_index()
    data2["model"]=model
    data3=data3.append(data2)
   data3=data3.groupby("date").mean()
   data3.plot(y="meas",ax=ax,label=name_flux[iff],color=color_flux[iff],linewidth=2)
###Add net emissions
  data3=pd.DataFrame()
  for model in models:
   for iff,flux in enumerate(list_flux):
    if flux == "OBS": continue
    print(dir_result+"/"+model+"/TIER11/"+compound+"_"+compound+"-"+flux+"-Month2_"+model+"_TIER1.nc")
    if not os.path.exists(dir_result+"/"+model+"/TIER11/"+compound+"_"+compound+"-"+flux+"-Month2_"+model+"_TIER1.nc"): continue
    data=xr.open_dataset(dir_result+"/"+model+"/TIER11/"+compound+"_"+compound+"-"+flux+"-Month2_"+model+"_TIER1.nc")
    data=data.to_dataframe()
    data=data[data.station==station]
    data=data.sort_values(by='time')
    data.set_index("time",inplace=True)
    if iff == 0:
     dataf=data.copy(deep=True)
    else:
     dataf.cos+=data.cos
   dataf.reset_index(inplace=True)
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
   data2=data2[(data2.date.dt.year<=endy)&(data2.date.dt.year>=begy)]
   data2=pd.DataFrame(data2)
   data2.set_index("date",inplace=True)
   data2=data2.resample("M").mean().dropna().reset_index()
   data2=data2.groupby(data2.date.dt.month).mean().reset_index()
   data2["model"]=model
   data3=data3.append(data2)
  data4=data3.groupby("date").mean()
  data4["std"]=data3.groupby("date").std().meas
  data4.plot(y="meas",ax=ax,label="Net",color="k",linewidth=2)
  plt.fill_between(data4.index, data4.meas.values-data4["std"].values, data4.meas.values+data4["std"].values,alpha=0.1,color="k",axes=ax)

  ax.set_xlabel("Month",color="grey",fontsize=16)
  ax.set_ylabel("COS [ppt]",fontsize=16,color="grey")
  ax.tick_params(direction='in', length=6, width=2, color='k',grid_color="grey",colors="grey",
                grid_alpha=0.5,labelsize=15)
  ax.set_xlim(1,12)
  ax.set_xticks(np.arange(1,13,1))
  plt.legend(frameon="False",fontsize=12)
  plt.savefig('/home/users/mremaud/PYTHON/COS/TRANSCOM/FIG/CycleS_OCE-Lenn_Atmosphere_'+station+'_'+str(begy)+str(endy)+'.png',format='png', bbox_inches='tight', dpi = 200)

