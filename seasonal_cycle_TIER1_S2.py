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

dir_result="/home/surface1/mremaud/COS/TRANSCOM/OUTPUTS/"
compound="COS"
list_flux=["BB-Stin","SOIL2-ORC",'VEG-ORC','OCEAN-LennartzPICSES','ANTHR-Zumkehr']

list_flux=["BB-UU","SOIL-SIB4",'VEG-SIB4','OCEAN-Lennartz','ANTHR-Zumkehr']
#list_flux=["BB-UU","SOIL-ORC",'VEG-ORC','OCEAN-LennartzPICSES','ANTHR-Zumkehr']
name_flux=["Burn","Soil","Plant","Ocean","Anth"]
color_flux=["red","darkgoldenrod","darkgreen","skyblue","orange"]
models=["TOMCAT_ULeic","LMDZ_LSCE","TM5_UU","TM3_MPI","MIROC4_JAMSTEC","NICAM5_TM","NICAM6_TM"]
#models=["LMDZ_LSCE"]
list_station=['ALT','SUM','MHD','LEF','HFM','NWR','THD',"MHD","GIF",'MLO','KUM','CGO','SMO','PSA','SPO']

f,ax=plt.subplots(5,3,figsize=(8.27, 11.69))
row=0
col=0
for station in list_station:
  ax[row,col].set_title(station,fontsize=15)
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
  data2.plot(y="meas",ax=ax[row,col],label="Obs",color="darkgrey",linewidth=2,linestyle="--")
  for iff,flux in enumerate(list_flux):
   if flux == "OBS": continue
   data3=pd.DataFrame()
   for model in models:
    print(dir_result+"/"+model+"/TIER1/"+compound+"_"+compound+"-"+flux+"-Month_"+model+"_TIER1.nc")
    if not os.path.exists(dir_result+"/"+model+"/TIER1/"+compound+"_"+compound+"-"+flux+"-Month_"+model+"_TIER1.nc"): continue
    data=xr.open_dataset(dir_result+"/"+model+"/TIER1/"+compound+"_"+compound+"-"+flux+"-Month_"+model+"_TIER1.nc")
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
    data2=data2[(data2.date.dt.year<=2018)&(data2.date.dt.year>=2011)]
    data2=pd.DataFrame(data2)
    data2.set_index("date",inplace=True)
    data2=data2.resample("M").mean().dropna().reset_index()
    data2=data2.groupby(data2.date.dt.month).mean().reset_index()
    data2["model"]=model
    data3=data3.append(data2)
   data3=data3.groupby("date").mean()
   data3.plot(y="meas",ax=ax[row,col],label=name_flux[iff],color=color_flux[iff],linewidth=2)
  data3=pd.DataFrame()
  for model in models:
   for iff,flux in enumerate(list_flux):
    if flux == "OBS": continue
    print(dir_result+"/"+model+"/TIER1/"+compound+"_"+compound+"-"+flux+"-Month_"+model+"_TIER1.nc")
    if not os.path.exists(dir_result+"/"+model+"/TIER1/"+compound+"_"+compound+"-"+flux+"-Month_"+model+"_TIER1.nc"): continue
    data=xr.open_dataset(dir_result+"/"+model+"/TIER1/"+compound+"_"+compound+"-"+flux+"-Month_"+model+"_TIER1.nc")
    data=data.to_dataframe()
    data=data[data.station==station]
    data=data.sort_values(by='time')
    data.set_index("time",inplace=True)
    if iff == 0:
     dataf=data.copy(deep=True)
    else:
     dataf.cos+=data.cos
   if dataf.empty: continue
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
   data2=data2[(data2.date.dt.year<=2018)&(data2.date.dt.year>=2011)]
   data2=pd.DataFrame(data2)
   data2.set_index("date",inplace=True)
   data2=data2.resample("M").mean().dropna().reset_index()
   data2=data2.groupby(data2.date.dt.month).mean().reset_index()
   data2["model"]=model
   data3=data3.append(data2)
  data4=data3.groupby("date").mean()
  data4["std"]=data3.groupby("date").std().meas
  data4.plot(y="meas",ax=ax[row,col],label="Net",color="k",linewidth=2,legend=None)
  ax[row,col].fill_between(data4.index, data4.meas.values-data4["std"].values, data4.meas.values+data4["std"].values,alpha=0.1,color="k")

  ax[row,col].set_xlabel(" ",color="grey",fontsize=16)
  if col==0: ax[row,col].set_ylabel("COS [ppt]",fontsize=12,color="grey")
  ax[row,col].tick_params(direction='in', length=6, width=2, color='k',grid_color="grey",colors="grey",
                grid_alpha=0.5,labelsize=10)
  ax[row,col].set_xlim(1,12)
  ax[row,col].set_xticks(np.arange(1,13,1))
  if (row==4)&(col==0):
   ax[row,col].legend(frameon="False",fontsize=12,ncol=3,bbox_to_anchor=(2.5, -0.3))
  else:
   ax[row,col].get_legend().remove()
  if col==2:
   row+=1
   col=0
  else:
   col+=1 
plt.subplots_adjust(wspace=0.25,hspace=0.4)
plt.savefig('/home/users/mremaud/PYTHON/COS/TRANSCOM/FIG/TIER1_CycleS_process.png',format='png', bbox_inches='tight', dpi = 200)

