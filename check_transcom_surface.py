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
models=["TOMCAT_ULeic","LMDZ_LSCE","TM5_UU","TM3_MPI","MIROC4_JAMSTEC","NICAM5_TM","NICAM6_TM"]
models=["LMDZ_LSCE"] #,"TOMCAT_ULeic"]
nmodels=["TOMCAT","LMDZ","TM5","TM3","MIROC4","NICAM5","NICAM6"]
colors=["darkgreen","red","blue","brown","violet","orange","saddlebrown"]
lines=["solid","dotted","dashdot","dashed","dashed"]
process={"Name":["Ocean","Soil","Opt","Anthro","Bioburn","Vegetation"],"process":["OCEA","SOIL","OPT","ANT","BB","VEG"]}
process=pd.DataFrame(process)

dir_result="/home/surface1/mremaud/COS/TRANSCOM/OUTPUTS/"
dir_flux="/home/surface1/mremaud/COS/INPUT_sflx/TRANSCOM_FLUX/COS_Flux_1x1/"

station="GIF"
components=os.listdir(dir_flux)

nn=1
for dd in components:
 if '._' in dd: continue
 nn+=1

f,ax=plt.subplots(3,2)
plt.show()
row=0 ; col=0
for dd in components:
 if ('._') in dd: continue
 list_flux=os.listdir(dir_flux+"/"+dd)
 print(list_flux)
 list_flux=[ff[:-12] for ff in list_flux if "phy"  in ff]
 #list_flux=[ff for ff in list_flux if not "Diur"  in ff]

 list_flux=[ff.replace("Diurn_","Diurn") for ff in list_flux]
 list_flux=[ff.replace("Month_","Month") for ff in list_flux]

 list_flux=np.unique(list_flux)
 compound=copy.copy(dd[:3])
 compteur_flux=0
 for flux in  list_flux:
  if '._' in flux: continue
  for pp in process.process.unique():
   if pp in flux[4:8]:
    name_flux=process[process.process==pp].Name.iloc[0]

  for im,mm in enumerate(models):
   print(dir_result+"/"+mm+"/TIER1/"+compound+"_"+flux+"_"+mm+"_TIER1.nc")
   if not os.path.exists(dir_result+"/"+mm+"/TIER1/"+compound+"_"+flux+"_"+mm+"_TIER1.nc"): continue
   print(dd)
   print(dir_result+"/"+mm+"/TIER1/"+compound+"_"+flux+"_"+mm+"_TIER1.nc")
   data=xr.open_dataset(dir_result+"/"+mm+"/TIER1/"+compound+"_"+flux+"_"+mm+"_TIER1.nc")
   data=data.to_dataframe()
   print(data.station.unique())
   data=data[data.station==station]
   data=data.sort_values(by='time')
   data.set_index("time",inplace=True)
   data=data.resample("M").mean().dropna()
   print(data,flux)
   namef=flux.replace(flux[4:8],"")
   if im == 1:
    ax[row,col].plot(data.index.to_pydatetime(),data.cos,label=flux[4:-6],color=colors[im],linestyle=lines[compteur_flux])
   else:
    ax[row,col].plot(data.index.to_pydatetime(),data.cos,color=colors[im],linestyle=lines[compteur_flux])
  ax[row,col].set_ylabel("COS [ppt]") 
  ax[row,col].set_xlim(datetime.datetime(2010,1,1),datetime.datetime(2019,12,31))
  ax[row,col].tick_params(rotation=60,labelsize=6,axis="x") 
  ax[row,col].set_title(name_flux)
  compteur_flux+=1
 ax[row,col].legend(fontsize=6)
 if  row<2: ax[row,col].set_xlabel(" ")
 if  col==1: ax[row,col].set_ylabel(" ")
 plt.show()  
 row+=1
 if row==3: 
   row=0
   col=1
plt.subplots_adjust(wspace=0.6,hspace=0.6)
plt.suptitle("COS - "+station)
plt.savefig('/home/users/mremaud/PYTHON/COS/TRANSCOM/FIG/check_surface.png',format='png', bbox_inches='tight', dpi = 200)

#############################DIURNAL###################################
compound="COS"
flux=["ORC","SIB4"]
list_station=['ALT','BRW','CGO','GIF','HFM','KUM','LEF','MHD','MLO','NWR','PSA','SMO','SPO','SUM','THD','WIS']
plt.figure()
f,ax=plt.subplots(4,4)
row=0; col=0
LEGEND_HEIGHT = 0.1
for stat in list_station:
 for im,mm in enumerate(models):
   for iff,ff in enumerate(flux):
     if not os.path.exists(dir_result+"/"+mm+"/TIER1/"+compound+"_COS-VEG-"+flux[iff]+"-Diurn_"+mm+"_TIER1.nc"): continue
     for itt,tt in enumerate(["VEG","SOIL"]):
      data_D2=xr.open_dataset(dir_result+"/"+mm+"/TIER1/"+compound+"_COS-"+tt+"-"+flux[iff]+"-Diurn_"+mm+"_TIER1.nc").to_dataframe()
      data_M2=xr.open_dataset(dir_result+"/"+mm+"/TIER1/"+compound+"_COS-"+tt+"-"+flux[iff]+"-Month_"+mm+"_TIER1.nc").to_dataframe()
      data_M2=data_M2[(data_M2.station==stat)&(data_M2.time.dt.year==2016)].copy(deep=True).set_index("time")
      data_D2=data_D2[(data_D2.station==stat)&(data_D2.time.dt.year==2016)].copy(deep=True).set_index("time")
      if itt==0:
       data_D=data_D2.copy(deep=True); data_M=data_M2.copy(deep=True)
      else:
       data_D.cos+=data_D2.cos; data_M.cos+=data_M2.cos
     if data_D.empty: continue
     data_D=data_D.resample("M").mean().dropna()
     data_M=data_M.resample("M").mean().dropna()
     ax[row,col].plot(np.arange(1,13),(data_D.cos-data_M.cos),color=colors[im],linestyle=lines[iff],label=mm+" "+ff)
 ax[row,col].set_xticks(np.arange(1,13,2))
 ax[row,col].set_xticklabels(np.arange(1,13,2))
 ax[row,col].set_xlim(1,12)
 ax[row,col].tick_params(rotation=60,labelsize=8,axis="x")
 if col==0 : ax[row,col].set_ylabel("COS [ppt]")
 if col != 0: ax[row,col].set_yticks([])
 if row !=3: ax[row,col].set_xticks([])

 ax[row,col].set_title(stat)
 ax[row,col].set_ylim(-18,18)
 col+=1
 if col==4: col=0; row+=1
 entries = collections.OrderedDict()
 for axes in ax.flatten():
  for handle, label in zip(*axes.get_legend_handles_labels()):
    entries[label] = handle
f.tight_layout(rect=(0, LEGEND_HEIGHT, 1, 1), h_pad=0.5, w_pad=0.5)
legend = f.legend( 
        entries.values(), entries.keys(), 
        loc='lower left', bbox_to_anchor=(0., LEGEND_HEIGHT-0.08),ncol=6)                                                                              


plt.savefig('/home/users/mremaud/PYTHON/COS/TRANSCOM/FIG/Diurnal.png',format='png', bbox_inches='tight', dpi = 200)

