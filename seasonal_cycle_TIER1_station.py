"""
Author: Marine Remaud
Smooth seasonal cycle of COS at each surface station for each ATM

"""
#Seasonal cycle at ALT and BRW
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


models=["OBS" ,"TOMCAT_ULeic","LMDZ_LSCE","LMDZ3_LSCE","TM5_UU","TM3_MPI","MIROC4_JAMSTEC","NICAM5_TM","NICAM6_TM"]
nmodels=["OBS","TOMCAT"      ,"LMDZ6"    ,"LMDZ3"     ,"TM5","TM3","MIROC4","NICAM5","NICAM6"]
#models=["LMDZ3_LSCE"]
#nmodels=["LMDZ3"]
colors=["k","darkgreen","red","magenta","blue","skyblue","violet","orange","saddlebrown"]
lines=["solid","dotted","dashdot","dashed","dashed"]
process={"Name":["Ocean","Soil","Opt","Anthro","Bioburn","Vegetation"],"process":["OCEA","SOIL","OPT","ANT","BB","VEG"]}
process=pd.DataFrame(process)

dir_result="/home/surface1/mremaud/COS/TRANSCOM/OUTPUTS/"

f,ax=plt.subplots(1,2,figsize=(8,3))
plt.show()
compound="COS"
list_flux=["BB-Stin","SOIL-SIB4",'VEG-SIB4','OCEAN-Lennartz','ANTHR-Zumkehr']
list_flux=["OPT-JM"]
begy=2011
endy=2018
row=0
for station in ['BRW','MLO']: #'ALT',"KUM",'BRW','SUM','MHD','LEF','MLO','CGO','HFM','THD','GIF','NWR','SMO','PSA','SPO']:
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
  data2=data2[(data2.date.dt.year>=begy)&(data2.date.dt.year<=endy)]
  data2=data2.groupby(data2.date.dt.month).mean()

  if (not data2.empty)&(mm!="LMDZ3_LSCE"): ax[row].plot(data2.index,data2.meas,color=colors[im],label=nmodels[im])
  if mm== "LMDZ3_LSCE": ax[row].plot(data2.index,data2.meas,color=colors[im],label=nmodels[im],linewidth=4)
 if row==0: ax[row].set_ylabel("COS [ppt]",fontsize=14) 
 ax[row].set_xlim(1,12)
 ax[row].tick_params(rotation=60,labelsize=6,axis="x") 
 ax[row].set_title(station)
 if row==1: ax[row].legend(fontsize=8,frameon=False)
# ax[row].set_xlabel("Month",fontsize=14)
 ax[row].tick_params(direction='in', length=3, color='k',grid_color="grey",colors="grey",
                grid_alpha=0.5,labelsize=13)  
 ax[row].set_xticks(np.arange(1,13,1))

 row+=1
plt.show()
plt.subplots_adjust(wspace=0.25)
plt.savefig('/home/users/mremaud/PYTHON/COS/TRANSCOM/FIG/TIER1_SC_ctl.png',format='png', bbox_inches='tight', dpi = 200)

