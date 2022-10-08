import pandas as pd
import numpy as np
import os
import xarray as xr
from netCDF4 import Dataset
import datetime
from useful import *
import matplotlib.pyplot as plt
import matplotlib as mpl
import copy


dir_plot="/home/surface1/mremaud/COS/TRANSCOM/OUTPUTS/"
list_models=["LMDZ_LSCE","TM5_UU","NICAM5_TM","NICAM6_TM","MIROC4_JAMSTEC","TOMCAT_ULeic"]
name_models=["LMDZ","    TM5",    "NICAM5",   "NICAM6",    "MIROC",         "TOMCAT"]
list_flux=["COS-ANTHR-Zumkehr-Month","COS-BB-UU-Month","COS-SOIL-SIB4-Month","COS-VEG-SIB4-Month","COS-OCEAN-LennartzPICSES-Month"]
#list_flux=["COS-OPT-JM-Month"]

reference="average"
start_year=2012
end_year=2018
season="Annual"
suptitre="Opt-"+season
suptitre=season
nplot=len(list_models)


timet=pd.date_range(start=datetime.datetime(2010,1,1),end=datetime.datetime(2018,12,31),freq='M')
if season=="win":
  sel_t=np.where((timet.month==12)|(timet.month==1)|(timet.month==2)&(timet.year>=2012)&(timet.year<=2018))[0]
elif season=="sum":
  sel_t=np.where((timet.month==6)|(timet.month==7)|(timet.month==8)&(timet.year>=2012)&(timet.year<=2018))[0]
else:
  sel_t=np.where((timet.year>=2012)&(timet.year<=2018))[0]

"""
1) Reference flux
"""

if reference=="average":
  nplot=nplot+1
  for imm,model in enumerate(list_models):
   for iff,flux in enumerate(list_flux):
    tmp=xr.open_dataset(dir_plot+'/'+model+"/TIER4/COS_"+flux+"_"+model+"_TIER4.nc2")
    if model == "TM5_UU": tmp=tmp.mean("lon")
    if (iff==0)&(imm==0):
     ref=tmp.copy(deep=True)
    else:
     ref.cos.values=ref.cos.values+tmp.cos.values
  ref.cos.values/=len(list_models)
  ref=ref.isel(time=sel_t)
  ref=ref.mean("time")
  ref=ref.isel(lev=np.arange(1,25))
nrow=int(nplot/3)+1
ncol=3 
f,ax=plt.subplots(nrow,ncol)
X,Y=np.meshgrid(ref.lat.values,ref.pressure.values[:,0])
l2=ax[0,0].contourf(X,ref.pressure.values[:,:]/100.,ref.cos.values[:,:]+473.33)
ax[0,0].invert_yaxis()
ax[0,0].set_title("Ref",pad=0.05)
ax[0,0].set_xticks(np.arange(-80,100,40))
ax[0,0].set_yticks(np.flip(np.arange(200,1100,200)))
ax[0,0].tick_params(axis="both",direction="in")
ax[0,0].set_ylabel("Pa (hPa)")
ax[0,0].xaxis.set_ticklabels([])
xmin = ax[0,0].get_position().xmin
xmax= ax[0,0].get_position().xmax
ymax = ax[0,0].get_position().ymax
ymin = ax[0,0].get_position().ymin
cb_ax = f.add_axes([xmax+(xmax-xmin)*0.1, ymin, 0.01, (ymax-ymin)])
cbar=f.colorbar(l2, ax=ax[0,0], cax=cb_ax)
cbar.set_label('COS [ppt]')

level_antd=np.arange(-24,25,4)
row=1;col=0
list_models=[mm for mm in list_models if mm not in reference]
for imm,mm in enumerate(list_models):
  for iff,flux in enumerate(list_flux):
    tmp=xr.open_dataset(dir_plot+"/"+mm+"/TIER4/COS_"+flux+"_"+mm+"_TIER4.nc2")
    if mm=="TM5_UU": tmp=tmp.mean("lon")
    if iff==0:
     tot=tmp.copy(deep=True)
    else:
     tot.cos.values=tot.cos.values+tmp.cos.values
  tot=tot.isel(time=sel_t)
  tot=tot.mean('time')
  tot=tot.isel(lev=np.arange(1,25))
  #Contour levels
  cmap=plt.cm.RdBu_r
  bounds = level_antd
  norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

  l2=ax[row,int(col)].contourf(X,tot.pressure/100.,tot.cos-ref.cos,level_antd,norm=norm,cmap=cmap,extend="both")
  ax[row,int(col)].set_title(name_models[imm]+"-ref",pad=0.05)
  ax[row,int(col)].invert_yaxis()
  ax[row,int(col)].set_xticks(np.arange(-80,100,40))
  ax[row,int(col)].set_yticks(np.flip(np.arange(200,1100,200)))
  ax[row,int(col)].tick_params(axis="both",direction="in")
  if col==0: ax[row,int(col)].set_ylabel("Pa (hPa)")
  if col>=1: ax[row,int(col)].yaxis.set_ticklabels([]) 
  if  row<2: ax[row,int(col)].xaxis.set_ticklabels([]) 
  if row==2: ax[row,int(col)].set_xlabel("latitude")
  #COLORBAR
  if (row==2)&(col==1): 
   #Position
   col=int(col)
   xmin = ax[1,2].get_position().xmin
   xmax=  ax[1,2].get_position().xmax
   ymax = ax[1,col].get_position().ymax
   ymin = ax[row,col].get_position().ymin
   cbar_ax = f.add_axes([xmax+(xmax-xmin)*0.03, ymin, 0.02,  (ymax-ymin)])
   cbar=f.colorbar(l2, ax=ax[1:,:],cax=cbar_ax,cmap=cmap,norm=norm,ticks=bounds, boundaries=bounds,extend="both")
   #Colorbar
   cbar.ax.tick_params(color='k', direction='in')
  col+=1
  if col==3: col=0.; row+=1
f.subplots_adjust(wspace=0.09)
f.subplots_adjust(hspace=0.2)

ax[0,1].set_axis_off()
ax[0,2].set_axis_off()
#f.suptitle(suptitre)
plt.savefig("FIG/"+suptitre)
plt.show()
