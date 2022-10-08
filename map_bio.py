#Maps of the prescribed fluxes


import numpy as np
import os
from netCDF4 import Dataset
import datetime
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import xarray as xr
plt.rcParams.update({'font.size': 16})


list_flux=["VEG"] #,"SOIL"]
name_flux=["COS_Vegetation_1x1"] #,"COS_Soil_1x1"]
list_model=["SIB4","ORC"]
dirflux='/home/surface1/mremaud/COS/INPUT_sflx/TRANSCOM_FLUX/COS_Flux_1x1/'
nomvar="flux"
titre=['a) SIB4','b) ORC-SIB4','c) |ORC-SIB4|/SIB4']

for iff,which in enumerate(name_flux):
  
  for yy in range(2011,2019):
   dataset = xr.open_dataset(dirflux+which+"/COS-"+list_flux[iff]+"-SIB4-Month_"+str(yy)+".nc")
   longitude1=dataset.lon.values
   latitude1=dataset.lat.values

   if yy == 2011:
    tmp1=np.copy(dataset.flux.values)
   else:
    tmp1=dataset.flux.values+tmp1
   tmp1=tmp1/(2018-2011+1)
   tmp1=np.squeeze(np.mean(tmp1,axis=0))
  if iff==0:
   SIB4=np.copy(tmp1)
  else:
   SIB4=SIB4+tmp1
plt.figure()
m = Basemap(projection='cyl',llcrnrlat=-70,urcrnrlat=90,\
              llcrnrlon=-180,urcrnrlon=178,resolution='c')
m.drawcoastlines()
n_graticules = 18
parallels = np.arange(-90., 90, 30)
meridians = np.arange(0., 360., 60)
lw = 0.5
dashes = [1,1] # 5 dots, 7 spaces... repeat
graticules_color = 'grey'
m.drawparallels(parallels, linewidth=0.5,labels=[1,0,0,0], dashes=dashes, color=graticules_color,zorder=20)
m.drawmeridians(meridians,linewidth=0.5, labels=[0,0,0,1],dashes=dashes, color=graticules_color,zorder=20)
Y, X=np.meshgrid(longitude1,latitude1) 
Y, X = m(Y,X)
v=np.arange(-1.75,1.8,0.5)
plt.title("a) SIB4")
ax=m.contourf(Y,X,SIB4*10**6,v,cmap=plt.cm.PiYG_r,extend='both')
cb=plt.colorbar(ax,orientation='horizontal',extend='both')
cb.set_label('$10^6$ mmol/m2/year')
plt.savefig("FIG/SIB4.png",format="png")

for  iff,which in enumerate(name_flux):

  for yy in range(2011,2019):
   dataset = xr.open_dataset(dirflux+which+"/COS-"+list_flux[iff]+"-ORC-Month_"+str(yy)+".nc")
   longitude1=dataset.lon.values
   latitude1=dataset.lat.values

   if yy == 2011:
    tmp1=np.copy(dataset.flux.values)
   else:
    tmp1=dataset.flux.values+tmp1
   tmp1=tmp1/(2018-2011+1)
   tmp1=np.squeeze(np.mean(tmp1,axis=0))
  if iff==0:
   ORC=np.copy(tmp1)
  else:
   ORC=ORC+tmp1
plt.figure()
m = Basemap(projection='cyl',llcrnrlat=-70,urcrnrlat=90,\
              llcrnrlon=-180,urcrnrlon=178,resolution='c')
m.drawcoastlines()
n_graticules = 18 
parallels = np.arange(-90., 90, 30)
meridians = np.arange(0., 360., 60)
lw = 0.5
dashes = [1,1] # 5 dots, 7 spaces... repeat
graticules_color = 'grey'
m.drawparallels(parallels, linewidth=0.5,labels=[1,0,0,0], dashes=dashes, color=graticules_color,zorder=20)
m.drawmeridians(meridians,linewidth=0.5, labels=[0,0,0,1],dashes=dashes, color=graticules_color,zorder=20)
Y, X=np.meshgrid(longitude1,latitude1)  
Y, X = m(Y,X)

plt.title("b) ORC-SIB4")
ax=m.contourf(Y,X,(ORC-SIB4)*10**6,v,cmap=plt.cm.PiYG_r,extend='both')
cb=plt.colorbar(ax,orientation='horizontal',extend='both')
cb.set_label('$10^6$ mmol/m2/year')
plt.savefig("FIG/ORC-SIB4.png",format="png")


plt.figure()
  #ax=m.pcolormesh(Y,X,diff,cmap=plt.cm.PiYG_r,vmin=-1.2,vmax=1.2)
m= Basemap(projection='cyl',llcrnrlat=-70,urcrnrlat=90,\
              llcrnrlon=-180,urcrnrlon=178,resolution='c')
m.drawcoastlines()
n_graticules = 18
parallels = np.arange(-90., 90, 30)
meridians = np.arange(0., 360., 60)
lw = 0.5
dashes = [1,1] # 5 dots, 7 spaces... repeat
graticules_color = 'grey'
m.drawparallels(parallels, linewidth=0.5,labels=[1,0,0,0], dashes=dashes, color=graticules_color,zorder=20)
m.drawmeridians(meridians,linewidth=0.5, labels=[0,0,0,1],dashes=dashes, color=graticules_color,zorder=20)
Y, X=np.meshgrid(longitude1,latitude1)
Y, X = m(Y,X)
v=np.arange(-100,100,10)
plt.title("c) (ORC-SIB4)/SIB4")
ax=m.contourf(Y,X,(ORC-SIB4)/np.abs(SIB4)*100,v,cmap=plt.cm.PiYG_r,extend='both')
cb=plt.colorbar(ax,orientation='horizontal',extend='both')
plt.savefig("FIG/ORC-SIB4_perc.png",format="png")


