#Maps of the prescribed fluxes


import numpy as np
import os
from netCDF4 import Dataset
import datetime
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import xarray as xr
plt.rcParams.update({'font.size': 16})


dirflux='/home/surface1/mremaud/COS/INPUT_sflx/TRANSCOM_FLUX/COS_Flux_1x1/COS_BioBurning_1x1/'
nomvar="flux"

for yy in range(2011,2016):
   dataset = xr.open_dataset(dirflux+"/COS-BB-Stin-Month_"+str(yy)+".nc")
   longitude1=dataset.lon.values
   latitude1=dataset.lat.values

   if yy == 2011:
    REF=np.copy(dataset.flux.values)
   else:
    REF=dataset.flux.values+REF
   REF=REF/(2016-2011+1)
   REF=np.squeeze(np.mean(REF,axis=0))

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
v=np.arange(-9.5,10.5,1)
plt.title("a) BB Stin")
ax=m.contourf(Y,X,REF*10**7,v,cmap=plt.cm.PiYG_r,extend='both')
cb=plt.colorbar(ax,orientation='horizontal',extend='both')
cb.set_label('$10^7$ mmol/m2/year')
plt.savefig("FIG/BB.png",format="png")


for yy in range(2011,2016):
   dataset = xr.open_dataset(dirflux+"/COS-BB-UU-Month_"+str(yy)+".nc")
   longitude1=dataset.lon.values
   latitude1=dataset.lat.values
   if yy == 2011:
    PIC=np.copy(dataset.flux.values)
   else:
    PIC=dataset.flux.values+PIC
   PIC=PIC/(2016-2011+1)
   PIC=np.squeeze(np.mean(PIC,axis=0))
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

v=np.arange(-1.1,1.15,0.2)
v=np.arange(-9.5,10.5,1)

plt.title("b) BB UU")
ax=m.contourf(Y,X,(PIC)*10**7,v,cmap=plt.cm.PiYG_r,extend='both')
cb=plt.colorbar(ax,orientation='horizontal',extend='both')
cb.set_label('$10^7$ mmol/m2/year')
plt.savefig("FIG/BB-diff.png",format="png")


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
plt.title("c) (DMS (PISCES)-DMS (Lennartz)) \n /DMS (Lennartz) ")
ax=m.contourf(Y,X,(PIC-REF)/np.abs(REF)*100,v,cmap=plt.cm.PiYG_r,extend='both')
cb=plt.colorbar(ax,orientation='horizontal',extend='both')
plt.savefig("FIG/diff-BB_perc.png",format="png")


