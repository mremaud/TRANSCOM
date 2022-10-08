#Map of the observations to be used in TRANSCOM
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import h5py
import datetime
from pyhdf.SD import SD, SDC
import pandas as pd
from scipy import stats
from sklearn.linear_model import LinearRegression
#os.environ['PROJ_LIB'] = r'/usr/local/install/python-3/pkgs/proj4-5.2.0-he6710b0_1/share/proj/'
from mpl_toolkits.basemap import Basemap, cm


dir_transcom="/home/surface1/mremaud/COS/TRANSCOM/TRANSCOM_OBS/DATA_TIER1/TIER1.csv"
obs=pd.read_csv(dir_transcom)
obs=obs.groupby("station").mean().reset_index()
obs=obs[obs.station!="IZA"]
#--------------MAP--------------------------------------------------------------------------------
plt.figure()
fig, ax = plt.subplots()
m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90, 
                llcrnrlon=-180,urcrnrlon=180,resolution='c')

zref=np.linspace(500,10500,5)
n_graticules = 18
parallels = np.arange(-90., 90, 30)
meridians = np.arange(-180, 180., 40)
lw = 0.5
dashes = [1,1] # 5 dots, 7 spaces... repeat
graticules_color = 'pink'

m.drawparallels(parallels, linewidth=0.01,labels=[1,0,0,0], dashes=dashes, color=graticules_color,zorder=1)
m.drawmeridians(meridians,linewidth=0.01, labels=[0,0,0,1],dashes=dashes, color=graticules_color,zorder=1)
#m.drawmapboundary(fill_color='white',zorder=20)
m.drawcoastlines(linewidth=0.5,zorder=10,color="grey")

station=obs.station
lons=obs.lon
lats=obs.lat

x, y = m(lons,lats)
l1=m.scatter(x,y,s=5,marker='s',color='red',zorder=80,alpha=1,label="NOAA stations")
#for label, xpt, ypt in zip(station, x, y):
#    plt.text(xpt, ypt, label,fontsize=9,color='r')

dir_transcom="/home/surface1/mremaud/COS/TRANSCOM/TRANSCOM_OBS/DATA_TIER3/TIER3.csv"
obs=pd.read_csv(dir_transcom)
obs=obs.groupby("station").mean().reset_index()
lons=obs.lon
lats=obs.lat

x, y = m(lons,lats)
l1=m.scatter(x,y,s=40,marker='*',label="FTIR stations",color='darkorange',zorder=80,alpha=1)

dir_transcom="/home/surface1/mremaud/COS/TRANSCOM/TRANSCOM_OBS/DATA_TIER2/TIER2.csv"
obs=pd.read_csv(dir_transcom)
obs["date (UTC)"]=pd.to_datetime(obs["date (UTC)"])
obs=obs[(obs["date (UTC)"].dt.year>=2010)&(obs["date (UTC)"].dt.year<=2019)]
obs=obs[(obs.flight=="HIPPO#3")|(obs.flight=="HIPPO#4")|(obs.flight=="HIPPO#5")]
print(obs.flight.unique())
obs=obs.groupby(["lon","lat"]).mean().reset_index()
lons=obs.lon
lats=obs.lat
x, y = m(lons,lats)
l1=m.scatter(x,y,s=4,label="HIPPO",marker='o',color='powderblue',zorder=40,alpha=1)

dir_transcom="/home/surface1/mremaud/COS/TRANSCOM/TRANSCOM_OBS/DATA_TIER2/TIER2.csv"
obs=pd.read_csv(dir_transcom) 
obs["date (UTC)"]=pd.to_datetime(obs["date (UTC)"])
obs=obs[(obs["date (UTC)"].dt.year>=2010)&(obs["date (UTC)"].dt.year<=2019)]
obs=obs[(obs.flight=="ATOM")]
print(obs.flight.unique())
obs=obs.groupby(["lon","lat"]).mean().reset_index()
lons=obs.lon
lats=obs.lat
x, y = m(lons,lats)
l1=m.scatter(x,y,s=4,marker='o',color='lightgreen',zorder=40,alpha=1,label="ATOM")

from mpl_toolkits.axes_grid.inset_locator import inset_axes

inset_axes = inset_axes(ax, 
                    width="40%", # width = 30% of parent_bbo
                    height=0.9, # height : 1 inch
                    loc=1)

m = Basemap(projection='cyl',llcrnrlat=0,urcrnrlat=80,
                llcrnrlon=-180,urcrnrlon=-40,resolution='c')

zref=np.linspace(500,10500,5)
n_graticules = 18
parallels = np.arange(-90., 90, 30)
meridians = np.arange(-180, 180., 40)
lw = 0.5
dashes = [1,1] # 5 dots, 7 spaces... repeat
graticules_color = 'pink'

m.drawparallels(parallels, linewidth=0.01,labels=[0,0,0,0], dashes=dashes, color=graticules_color,zorder=1)
m.drawmeridians(meridians,linewidth=0.01, labels=[0,0,0,0],dashes=dashes, color=graticules_color,zorder=1)
#m.drawmapboundary(fill_color='white',zorder=20)
m.drawcoastlines(linewidth=0.5,zorder=10,color="grey")

dir_transcom="/home/surface1/mremaud/COS/TRANSCOM/TRANSCOM_OBS/DATA_TIER2/TIER2.csv"
obs=pd.read_csv(dir_transcom)
obs["date (UTC)"]=pd.to_datetime(obs["date (UTC)"])
obs=obs[(obs["date (UTC)"].dt.year>=2010)&(obs["date (UTC)"].dt.year<=2019)]
obs=obs[(obs.flight!="ACT")&(obs.flight!="ATOM")&(obs.flight!="HIPPO#3")&(obs.flight!="HIPPO#4")&(obs.flight!="HIPPO#5")]
#print(obs.flight.unique())
obs=obs.groupby(["lon","lat"]).mean().reset_index()

lons=obs.lon
lats=obs.lat
x, y = m(lons,lats)
l1=m.scatter(x,y,s=4,marker='o',color='mediumpurple',zorder=40,alpha=1,label="Airborne stations")
inset_axes.legend(frameon=False,bbox_to_anchor=(1., -2.2))
#plt.title("Airborne stations")
ax.legend(ncol=4,frameon=False,bbox_to_anchor=(1., -0.07))

plt.savefig('FIG/map_obs.jpg',bbox_inches='tight')





