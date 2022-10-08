"""
Author: Marine Remaud
Goal:   Script to compute the average total column of OCS simulated by the transport models that participated in the TRANSCOM COS experiment.
1) The simulated profiles are interpolated vertically from the model grid to the retrival vertical resolution.
2) Removing of the stratospheric part (above 9.8 km)
3) Convolution of the vertical profile by the pressure (create_h)
4) Average of the vertical profil of COS

"""

import h5py
from pyhdf.SD import SD, SDC
from sklearn.linear_model import LinearRegression
from useful import *
import os
import copy
import pandas
import  xarray as xr
import datetime
from dateutil.relativedelta import relativedelta
from scipy.interpolate import interp1d
import numpy as np
from scipy import stats
from utils import *

dir_transcom="/home/surface1/mremaud/COS/TRANSCOM/OUTPUTS/"
dirobs="/home/surface1/mremaud/COS/XCOS/RETRIEVAL/"
remove_strat=0

#model="TOMCAT_ULeic"
#model="TM5_UU"
#model="LMDZ_LSCE"
model="NICAM5_TM"
toplevel= 9.8 #Tropopause level
#1) TM5_UU
list_files=os.listdir(dir_transcom+"/"+model+"/TIER3")
xcos={'xcos':[],'date':[],'obs':[],'stat':[]}
list_station_xocs=["name_station"] #["paramaribo","kiruna","bremen","ny.alesund","zugspitze","izana"]
list_station= [ff[:3].upper() for ff in  list_station_xocs]   

#Read the NDAC data per file
e=os.listdir(dirobs)
for iss,stat in enumerate(list_station_xocs):
  #Loop over the stations
  print(stat)
  for FILE_NAME in e:
   if "SAV" in FILE_NAME: continue
   if not stat in FILE_NAME: continue
   # Open file.
   print(FILE_NAME)
   if FILE_NAME[-1]=='5' :
    f = h5py.File(dirobs+FILE_NAME, 'r')
    a_group_key = list(f.keys())[0]
    height=f['ALTITUDE'][:]
    alt=f['ALTITUDE.INSTRUMENT'][:]
    pcol=f['OCS.MIXING.RATIO.VOLUME_ABSORPTION.SOLAR']
    punit="ppmv"
    pcol=pcol[:]
    p_surface=f['SURFACE.PRESSURE_INDEPENDENT'][:]
    pressure_o=f['PRESSURE_INDEPENDENT'][:]
    longitude=f['LONGITUDE.INSTRUMENT'][:][0]
    latitude =f['LATITUDE.INSTRUMENT'][:][0]
    time=f['DATETIME'][:]
   else:
    hdf = SD(dirobs+FILE_NAME, SDC.READ)
    pcol = hdf.select('OCS.MIXING.RATIO.VOLUME_ABSORPTION.SOLAR')
    punit=pcol.units
    pcol=pcol.get()
    p_surface=hdf.select('SURFACE.PRESSURE_INDEPENDENT')
    p_surface=p_surface.get()
    var=hdf.select('PRESSURE_INDEPENDENT')
    pressure_o=var.get()
    var=hdf.select('LONGITUDE.INSTRUMENT')
    longitude=var.get()[0]
    var=hdf.select('LATITUDE.INSTRUMENT')
    latitude=var.get()[0]
    var=hdf.select('ALTITUDE')
    height=var.get()
    var=hdf.select('ALTITUDE.INSTRUMENT')
    alt=var.get()[0]
    var=hdf.select( 'OCS.COLUMN_ABSORPTION.SOLAR')
    xcos=var.get()
    var=hdf.select( 'DATETIME')
    time=var.get()

   date=[datetime.datetime(2000,1,1)+datetime.timedelta(days=int(ii)) for ii in time]
   date=[date[ii]+datetime.timedelta(hours=round((time[ii]-int(time[ii]) )*24.)) for ii in range(len(time))]
   date=pd.to_datetime(date)
   begy=date[0].year;endy=date[-1].year
   ndata=len(date)
   #Pressure, and height of the cutoff level
   pcol_98=np.zeros(ndata)
   pressure_98=np.zeros(ndata)
   if  height.ndim == 1 :
     if height[0]<height[1]:
      height=height[::-1]
      pcol=pcol[:,::-1]
      pressure_o =pressure_o[:,::-1]
     ih=np.abs(height-toplevel).argmin() #Take only the troposphere
     ih_min=(height[:]>toplevel).argmin()-1
     ih_max=(height[:]<toplevel).argmax()
     X=np.asarray([height[ih_max],height[ih_min]])[:,np.newaxis]
     for i in range(ndata):
      Y=[pressure_o[i,ih_max],pressure_o[i,ih_min]]
      reg = LinearRegression().fit(X, Y)
      pressure_98[i]=reg.intercept_+toplevel*reg.coef_
      Y=[pcol[i,ih_max],pcol[i,ih_min]]
      reg = LinearRegression().fit(X, Y)
      pcol_98[i]=reg.intercept_+toplevel*reg.coef_
   else:
     if height[0,0]<height[0,1]:
      height=height[:,::-1]
      pcol=pcol[:,::-1]
      pressure_o =pressure_o[:,::-1]
     ih=[(np.abs(height[i,:]-toplevel)).argmin() for i in range(ndata)]
     ih_min=[(height[i,:]>toplevel).argmin()-1 for i in range(ndata)]
     ih_max=[(height[i,:]<toplevel).argmax() for i in range(ndata)]
     ih=int(np.mean(ih)); ih_min=int(np.mean(ih_min)); ih_max=int(np.mean(ih_max))
     for i in range(ndata):
      X=np.asarray([height[i,ih_max],height[i,ih_min]])[:,np.newaxis]
      Y=[pressure_o[i,ih_max],pressure_o[i,ih_min]]
      reg = LinearRegression().fit(X, Y)
      pressure_98[i]=reg.intercept_+toplevel*reg.coef_
      Y=[pcol[i,ih_max],pcol[i,ih_min]]
      reg = LinearRegression().fit(X, Y)
      pcol_98[i]=reg.intercept_+toplevel*reg.coef_
   #Insert the level 9.8 
   if remove_strat:
    pressure_o=np.insert(pressure_o,ih_min+1,pressure_98,axis=1)
    pcol=np.insert(pcol,ih_min+1,pcol_98,axis=1)
   ntime=len(date)
   nz=len(pressure_o[0,:])
   ph=np.zeros((ntime,nz-1))
   ph=(pressure_o[:,:-1]+pressure_o[:,1:])/2
   h=create_h(pressure_o)
   if remove_strat:
    xocs_o = h*pcol/ (p_surface-ph[:,ih_min+1])[:,np.newaxis]
    xocs_o=np.sum(xocs_o[:,ih_min+1:],axis=1)
   else:
    xocs_o = h*pcol/ (p_surface-ph[:,0])[:,np.newaxis]
    xocs_o=np.sum(xocs_o,axis=1)

   if (punit == "ppmv") & (stat != 'izana') :  xocs_o=xocs_o*10**6
   if (punit == "ppbv") | (stat == 'izana'):  xocs_o=xocs_o*10**3
   if stat == "zugspitze": xocs_o*=10**(-3)
   for iff,ff in enumerate(list_files):
    if "XCOS" in ff: continue
    if (not 'nc' in ff[-2:]): continue
    if  not "COS" in ff:continue
    if "Diurn" in ff: continue
   # if not 'OPT-JM' in ff: continue
    print(ff,FILE_NAME)
    xcos={"obs":[],"cos":[],"lat":[],"lon":[],"station":[],"date":[]}
    Data=xr.open_dataset(dir_transcom+"/"+model+"/TIER3/"+ff,decode_times=True).to_dataframe().reset_index()
    name_var="cos"
    if not "station" in  Data.columns: 
     #del Data["station"]
     Data=Data.rename(columns={"id":"station"})
    Data["time"]=[datetime.datetime(tt.year,tt.month,tt.day,tt.hour,0,0) for tt in Data.time]
    Data=Data[(Data.station==list_station[iss])&(Data.time.dt.year>=begy)&(Data.time.dt.year<=endy)]
    if Data.empty: continue
    Data=Data.groupby(["time","nlev"]).mean().reset_index()
    Data_tot=pd.DataFrame()
    for kk in Data.nlev.unique():
     Data_kk=Data[Data.nlev==kk].set_index("time").copy(deep=True)
     date=[dd for dd in date if dd in Data_kk.index]
     Data_kk=Data_kk.loc[date]
     Data_tot=Data_tot.append(Data_kk.reset_index())
    for  itt,tt in enumerate(date):
     Data_t=Data_tot[Data_tot.time==tt].copy(deep=True)
     Data_t=Data_t.groupby(["time","nlev"]).mean().reset_index()
     pressure_m=np.copy(Data_t.pres.values)/100.
     cos_m=np.copy(Data_t.cos.values)
     if pressure_m[0]>pressure_m[1]:
      pressure_m=pressure_m[::-1]
      cos_m=cos_m[::-1]
     #Interpolation to the NDAC retrieval grid
     f= interp1d(pressure_m,cos_m, kind='linear',fill_value='extrapolate')
     cos_m=f(pressure_o[itt,:])
     xocs_m = h[itt,:]*cos_m/ (p_surface[itt]-ph[itt,ih-1])
     xocs_m=  np.sum(xocs_m[ih:])
     xcos["date"].append(tt)
     xcos["station"].append(list_station[iss])
     xcos["lon"].append(longitude)
     xcos["lat"].append(latitude)
     xcos["cos"].append(xocs_m)
     xcos["obs"].append(xocs_o[itt])
    xcos=pd.DataFrame(xcos)
    name_file=dir_transcom+"/"+model+"/TIER3/"+ff+FILE_NAME
    os.system("rm -f "+name_file)
    xcos.to_xarray().to_netcdf(name_file)

list_flux=["COS_COS-ANTHR-Zumkehr-Month","COS_COS-BB-UU-Month","COS_COS-OCEAN-LennartzPICSES-Month","COS_COS-BB-Stin-Month","COS_COS-OCEAN-Lennartz-Month","COS_COS-OPT-LSCE-Month","COS_COS-OPT-JM-Month","COS_COS-VEG-SIB4-Month","COS_COS-VEG-ORC-Month","COS_COS-SOIL-SIB4-Month","COS_COS-SOIL-ORC-Month"]
#list_station_xocs=["paramaribo","kiruna","bremen","ny.alesund","zugspitze","izana"]
#list_station=     ['PAR',       'KIR',  'BRE',   'NY.',       'ZUG','IZA']

dir_model=dir_transcom+"/"+model+"/TIER3/"
for ff in list_flux:
 if remove_strat:
  os.system("rm  -f "+dir_model+ff+"_name_station"+"_XCOS.nc")
 else:
  os.system("rm  -f "+dir_model+ff+"_name_station"+"_XCOSt.nc")
 data=pd.DataFrame()
 for ss in  list_station_xocs:
  listfile=os.listdir(dir_model)
  listfile=[files for files in listfile if (ss in files)&(ff in files) ]
  for iff,files in enumerate(listfile):
   tmp=xr.open_dataset(dir_model+"/"+files).to_dataframe()
   data=data.append(tmp)
#   os.system('rm  -f '+dir_model+"/"+files)
  data=data.sort_values("date")
  data.set_index("date",drop=True,inplace=True) 
  #Remove erroneous values
  data=data[data.obs>100]
  data['z_score']=stats.zscore(data['obs'])
  data=data[data.obs>300*10**(-6)]
  data=data[data.z_score.abs()<3]
  print(dir_model+ff+ss+"_XCOS.nc")
 if remove_strat:
  output_name=dir_model+ff+"_name_station"+"_XCOS.nc"
 else:
  output_name=dir_model+ff+"_name_station"+"_XCOSt.nc"
 data.to_xarray().to_netcdf(output_name) 

