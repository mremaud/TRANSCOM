from calendar import isleap
import sys
import string, os
import numpy, math
from scipy.io import netcdf
from netCDF4 import Dataset
import pandas as pd
from pandas import *
import datetime
#flag_sel=0: donnee de lapres midi
year_ref=2006
month_ref=4
import numpy as np

def ccgv(data,seas):
     data.reset_index(inplace=True)
     data['frac']=data.apply(lambda row: fraction_an(row),axis=1)
     file_out='station.txt'
     os.system('rm -f '+file_out)
     np.savetxt(file_out,data[['frac','cos']],fmt='%4.8f %3.3f')
     file_fit='station2.txt'
     #Model: without cal and ori
     os.system('./ccgcrv -cal -interv 7 -short 80 -equal -smooth -npoly 2 -nharm 4  -func -trend  -f  '+file_fit+' '+file_out)
     cgv=np.loadtxt(file_fit)
     data2=pd.DataFrame()
     if seas: 
      data2['meas']=np.copy(cgv[:,4] -cgv[:,5])
     else:
      data2['meas']=np.copy(cgv[:,4])
     data2['date']=[datetime.datetime(int(cgv[ii,0]),int(cgv[ii,1]),int(cgv[ii,2])) for ii in range(len(cgv))]
     return data2

def fraction_an(row):
        hour=0
        date=row['time']
        day_in_year = isleap(date.year) and 366 or 365
        mon_in_year=isleap(date.year) and [31,29,31,30,31,30,31,31,30,31,30,31] or [31,28,31,30,31,30,31,31,30,31,30,31]
        start_date=datetime.datetime(date.year,1,1,0,0)
        end_date = datetime.datetime(date.year,date.month,date.day,date.hour,0,0)
        result=date.year+(end_date-start_date).total_seconds()/(day_in_year*86400.)
        return result


def numbers_to_strings(argument):
    switcher = {
        "COS-VEG-ORC-Month"     : "Exp1a",
        "COS-VEG-SIB4-Month"  : "Exp1b",
        "COS-VEG-ORC-Diurn"    : "Exp1c",
        "COS-VEG-SIB4-Diurn" : "Exp1d",
        "COS-BB-Stin" : "Exp2a",
        "COS-BB-UU" : "Exp2b",
        "COS-ANTHR-Zumkehr-Month" : "Exp3a",
        "COS-SOIL-ORC-Month": "Exp4a",
        "COS-SOIL-SIB4-Month": "Exp4b",
        "COS-SOIL-oxicLaunois_clim": "Exp4c",
        "COS-SOIL-ORC-Diurn": "Exp4e",
        "COS-SOIL-SIB4-Diurn": "Exp4d",
        "COS-OCEAN-Lennartz-Month": "Exp5a",
        "COS-OCEAN-LennartzPICSES-Month": "Exp5b",
        "COS-OPT-JM-Month": "Exp6a",
        "COS-OPT-LSCE-Month": "Exp6b",
        "CO2-GPP-ORC-Month"     : "Exp11a",
        "CO2-GPP-SIB4-Month"  : "Exp11b",
        "CO2-TER-ORC-Month"    : "Exp11c",
        "CO2-TER-SIB4-Month" : "Exp11d",
        "CO2-GPP-ORC-Diurn"     : "Exp11e",
        "CO2-GPP-SIB4-Diurn"  : "Exp11f",
        "CO2-TER-ORC-Diurn"    : "Exp11g",
        "CO2-TER-SIB4-Diurn" : "Exp11h",

    }
    return switcher.get(argument, "nothing")

def numbers_to_strings_inv(argument):
    switcher = {
        "Exp1a" : "COS-VEG-ORC-Month",
        "Exp1b" : "COS-VEG-SIB4-Month",
        "Exp1c" : "COS-VEG-ORC-Diurn",
        "Exp1d"  : "COS-VEG-SIB4-Diurn",
        "Exp2a" : "COS-BB-Stin",
        "Exp2b" : "COS-BB-UU",
        "Exp3a" : "COS-ANTHR-Zumkehr",
        "Exp4a" : "COS-SOIL-ORC-Month",
        "Exp4b" : "COS-SOIL-SIB4-Month",
        "Exp4c" : "COS-SOIL-oxicLaunois_clim",
        "Exp4e" : "COS-SOIL-ORC-Diurn",
        "Exp4d" : "COS-SOIL-SIB4-Diurn",
        "Exp5a" : "COS-OCEAN-Lennartz-Month",
        "Exp5b" : "COS-OCEAN-LennartzPICSES-Month",
        "Exp6a" : "COS-OPT-JM-Month",
        "Exp6b" : "COS-OPT-LSCE-Month",
        "Exp11a" : "CO2-GPP-ORC-Month",
        "Exp11b" : "CO2-GPP-SIB4-Month",
        "Exp11c" : "CO2-TER-ORC-Month",
        "Exp11d" : "CO2-TER-SIB4-Month",
        "Exp11e" : "CO2-GPP-ORC-Diurn",
        "Exp11f" :"CO2-GPP-SIB4-Diurn" ,
        "Exp11g" : "CO2-TER-ORC-Diurn",
        "Exp11h" : "CO2-TER-SIB4-Diurn",
    }
    return switcher.get(argument, "nothing")



def psol():
  year_ref=2008
  month_ref=6
  gasct = 8.314    # J K-1 mol-1
  mmol  = 0.02894  # kg mol -1
  rg = 9.80665     # m s-2
  dry_mass= 28.966/1000.    # dry air molar mass, kg/mol

  # surface orography
  DIRRESTART = '/home/satellites1/fcheval/LMDZ5/'
  sgeop = DIRRESTART + 'phystoke.an'+str(year_ref)+'.m'+"%2.2i"%month_ref+'.nc'
  geopfile=netcdf.netcdf_file(sgeop,'r')
  zsurf = geopfile.variables['phis'][:]
  tmp = geopfile.variables['t'][:]
  nlon = len( tmp[0,0,0,:] ) + 1  # file phystoke does not contain twin grid point
  nlat = len( tmp[0,0,:,0] )
  nlev = len( tmp[0,:,0,0] )
  t = numpy.zeros( (nlev,nlat,nlon) )
  t[:,:,:nlon-1] = tmp[0,:,:,:] # use first time step
  t[:,:,nlon-1] = t[:,:,0]
  soro = numpy.zeros( (nlat,nlon) )
  soro[:,:nlon-1] = zsurf / rg
  soro[:,nlon-1] = soro[:,0]
  geopfile.close()
  thefile = DIRRESTART + 'fluxstoke.an'+str(year_ref)+'.m'+"%2.2i"%month_ref+'.nc'
  f=netcdf.netcdf_file(thefile,'r')
  masse=f.variables['masse'][:]
  masse = masse[0,:,:,:] # use first time step
  aire=f.variables['aire'][:]
  f.close()

  # Compute surface pressure
  smass = numpy.zeros((nlat,nlon))
  for ilev in range(nlev):
    smass[:,:] += masse[ilev,:,:]
  sp = smass[:,:] / aire[:,:] * rg

  ap39 = [ 0., 281.706852929675, 636.375366859432, 1136.4155943717,
  1849.44737296152, 2837.9280405897, 4157.50869711024, 5852.26417553552,
  7944.82760257521, 10419.7321060548, 13199.5402485757, 16116.838631109,
  18892.5506083084, 21142.7861688347, 22446.8782335346, 22497.7595057074,
  21289.4514040892, 19186.4448521733, 16726.0167170181, 14296.6075398689,
  12045.5252156898, 10008.0699708566, 8195.30688158518, 6609.57734809348,
  5246.06241645797, 4093.97069658929, 3137.87110503367, 2359.06747638525,
  1736.92080010302, 1250.0441331629, 877.317766066104, 598.69518472359,
  395.791432392045, 252.262962812769, 154.001062228325, 89.1691757654163,
  48.118336844774, 23.2150838051844, 8.61367138231894, 0. ]
  bp39 = [ 1., 0.988762998564637, 0.974701632524877, 0.955030180542262,
  0.927265686818222, 0.889271183017204, 0.839314287303751,
  0.776176963345792, 0.699355000528228, 0.609380929528209,
  0.508279179863367, 0.400093029287433, 0.291276044521652, 0.1905059741263,
  0.107270901645713, 0.0488076380330498, 0.0163000324948695,
  0.00344452582803985, 0.000364568302239398, 1.32673501351578e-05,
  8.94025963956579e-08, 3.93038386309706e-11, 1.85617183702948e-16,
  3.81418691876566e-25, 9.62385916256612e-40, 4.5206602224536e-65,
  1.46258410571586e-110, 2.18996926640889e-195, 0., 0., 0., 0., 0., 0., 0., 0., 0.,
  0., 0., 0. ]
  ap_lmdz = ap39
  bp_lmdz = bp39

  # LMDZ pressure grid (Pa)
  plmdzf = numpy.zeros((nlev+1,nlat,nlon))
  plmdzh = numpy.zeros((nlev,nlat,nlon))
  for i in range( nlev+1 ) :
     plmdzf[i,:,:] = ap_lmdz[i] + bp_lmdz[i] * sp[:,:]
  plmdzh = 0.5*( plmdzf[0:nlev,:,:] + plmdzf[1:nlev+1,:,:] )
  dp=plmdzf[0:nlev,:,:] - plmdzf[1:nlev+1,:,:]
  return plmdzh,dp

def defheight(clon,clat,alt,h,hsol):
        hh=np.squeeze(h[:,clat,clon])
        hhsol=np.squeeze(hsol[clat,clon])
        altf=max(alt,hhsol)
        altf=(np.abs(hh-altf)).argmin()
        return altf


