#From https://ecowatch.ncddc.noaa.gov/thredds/dodsC/ncom/ncom_reg1_agg/NCOM_Region_1_Aggregation_best.ncd.html
#Convert the .nc files of HYCOM to .npz and preprocess the data

import datetime as dt
from datetime import timedelta  
import numpy as np
from scipy.io import netcdf


f = netcdf.netcdf_file("P1.nc", 'r')

water_v=np.copy(f.variables["water_v"].data*0.001) #0.001 is the scale of the speed data we are converting to m/s
water_u=np.copy(f.variables["water_u"].data*0.001)
Time=np.copy(f.variables["time"].data)
Lat=np.copy(f.variables["lat"].data)
Long=np.copy(f.variables["lon"].data)

f.close()

for i in range(2,14,1):
    f = netcdf.netcdf_file("P"+str(i)+".nc", 'r')
    water_v=np.concatenate((water_v, np.copy(f.variables["water_v"].data*0.001)),axis=0) #0.001 is the scale of the data
    water_u=np.concatenate((water_u, np.copy(f.variables["water_u"].data*0.001)), axis=0)
    Time=np.concatenate((Time, np.copy(f.variables["time"].data)), axis=0)
    f.close()


#Defined based on the data we donwloaded from HYCOM website
DateStart=dt.datetime(2009, 1, 1)
DateEnd=dt.datetime(2014, 12, 1,12)

LatLong=[]
Speed_Temp=[]
for IdxLat in range(len(Lat)):
    for IdxLong in range(len(Long)):
        if water_v[1, 0, IdxLat, IdxLong]!=-30: #Missing Value 
            Speed_Temp.append((water_v[:, 0, IdxLat, IdxLong]**2+water_u[:, 0, IdxLat, IdxLong]**2)**0.5)
            LatLong.append([Lat[IdxLat], Long[IdxLong]])

Speed_Temp=np.array(Speed_Temp)
LatLong=np.array(LatLong)


Speed=[]
DateTime=[]
NumIndex=(DateEnd-DateStart).days*8#number of hours in the interval

HistoricalDateIdx=-1
NuberOfSubstitutions=0


#Fill missing days and other ordinary missing values
#Attention in the hycom Time fromat:
#dt.datetime(2000,1,1,0) + timedelta(hours=int(Time[0]))=dt.datetime(2009,1,1,0)

for i in range(NumIndex):
    
    DateTime.append(DateStart+dt.timedelta(hours=3*i))
    
    HistoricalDateIdx=HistoricalDateIdx+1
    
    if (dt.datetime(2000,1,1,0) + timedelta(hours=int(Time[HistoricalDateIdx])))==(DateStart+dt.timedelta(hours=3*i)):
        Speed.append(Speed_Temp[:,HistoricalDateIdx])
    
    #If idx does not exist represent it is as the previous day, we used this approach because there are very few missing data
    #There are 246 missing data points in a set of 17280 elements
    else:
        HistoricalDateIdx=HistoricalDateIdx-1
        Speed.append(Speed_Temp[:,HistoricalDateIdx])
        NuberOfSubstitutions=NuberOfSubstitutions+1


Speed=np.array(Speed)
Speed=Speed.T

DateTime=np.array(DateTime)
SpeedByDay_AllLatLong=[]

#Organize the current speed data in a daily format
for TargetIdxLatLong in range(Speed.shape[0]):
    CountHour_Day=0
    Temp_SpeedByDay=np.empty((8))
    Temp_SpeedByDay[:]=np.NaN  
    
    SpeedByDay=[]
    
    for IdxHour in range(Speed.shape[1]):
        
        Temp_SpeedByDay[CountHour_Day]=Speed[TargetIdxLatLong, IdxHour]
        
        CountHour_Day=CountHour_Day+1
        if CountHour_Day==8:
            CountHour_Day=0
            SpeedByDay.append(Temp_SpeedByDay)
            Temp_SpeedByDay=np.empty((8))
            Temp_SpeedByDay[:]=np.NaN          
    
    SpeedByDay=np.array(SpeedByDay)
    SpeedByDay_AllLatLong.append(SpeedByDay)

SpeedByDay_AllLatLong=np.array(SpeedByDay_AllLatLong)#LatLong, Day, Time
np.savez("CurrentSpeedHYCOM2009_2014.npz",Speed=Speed, SpeedByDay_AllLatLong=SpeedByDay_AllLatLong, DateTime=DateTime,LatLong=LatLong  )