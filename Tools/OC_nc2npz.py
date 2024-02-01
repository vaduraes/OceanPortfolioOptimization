#py37 env
import numpy as np
from datetime import datetime, timedelta
import sys
import numpy as np
import datetime as dt
from dateutil.relativedelta import relativedelta
import os
import urllib.request
from multiprocessing import Pool
from tqdm import tqdm  
import datetime as dt
from datetime import timedelta  
import numpy as np
from scipy.io import netcdf_file

def ReadHYCOMData(SavePath=None):
    PathInputData="./InputData/OceanCurrent/Hycom/"
    ReferenceTime=datetime(2000, 1, 1, 0, 0, 0) # According to HYCOM

    ListFiles=os.listdir(PathInputData) #Get list of all files in the directory
    f = netcdf_file(PathInputData+ListFiles[0], 'r',mmap=False)
    lat=f.variables['lat'].data
    lon=f.variables['lon'].data
    depth=f.variables['depth'].data
    f.close()

    LatLong=np.array(np.meshgrid(lon,lat)).reshape(2, -1).T[:,[1,0]]

    DateTimeS1=[]
    DateTimeS2=[]
    for files in ListFiles:
        DateTimeS1.append(datetime.strptime((files.split("_")[0]), "%Y-%m-%dT%H"))
        DateTimeS2.append(datetime.strptime((files.split("_")[1].split(".")[0]), "%Y-%m-%dT%H"))
        
    DateTimeS1=np.array(DateTimeS1)
    DateTimeS2=np.array(DateTimeS2)


    EndDate=np.max(DateTimeS2)
    StartDate=np.min(DateTimeS1)
    NumTimeSteps=(int((EndDate-StartDate).total_seconds()/3600))/3+1
    TimeList=[StartDate + timedelta(hours=3*x) for x in range(0, int(NumTimeSteps))]

    U_speed=np.ones((len(TimeList), len(depth), int(len(lat)*len(lon))))*-30000
    V_speed=np.ones((len(TimeList), len(depth), int(len(lat)*len(lon))))*-30000
    U_speed=U_speed.astype('float32')
    V_speed=V_speed.astype('float32')

    for file in tqdm(ListFiles):
        f = netcdf_file(PathInputData+file, 'r',mmap=False)
        water_u=(f.variables['water_u'].data*0.001).astype('float32') #0.001 is the scale of the speed data we are converting to m/s
        water_v=(f.variables['water_v'].data*0.001).astype('float32') #0.001 is the scale of the speed data we are converting to m/s
        time=f.variables['time'].data
        TimeListTmp=[ReferenceTime + timedelta(hours=x) for x in time]
        IdxTimeList=[int((time-StartDate).total_seconds()/3600/3) for time in TimeListTmp]
        
        U_speed[IdxTimeList,:,:]=np.reshape(water_u, (len(IdxTimeList), len(depth), int(len(lat)*len(lon))))
        V_speed[IdxTimeList,:,:]=np.reshape(water_v, (len(IdxTimeList), len(depth), int(len(lat)*len(lon)))) 
        
        f.close()

    OCSpeed=(U_speed**2+V_speed**2)**0.5
    OCSpeed[OCSpeed>10]=-1 #missing values are set to -1
    
    #Weep sites that have at least one data point not missing
    IdxSitesWithData=~((np.sum(OCSpeed==-1,axis=0)==OCSpeed.shape[0]).sum(axis=0)==OCSpeed.shape[1]) #Sites with at least one data point
    LatLong=LatLong[IdxSitesWithData,:]
    OCSpeed=OCSpeed[:,:,IdxSitesWithData]
    
    if SavePath!=None:
        np.savez(SavePath, OCSpeed=OCSpeed, LatLong=LatLong, depth=depth, TimeList=TimeList)
    else:  
        return OCSpeed, LatLong, depth, TimeList