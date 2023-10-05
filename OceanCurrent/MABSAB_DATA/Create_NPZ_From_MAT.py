#Organize data from MABSAB- reduce its size

import numpy as np
import datetime as dt
import csv
from scipy.io import loadmat
from ReadTurbineData import ReadTurbineData
from ConcatenateDataForMultipleYears import ConcatenateData


YearsAnalysed=[2009,2010,2011,2012,2013,2014]

#Get coastline data
CoastLine=[]

File_CoastLine=open('GEO_data/Coastline_NC.csv', "r")
File_CoastLine_csv=csv.reader(File_CoastLine,delimiter=',')

for EachLine in File_CoastLine_csv:

    if File_CoastLine_csv.line_num > 1:
        CoastLine.append([float(EachLine[1]), float(EachLine[0])] ) #LatLong

CoastLine=np.array(CoastLine)


#Distance between two lat long points
def DistanceToShore (CoastLine, LatLong1): #Compute distance to shore in km of a lat long point
   
    CoastLine=CoastLine*2*np.pi/360
    LatLong1=LatLong1*2*np.pi/360
    
    LatLong1=np.reshape(LatLong1,(1,2))
    dLat=LatLong1[:,0]-CoastLine[:,0]
    dLong=LatLong1[:,1]-CoastLine[:,1]
    
    a=np.power(np.sin(dLat/2),2)+np.cos(CoastLine[:,0])*np.cos(LatLong1[:,0])*np.power(np.sin(dLong/2),2)
    c=2*np.arcsin(np.minimum(1,np.sqrt(a)))
    d=6367*c
    
    Distance=np.min(d) #Minimum distance to shore km
    
    return Distance


#Set of years we are interested in preprocess the data
for YearOfSimulation in YearsAnalysed:
    print('Year: '+str(YearOfSimulation))
 
    FileSaveName='CurrentSpeedMABSAB'
    
    TurbineFile='RM4Sandia.txt'
    
    Turbine=ReadTurbineData(TurbineFile)
    RotorDepth=Turbine["RotorDepth"]
    MinDepth=Turbine["MinDepth"]
    MaxDepth=Turbine["MaxDepth"]
    RatedPower=Turbine["RatePower"]
    
    DepthDomain= loadmat('depth_domain.mat')
    SiteDepth=DepthDomain["h"]
    DepthLayers=DepthDomain["depth_rho0"]
    
    #Get the feasible site locations based on the maximum and minimum depths
    IndexFDepth=((SiteDepth<=MaxDepth) & (SiteDepth>=MinDepth))
    IndexFDepth=(SiteDepth>=MinDepth) 
    
    #Get index closest to the rated depth of operation for the turbine (RotorDepth)
    IndexTurbineDepth=(np.argmin(np.abs(DepthLayers-RotorDepth),axis=2))[IndexFDepth]
    
    CurrentRawData= loadmat("MABSAB_"+str(YearOfSimulation)+'.mat')
    
    ocean_time=CurrentRawData['ocean_time']
    
    #Get Time 
    OceanDateTime=[]
    for i in range(ocean_time.shape[0]):
        Date=dt.timedelta(seconds=ocean_time[i,0]) + dt.datetime(1858,11,17,0,0,0)
        if Date.year==YearOfSimulation:
            OceanDateTime.append(Date)
    
    lon_range=CurrentRawData['lon_range']
    lat_range=CurrentRawData['lat_range']
    
    #Ocean current speed filtered for feasible site locations
    #We ignore the last day in the speed data because it corresponds to the first day of the next year
    udata=CurrentRawData['udata'][IndexFDepth,:,0:len(OceanDateTime)]#speed towards east (Zonal velocity)
    vdata=CurrentRawData['vdata'][IndexFDepth,:,0:len(OceanDateTime)]#speed towards north (Meridional velocity)
    
    #Filter speed to get only the depth closest to the rated depth
    u_filtered=np.zeros((udata.shape[0],udata.shape[2]),dtype=float)
    v_filtered=np.zeros((udata.shape[0],udata.shape[2]),dtype=float)
    
    for site in range(udata.shape[0]):
        u_filtered[site,:]=udata[site,IndexTurbineDepth[site],:]
        v_filtered[site,:]=vdata[site,IndexTurbineDepth[site],:]
    
    #Get the magnitude of the speed vector
    Speed=np.sqrt(u_filtered**2+v_filtered**2)
    
    #Release RAM memory
    del(CurrentRawData, udata, vdata, u_filtered, v_filtered)
    
    DimLat=lat_range.shape[1]
    DimLong=lon_range.shape[0]
    
    #Organize the data 
    lon_range=np.repeat(lon_range,DimLat,axis=1).T
    lat_range=np.repeat(lat_range,DimLong,axis=0).T
    
    LatLong=np.zeros((np.sum(IndexFDepth),2),dtype=float)
    
    LatLong[:,0]=lat_range[IndexFDepth]
    LatLong[:,1]=lon_range[IndexFDepth]
    

    DepthSites=SiteDepth[IndexFDepth]
    
    ShoreDistance=[]
    for LatLong1 in LatLong:
        Distance=DistanceToShore(CoastLine, LatLong1)
        ShoreDistance.append(Distance)

    ShoreDistance=np.array(ShoreDistance)
    
    ReadMe='\
    CurrentSpeed: Ocean current speed\
    LatLong: Latitude,Logitude data\n\
    1) The data is in Daily discretization starting at 1/1/'+str(YearOfSimulation)+' and going up to \
    12/31/'+str(YearOfSimulation)
    
    np.savez(FileSaveName+str(YearOfSimulation)+'.npz',ReadMe=ReadMe, CurrentSpeed=Speed, DepthSites=DepthSites, ShoreDistance=ShoreDistance, RatedPower=RatedPower,LatLong=LatLong, OceanDateTime=OceanDateTime)
    
    
ConcatenateData(YearsAnalysed, FileSaveName, RotorDepth) #Concatenate multiple years of data in a single .npz file