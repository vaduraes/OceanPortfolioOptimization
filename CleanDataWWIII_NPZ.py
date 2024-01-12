#Preprocess the data from WWIII
#Filter data to cover only the depth region that the turbine can operate
#Fix missing data

import numpy as np
import xarray as xr
import csv


WWIII_Data=np.load("WWIII_RawData.npz")
Depth_NETCDF = xr.open_dataset("./Depths.nc")

#missing values are equal to 9999
LatLong=WWIII_Data["LatLong"]
LatLong[:,1]=LatLong[:,1]-360

Hs=WWIII_Data["Hs"]
Tp=WWIII_Data["Tp"]

##Constrain the data to 2009-2013. Use if the raw data is between 2009-2015
#Hs=Hs[:,0:14607]
#Tp=Tp[:,0:14607]

#Compute depth in a given location
def GetDepth(Lat,Long):
    
    I_lat=np.argmin(np.square(Depth_NETCDF.lat.data-Lat))
    I_lon=np.argmin(np.square(Depth_NETCDF.lon.data-Long))
    
    distError=(np.square(Depth_NETCDF.lat.data[I_lat]-Lat)+np.square(Depth_NETCDF.lon.data[I_lon]-Long))**0.5
    
    depth=Depth_NETCDF.elevation.data[I_lat,I_lon]  
    
    return depth, distError

DepthAllSites=[]
for i in range(LatLong.shape[0]):
    Depth, Error=GetDepth(LatLong[i,0],LatLong[i,1])#Get depth for each .mat site location
    
    DepthAllSites.append(Depth)
    
DepthAllSites=np.array(DepthAllSites)

#IdxFeasibleSites=(DepthAllSites<=-40)*(DepthAllSites>=-100)
IdxFeasibleSites=(DepthAllSites<=-1)*(DepthAllSites>=-200)

LatLong=LatLong[IdxFeasibleSites,:]
Tp=Tp[IdxFeasibleSites,:]
Hs=Hs[IdxFeasibleSites,:]

#Fill the missing values with the value in the previous 3h
#We used this simple procedure because there are very few missing values
for i in range(Tp.shape[0]):
    for j in range(Tp.shape[1]):
        if Tp[i,j]==9999:
            Tp[i,j]=Tp[i,j-1]
            
        if Hs[i,j]==9999:
            Hs[i,j]=Hs[i,j-1]

WWIII_Data.close()

#Get coastline data
CoastLine=[]

File_CoastLine=open('./GEO_data/Coastline_NC.csv', "r")
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

#Minimum distance from shore
ShoreDistance=[]
for LatLong1 in LatLong:
    Distance=DistanceToShore(CoastLine, LatLong1)
    ShoreDistance.append(Distance)

ShoreDistance=np.array(ShoreDistance)
    
np.savez('WaveHsTp_WWIII.npz', HsNC=Hs, TpNC=Tp, LatLong=LatLong, Depth=(DepthAllSites[IdxFeasibleSites]*-1), ShoreDistance=ShoreDistance)