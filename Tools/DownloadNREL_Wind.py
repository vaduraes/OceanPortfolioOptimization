#Functions to download NREL wind data
#run on geo_env
#Download wind data from NREL
#Run on geo_env
import h5pyd
import numpy as np
import xarray as xr
import time

from GeneralGeoTools import GetDepth

def DonwloadNREL_WindData(InputDataPath, SavePath, Data2Download="windspeed_100m", LatMinMax=(33, 37), LongMinMax=(-81,-73), DepthMinMax=(0,1000)):
    #LatMinMax: Latitude Min and Max for download
    #LongMinMax: Longitude Min and Max for download


    #Hourly Wind Data from NREL
    #Information available
    ListOfAvailableData=["winddirection_100m","winddirection_10m","winddirection_120m","winddirection_140m","winddirection_160m",
                        "winddirection_200m","winddirection_40m","winddirection_60m","winddirection_80m",
                        "windspeed_100m","windspeed_10m","windspeed_120m","windspeed_140m","windspeed_160m","windspeed_200m",
                        "windspeed_40m","windspeed_60m","windspeed_80m"]

    if Data2Download not in ListOfAvailableData:
        print("Data not available")
        return


    NREL = h5pyd.File("/nrel/wtk-us.h5", 'r') #Open conection with NREL server


    #The original projection on NREL data is a modified Lambert Conic
    Coordinates=NREL["coordinates"][:,:]#Coordinates in latitude longitude for each point y,x of the original data

    XY_NC=[]#coordinates we are interested in downloading
    #Get coordinates we are interested in downloading. 
    #We group these coordinates such that we can decrease the total size of the file downloaded

    NREL_LatLong=[]#Coordinates in latitude longitude for each point y,x of the original data
    for X in range(Coordinates.shape[0]):
        for Y in range(Coordinates.shape[1]):
            NREL_LatLong.append([X,Y,Coordinates[X,Y]['lat'],Coordinates[X,Y]['lon']])
    NREL_LatLong=np.array(NREL_LatLong)

    IdxIn = (NREL_LatLong[:,2]>=LatMinMax[0]) * (NREL_LatLong[:,2]<=LatMinMax[1]) \
            *(NREL_LatLong[:,3]>=LongMinMax[0]) * (NREL_LatLong[:,3]<=LongMinMax[1])
    NREL_LatLong=NREL_LatLong[IdxIn,:]
    Depths=GetDepth(InputDataPath, NREL_LatLong[:,2:4])

    IdxIn=(Depths>=DepthMinMax[0]) * (Depths<=DepthMinMax[1])
    NREL_LatLong=NREL_LatLong[IdxIn,:]

    UniqueX=np.unique(NREL_LatLong[:,0])


    MaxGroupSize=100
    for X in UniqueX:
        YValues=NREL_LatLong[NREL_LatLong[:,0]==X,1]       
        YMin=np.min(YValues)
        YMax=np.max(YValues)+1
        
        #Divide the Ymin and Ymax in groups of MaxGroupSize
        Steps=np.arange(YMin,YMax,MaxGroupSize)
        if Steps[-1]<YMax:
            Steps=np.append(Steps,YMax)

        for i in range(len(Steps)-1):
            XY_NC.append([int(X),int(Steps[i]),int(Steps[i+1])])

    XY_NC=np.array(XY_NC,dtype=int)

    # #Downloading the data
    print("Downloading NREL data\n")
    i=-1
    while i!=(len(XY_NC)-1):
        error=0
        i=i+1
        try:
            if i==0:
                #Create initial windspeed matrix and concatenate future wind data (windspeedTemp) on this same matrix
                windspeed=NREL[Data2Download][:,XY_NC[i,0],XY_NC[i,1]:XY_NC[i,2]]
                
                if windspeed.ndim==1:
                    windspeed=np.reshape(windspeed,(windspeed.shape[0],1))
                
            
            else:
                windspeedTemp=NREL[Data2Download][:,XY_NC[i,0],XY_NC[i,1]:XY_NC[i,2]]
                if windspeedTemp.ndim==1:
                    windspeedTemp=np.reshape(windspeedTemp,(windspeedTemp.shape[0],1))
                
            
            Lat=Coordinates['lat'][XY_NC[i,0],XY_NC[i,1]:XY_NC[i,2]]
            Lat=np.reshape(Lat,(len(Lat),1))
            
            Long=Coordinates['lon'][XY_NC[i,0],XY_NC[i,1]:XY_NC[i,2]]
            Long=np.reshape(Long,(len(Long),1))

        except:
            error=1
            i=i-1
            print("Error in the NREL server, too many requests- Waiting access release (1h Waiting))")
            time.sleep(61*60)
            NREL = h5pyd.File("/nrel/wtk-us.h5", 'r') #Open conection with NREL server
            break
            
            
        if error==0:      
            if i==0:
                LatLong=np.concatenate((Lat,Long),axis=1)      
            else:
                LatLong=np.concatenate((LatLong,np.concatenate((Lat,Long),axis=1)),axis=0)
                windspeed=np.concatenate((windspeed,windspeedTemp),axis=1)

        print('Download----- %.1f%%'% ((i+1)/len(XY_NC)*100))
        
        
    ReadMe='\
    windspeed: m³\s \n\
    LatLong: Latitude,Logitude data\n\
    1) The data is in hourly discretization starting at 1/1/2007 and going up to\
    12/31/2013 23:00'

    np.savez(SavePath, ReadMe=ReadMe, windspeed=windspeed, LatLong=LatLong)
