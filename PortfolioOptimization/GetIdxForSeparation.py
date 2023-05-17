import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt

def DistanceBetweenLatLong(LatLong1, LatLong2):
    LatLong1=LatLong1*2*np.pi/360
    LatLong2=LatLong2*2*np.pi/360
    
    dLat=np.reshape(LatLong1[:,0],(len(LatLong1[:,0]),1))-np.reshape(LatLong2[:,0],(1,len(LatLong2[:,0])))
    dLong=np.reshape(LatLong1[:,1],(len(LatLong1[:,1]),1))-np.reshape(LatLong2[:,1],(1,len(LatLong2[:,1])))
    
    P1=np.reshape(np.cos(LatLong1[:,0]),(LatLong1.shape[0],1))
    P2=np.reshape(np.cos(LatLong2[:,0]),(1,LatLong2.shape[0]))
    
    a=np.power(np.sin(dLat/2),2) + (P1*P2)*np.power(np.sin(dLong/2),2)
    c=2*np.arcsin(np.minimum(1,np.sqrt(a)))
    Distance=6367*c #[km]
    
    return Distance

def GetWindIDXs(Radious, WindLatLong):

    NumWindSites=WindLatLong.shape[0]
    
    Set1IDX=np.arange(0,int(NumWindSites/2),1)
    Set2IDX=np.arange(int(NumWindSites/2),NumWindSites,1)
    
    LatLongSet1=WindLatLong[Set1IDX,:]
    LatLongSet2=WindLatLong[Set2IDX,:]
    
    Distance=DistanceBetweenLatLong(LatLongSet1, LatLongSet2)
    
    MinDistance2_Set1=np.min(Distance,axis=1)
    MinDistance2_Set2=np.min(Distance,axis=0)
    
    Set3IDX=np.concatenate((np.where(MinDistance2_Set1<Radious)[0], int(NumWindSites/2)+np.where(MinDistance2_Set2<Radious)[0]))

    return Set1IDX, Set2IDX, Set3IDX

def GetOceanIDXs(Radious, OceanLatLong):

    NumOceanSites=OceanLatLong.shape[0]
    
    Set1IDX=np.arange(0,int(NumOceanSites/2),1)
    Set2IDX=np.arange(int(NumOceanSites/2),NumOceanSites,1)
    
    LatLongSet1=OceanLatLong[Set1IDX,:]
    LatLongSet2=OceanLatLong[Set2IDX,:]
    
    Distance=DistanceBetweenLatLong(LatLongSet1, LatLongSet2)
    
    MinDistance2_Set1=np.min(Distance,axis=1)
    MinDistance2_Set2=np.min(Distance,axis=0)
    
    Set3IDX=np.concatenate((np.where(MinDistance2_Set1<Radious)[0], int(NumOceanSites/2)+np.where(MinDistance2_Set2<Radious)[0]))

    return Set1IDX, Set2IDX, Set3IDX



#min_longitude=-77.1
#max_longitude=-74
#
#min_latitude=32.8
#max_latitude=36.1
#
#xlim =[min_longitude,max_longitude]
#ylim=[min_latitude, max_latitude]
#
#df = gpd.read_file(ShapeFileCoast)
#df1 = gpd.read_file(ShapeFileStates)
#
#fig, ax = plt.subplots(figsize  = None)
#
#df.plot(color='black',linewidth=1,ax=ax)
#df1.plot(color='black',linewidth=1,ax=ax)
#
#plt.scatter(OceanLatLong[:,1],OceanLatLong[:,0],c='r',s=2)#Wind
#
#plt.scatter(OceanLatLong[Set3IDX,1], OceanLatLong[Set3IDX,0],s=5, marker='+',c='k',linewidth=0.3)#Wind
#
#ax.set_xlim(xlim)
#ax.set_ylim(ylim)
#plt.xlabel("Longitude")
#plt.ylabel("Latitude")
#
#plt.title("Optimal Wind Energy Portfolio, LCOE=")

    