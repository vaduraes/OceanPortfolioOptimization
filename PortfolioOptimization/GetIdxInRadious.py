#Get indexes of plants out of a certain range
import numpy as np

#Distance between lat long points
def DistanceBetweenLatLong(LatLong):
    LatLong=LatLong*2*np.pi/360
    LatLong1=np.copy(LatLong)
    
    dLat=np.reshape(LatLong[:,0],(len(LatLong[:,0]),1))-np.reshape(LatLong1[:,0],(1,len(LatLong[:,0])))
    dLong=np.reshape(LatLong[:,1],(len(LatLong[:,1]),1))-np.reshape(LatLong1[:,1],(1,len(LatLong[:,1])))
    
    P1=np.repeat(np.reshape(np.cos(LatLong[:,0]),(LatLong.shape[0],1)),LatLong.shape[0],axis=1)
    P2=np.repeat(np.reshape(np.cos(LatLong1[:,0]),(LatLong1.shape[0],1)),LatLong1.shape[0],axis=1)
    
    a=np.power(np.sin(dLat/2),2) + P1*P2*np.power(np.sin(dLong/2),2)
    c=2*np.arcsin(np.minimum(1,np.sqrt(a)))
    Distance=6367*c #[km]
    
    return Distance

#Get the index of all site locations inside a specific radious (each energy resource is treated separately, we do not consider
#the proximity of different energy resources here)
#Each latlong (latlong1, 2 ,3) means a specific energy resource
def GetIdxInRadious(Radious, LatLong1, LatLong2, LatLong3):
    Distance1=DistanceBetweenLatLong(LatLong1)
    Idx=Distance1<=Radious
    
    IdxR=np.zeros((Idx.shape[0], LatLong2.shape[0]), dtype=bool)
    IdxL=np.zeros((LatLong2.shape[0], Idx.shape[1]), dtype=bool)
    
    Idx1=np.concatenate((Idx,IdxR),axis=1)
    
    Distance2=DistanceBetweenLatLong(LatLong2)
    Idx2=Distance2<=Radious
    
    Idx2=np.concatenate((IdxL,Idx2),axis=1)
    
    Idx=np.concatenate((Idx1,Idx2),axis=0)
    
    IdxR=np.zeros((Idx.shape[0], LatLong3.shape[0]), dtype=bool)
    IdxL=np.zeros((LatLong3.shape[0], Idx.shape[1]), dtype=bool)    
    
    Distance3=DistanceBetweenLatLong(LatLong3)
    Idx3=Distance3<=Radious
    
    Idx1=np.concatenate((Idx,IdxR),axis=1)   
    Idx3=np.concatenate((IdxL,Idx3),axis=1)
    
    Idx=np.concatenate((Idx1,Idx3),axis=0)
    
    return Idx    
	
	
def GetIdxOutRadious(Radious, LatLong1, LatLong2, LatLong3):
    Distance1=DistanceBetweenLatLong(LatLong1)
    Idx=Distance1>=Radious
    
    IdxR=np.zeros((Idx.shape[0], LatLong2.shape[0]), dtype=bool)
    IdxL=np.zeros((LatLong2.shape[0], Idx.shape[1]), dtype=bool)
    
    Idx1=np.concatenate((Idx,IdxR),axis=1)
    
    Distance2=DistanceBetweenLatLong(LatLong2)
    Idx2=Distance2>=Radious
    
    Idx2=np.concatenate((IdxL,Idx2),axis=1)
    
    Idx=np.concatenate((Idx1,Idx2),axis=0)
    
    IdxR=np.zeros((Idx.shape[0], LatLong3.shape[0]), dtype=bool)
    IdxL=np.zeros((LatLong3.shape[0], Idx.shape[1]), dtype=bool)    
    
    Distance3=DistanceBetweenLatLong(LatLong3)
    Idx3=Distance3>=Radious
    
    Idx1=np.concatenate((Idx,IdxR),axis=1)   
    Idx3=np.concatenate((IdxL,Idx3),axis=1)
    
    Idx=np.concatenate((Idx1,Idx3),axis=0)
    
    return Idx