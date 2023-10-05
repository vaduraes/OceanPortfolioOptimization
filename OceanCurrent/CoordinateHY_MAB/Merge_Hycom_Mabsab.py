#Merge the data from HYCOM and MABSAB, keep the better resolution of MABSAB but allow the representation of
#hourly variability in the ocean current resource

import numpy as np
import csv

MABSAB_Dictionary=np.load("../MABSAB_DATA/CurrentSpeedMABSAB2009_2014_50m.npz", allow_pickle=True)

LatLong_MABSAB = MABSAB_Dictionary["LatLong"]
CurrentSpeed_MABSAB = MABSAB_Dictionary["CurrentSpeed"]
RatedPower_MABSAB = MABSAB_Dictionary["RatedPower"]
DepthSites_MABSAB = MABSAB_Dictionary["DepthSites"]
ShoreDistance_MABSAB = MABSAB_Dictionary["ShoreDistance"]
OceanDateTime_MABSAB = MABSAB_Dictionary["OceanDateTime"]


HYCOM_Dictionary=np.load("../HYCOM_RAW_DATA/CurrentSpeedHYCOM2009_2014.npz", allow_pickle=True)
CurrentSpeedByDay_HYCOM=HYCOM_Dictionary['SpeedByDay_AllLatLong']
CurrentSpeed_HYCOM=HYCOM_Dictionary['Speed']
LatLong_HYCOM=HYCOM_Dictionary['LatLong']
DateTime_HYCOM=HYCOM_Dictionary["DateTime"]


#Get coastline data
CoastLine=[]

File_CoastLine=open('../GEO_data/Coastline_NC.csv', "r")
File_CoastLine_csv=csv.reader(File_CoastLine,delimiter=',')

for EachLine in File_CoastLine_csv:

    if File_CoastLine_csv.line_num > 1:
        CoastLine.append([float(EachLine[1]), float(EachLine[0])] ) #LatLong

CoastLine=np.array(CoastLine)

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

#For a deplyment of 200MVA
def EfficiencyTransmisison(LatLong):
    E_Transmission=[]
    for LatLong_i in LatLong:
        D=DistanceToShore(CoastLine, LatLong_i)
        if D<66:
            E_Transmission.append(-0.0362*D + 98.804)#HVAC
        else:
            E_Transmission.append(-0.0206*D + 96.432)#HVDC
            
    E_Transmission=np.array(E_Transmission)/100
    E_Transmission=np.reshape(E_Transmission,(len(E_Transmission),1))
    
    return E_Transmission

def DistanceBetweenLatLong(LatLong1, LatLong2):
    LatLong1=LatLong1*2*np.pi/360
    LatLong2=LatLong2*2*np.pi/360
    
    dLat=np.reshape(LatLong1[:,0],(len(LatLong1[:,0]),1))-np.reshape(LatLong2[:,0],(1,len(LatLong2[:,0])))
    dLong=np.reshape(LatLong1[:,1],(len(LatLong1[:,1]),1))-np.reshape(LatLong2[:,1],(1,len(LatLong2[:,1])))
    
    P1=np.repeat(np.reshape(np.cos(LatLong1[:,0]),(LatLong1.shape[0],1)),LatLong2.shape[0],axis=1)
    P2=np.repeat(np.reshape(np.cos(LatLong2[:,0]),(1,LatLong2.shape[0])),LatLong1.shape[0],axis=0)
    
    a=np.power(np.sin(dLat/2),2) + P1*P2*np.power(np.sin(dLong/2),2)
    c=2*np.arcsin(np.minimum(1,np.sqrt(a)))
    Distance=6367*c #[km]
    
    return Distance

#Making hycom and mabsab consistent in their time scale. They both start at 1/1/2009 but MABSAB end later than HYCOM
OceanDateTime_MABSAB=OceanDateTime_MABSAB[0:CurrentSpeedByDay_HYCOM.shape[1]]
CurrentSpeed_MABSAB=CurrentSpeed_MABSAB[:,0:CurrentSpeedByDay_HYCOM.shape[1]]
DailyAverage_Hycom=np.mean(CurrentSpeedByDay_HYCOM,axis=2)


np.savez("HycomMabsab", OceanDateTime_MABSAB=OceanDateTime_MABSAB, CurrentSpeed_MABSAB=CurrentSpeed_MABSAB,
         LatLong_MABSAB=LatLong_MABSAB,DepthSites_MABSAB=DepthSites_MABSAB, ShoreDistance=ShoreDistance_MABSAB,
         DateTime_HYCOM=DateTime_HYCOM, DailyAverage_Hycom=DailyAverage_Hycom,CurrentSpeedByDay_HYCOM=CurrentSpeedByDay_HYCOM,
         LatLong_HYCOM=LatLong_HYCOM)