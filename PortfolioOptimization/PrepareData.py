#Prepare the data for the portfolio optimization
#Set time range and peculiarities of the simulation
#Remember to convert the pu values to a common pu base or the actual values

import numpy as np
import datetime as dt


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


def PrepareData():
    #The wind data starts at 1/1/2007 and ends at 12/31/2013
    WindE_File=np.load("./WindEnergyNREL_100m_Haliade150_6MW.npz")# Wind energy File
    WindEnergy=WindE_File['WindEnergy']
    WindLatLong=WindE_File['LatLong']
    AnnualizedCostWind=WindE_File['AnnualizedCostWind']
    RatedPowerWind=float(WindE_File['RatedPower'])
    
    #The wave data starts at 1/1/2009 and ends at 12/31/2013
    WaveE_File=np.load('./WaveEnergy_Pelamis_2009_2013.npz')# Wave energy File
    WaveEnergy=WaveE_File['Energy_pu']
    WaveLatLong=WaveE_File['LatLong']
    RatedPowerWave=float(WaveE_File['RatedPower'])
    AnnualizedCostWave=WaveE_File['AnnualizedCostWave']
    
    
#-------------------------- WIND-----------------#
    SIdxWind=(dt.datetime(2009,1,1,0)-dt.datetime(2007,1,1,0)).days*24# Start index
    EIdxWind=(dt.datetime(2014,1,1,0)-dt.datetime(2007,1,1,0)).days*24# Last index
    
    Num3HourIntervals=int((EIdxWind-SIdxWind)/3)
    WindEnergy_Temp=np.zeros((WindEnergy.shape[0],Num3HourIntervals),dtype=float)
    
    for H3_Steps in range(Num3HourIntervals):
        
        StartIdx=SIdxWind + 3*(H3_Steps)
        
        WindEnergy_Temp[:,H3_Steps]=(WindEnergy[:,StartIdx])*RatedPowerWind
        
    WindEnergy=WindEnergy_Temp
    
    #Filter sites of wind energy with very low capacity factors
    AvgCFWind=np.average(WindEnergy,axis=1)/(RatedPowerWind)
    WindEnergy=WindEnergy[AvgCFWind>=0.20,:]
    WindLatLong=WindLatLong[AvgCFWind>=0.20,:]
    AnnualizedCostWind=AnnualizedCostWind[AvgCFWind>=0.20]
        
#-------------------------- WAVE-----------------#
    SIdxWave=(dt.datetime(2009,1,1,0)-dt.datetime(2009,1,1,0)).days*8
    EIdxWave=((dt.datetime(2014,1,1,0)-dt.datetime(2009,1,1,0)).days)*8
    WaveEnergy=WaveEnergy[:,SIdxWave:EIdxWave]*RatedPowerWave

    #Filter sites of wave energy with very low capacity factors
    AvgCFWave=np.average(WaveEnergy,axis=1)/(RatedPowerWave)
    
    WaveEnergy=WaveEnergy[AvgCFWave>=0.05,:]
    WaveLatLong=WaveLatLong[AvgCFWave>=0.05,:]
    AnnualizedCostWave=AnnualizedCostWave[AvgCFWave>=0.05]
    
#-------------------------- OCEAN CURRENT-----------------#   
#The wave data starts at 1/1/2009 and ends at 11/30/2014
    OceanE_File=np.load('./OceanCurrentEnergyRM4.npz')# Ocean current energy File
    OceanEnergy=OceanE_File['CurrentEnergy_pu']
    OceanLatLong=OceanE_File['LatLong']
    
    
    
    tempKite=np.load('./tempKite.npz')# Ocean current energy File
    
    Power         = tempKite["Power"] #[Kw]
    AnualizedCost = tempKite["AnualizedCost"] # $/year
    LatLong_Kit   = tempKite["LatLong_Kit"]
    
    
    
    Distance=DistanceBetweenLatLong(LatLong_Kit, OceanLatLong)
    
    #index of closest site location of hycom from a specific mabsab site location
    IdxLatLong=np.argmin(Distance, axis=1)
    OceanEnergy=OceanEnergy[IdxLatLong,:]
    
    AvgCF=np.mean(OceanEnergy,axis=1)
    RatedPower=Power/AvgCF
    RatedPowerOcean=np.reshape(RatedPower,(RatedPower.shape[0],1))*1000
    
    
    OceanEnergy=OceanEnergy*RatedPowerOcean
    OceanLatLong=LatLong_Kit
    AnnualizedCostOcean=AnualizedCost

    SIdxCurrent=(dt.datetime(2009,1,1,0)-dt.datetime(2009,1,1,0)).days*8
    EIdxCurrent=((dt.datetime(2014,1,1,0)-dt.datetime(2009,1,1,0)).days)*8
    OceanEnergy=OceanEnergy[:,SIdxCurrent:EIdxCurrent]
    

	
    Data={"WindEnergy":WindEnergy,
          "WindLatLong":WindLatLong,
          "AnnualizedCostWind":AnnualizedCostWind,
          "WaveEnergy":WaveEnergy,
          "WaveLatLong":WaveLatLong,
          "AnnualizedCostWave":AnnualizedCostWave,
          "OceanEnergy":OceanEnergy,
          "OceanLatLong":OceanLatLong,
          "AnnualizedCostOcean":AnnualizedCostOcean,
          "RatedPowerWind":RatedPowerWind,
          "RatedPowerWave":RatedPowerWave,
          "RatedPowerOcean":RatedPowerOcean
          }
    
    return Data
