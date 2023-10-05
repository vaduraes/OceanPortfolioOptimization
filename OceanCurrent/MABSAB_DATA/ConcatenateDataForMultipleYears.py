#Concatenate in one npz file the data from multiple years
#The data needs to have the same format and same site locations
import numpy as np

#YearsToConcatenate=[2009,2010,2011,2012,2013,2014]
#FileSaveName='CurrentEnergyRM4_'

def ConcatenateData(YearsToConcatenate,FileSaveName, RotorDepth):

    YearIdx=0
    for YearIdx in range(len(YearsToConcatenate)):
        
        Data=np.load(FileSaveName+str(YearsToConcatenate[YearIdx])+'.npz',allow_pickle=True)
        
        if YearIdx==0:
            CurrentSpeed=Data['CurrentSpeed']
            RatedPower=Data['RatedPower']
            LatLong=Data['LatLong']
            OceanDateTime=Data['OceanDateTime']
            DepthSites=Data['DepthSites']
            ShoreDistance=Data['ShoreDistance']
            
        
        if YearIdx>=1:
            OceanDateTime=np.concatenate((OceanDateTime,Data['OceanDateTime']))
            CurrentSpeed=np.concatenate((CurrentSpeed, Data['CurrentSpeed']), axis=1)
    
        ReadMe='\
        CurrentSpeed: Ocean current speed\
        LatLong: Latitude,Logitude data\n\
        1) The data is in Daily discretization starting at 1/1/'+str(YearsToConcatenate[0])+' and going up to \
        12/31/'+str(YearsToConcatenate[len(YearsToConcatenate)-1])       
        
    np.savez(FileSaveName+str(YearsToConcatenate[0])+'_'+str(YearsToConcatenate[len(YearsToConcatenate)-1])+"_"+str(int(RotorDepth))+"m"+'.npz',ReadMe=ReadMe,CurrentSpeed=CurrentSpeed,
             RatedPower=RatedPower,DepthSites=DepthSites,ShoreDistance=ShoreDistance,LatLong=LatLong,OceanDateTime=OceanDateTime)
