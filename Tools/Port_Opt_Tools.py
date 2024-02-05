#env Gurobi
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
from tqdm import tqdm
import sys

#Determine the areas that are overlapping between two sets of lat/long coordinates.
#We need this to limit the number of turbines that are placed in the same area.
def GetOverlaps_Idx_Area(LatLong1, ResolutionKm1, ResolutionDegrees1,MaxNumTurbinesPerSite1,
                         LatLong2, ResolutionKm2, ResolutionDegrees2,MaxNumTurbinesPerSite2, SameTech=1, PrintName=""):
    IdxOverlap=[]
    AreaOverlap=[]
    AreaRef1Ref2=[]
    MaxTurbinesRef1Ref2=[]
    PercentageOverlap=[]
    
    if len(LatLong2)!=0 and len(LatLong1)!=0:
        print('Finding Overlap Site Locations for '+PrintName)
        for i in tqdm(range(len(LatLong1))):
            LatLong1_tmp=LatLong1[i,:]
            LatKmPerDegree1= 111.32 #How many km per lat degrees
            LongKmPerDegree1= (111.320 * np.cos(LatLong1_tmp[0] * np.pi / 180)) #How many km per long degrees
            
            if ResolutionKm1[i]!=-1: #Convert to degrees
                ResolutionDegrees1_Lat=ResolutionKm1[i]/LatKmPerDegree1
                ResolutionDegrees1_Long=ResolutionKm1[i]/LongKmPerDegree1
            else:
                ResolutionDegrees1_Lat=ResolutionDegrees1[i]
                ResolutionDegrees1_Long=ResolutionDegrees1[i]
                
            #Box Limiting the area of the site
            Lat1_1=LatLong1_tmp[0]-ResolutionDegrees1_Lat/2
            Lat1_2=LatLong1_tmp[0]+ResolutionDegrees1_Lat/2
            
            Long1_1=LatLong1_tmp[1]-ResolutionDegrees1_Long/2
            Long1_2=LatLong1_tmp[1]+ResolutionDegrees1_Long/2
            TotalAreaDeg1=(Lat1_2-Lat1_1)*LatKmPerDegree1*(Long1_2-Long1_1)*LongKmPerDegree1
            
            
            #Filter Sites that are  likely not overlapping (important for speed). This is an approximation, and if the difference
            # in resolution between sites is more than 10x this MAY not capture all the sites
            FilterIdx=(LatLong2[:,0]>(-10*ResolutionDegrees1_Lat+Lat1_1)) * (LatLong2[:,0]<(10*ResolutionDegrees1_Lat+Lat1_2)) \
                *(LatLong2[:,1]>(-10*ResolutionDegrees1_Long+Long1_1)) * (LatLong2[:,1]<(10*ResolutionDegrees1_Long+Long1_2))

            JsIn=np.where(FilterIdx==True)[0]#Index of sites that are in the box
            
            if SameTech==1:
                JsIn=JsIn[JsIn>i]  #Only check sites that have not been checked before
            
            if len(JsIn)!=0:
                for j in JsIn:
                    LatLong2_tmp=LatLong2[j,:]
                    LatKmPerDegree2= 111.32 #How many km per lat degrees
                    LongKmPerDegree2= (111.320 * np.cos(LatLong2_tmp[0] * np.pi / 180)) #How many km per long degrees
                    
                    if ResolutionKm2[j]!=-1: #Convert to degrees
                        ResolutionDegrees2_Lat=ResolutionKm2[j]/LatKmPerDegree2
                        ResolutionDegrees2_Long=ResolutionKm2[j]/LongKmPerDegree2
                    else:
                        ResolutionDegrees2_Lat=ResolutionDegrees2[j]
                        ResolutionDegrees2_Long=ResolutionDegrees2[j]
                        
                
                    #Box Limiting the area of the site
                    Lat2_1=LatLong2_tmp[0]-ResolutionDegrees2_Lat/2
                    Lat2_2=LatLong2_tmp[0]+ResolutionDegrees2_Lat/2
                    
                    Long2_1=LatLong2_tmp[1]-ResolutionDegrees2_Long/2
                    Long2_2=LatLong2_tmp[1]+ResolutionDegrees2_Long/2
                    
                    TotalAreaDeg2=(Lat2_2-Lat2_1)*LatKmPerDegree2*(Long2_2-Long2_1)*LongKmPerDegree2

                    if max(Lat1_1, Lat2_1) <= min(Lat1_2, Lat2_2) and max(Long1_1, Long2_1) <= min(Long1_2, Long2_2):
                        
                        #Overlap (Here LatKmPerDegree2 is very close to LatKmPerDegree1 se we can use either)
                        AreaOverlap_tmp=(min(Lat1_2, Lat2_2)-max(Lat1_1, Lat2_1))*LatKmPerDegree2*(min(Long1_2, Long2_2)-max(Long1_1, Long2_1))*LongKmPerDegree2
                        PercentageOverlap_tmp=AreaOverlap_tmp/TotalAreaDeg1

                        # if overlap is less than 5% of the total area of the site ignore
                        if PercentageOverlap_tmp>0.05:
                            IdxOverlap.append([i,j])
                            AreaOverlap.append(AreaOverlap_tmp)
                            AreaRef1Ref2.append([TotalAreaDeg1,TotalAreaDeg2])
                            MaxTurbinesRef1Ref2.append([MaxNumTurbinesPerSite1[i],MaxNumTurbinesPerSite2[j]])
                            PercentageOverlap.append(PercentageOverlap_tmp)
                                
    if len(LatLong2)==0 or len(LatLong1)==0:
        IdxOverlap=[]
        AreaOverlap=[]
        AreaRef1Ref2=[]
        MaxTurbinesRef1Ref2=[]
        PercentageOverlap=[]
        
    # For wind data some sites were duplicated on the donwload, and they show here as overlapping
    #This is not a problem as they are properly constrained in the optimization
    IdxOverlap=np.array(IdxOverlap)
    AreaOverlap=np.array(AreaOverlap)
    AreaRef1Ref2=np.array(AreaRef1Ref2)
    MaxTurbinesRef1Ref2=np.array(MaxTurbinesRef1Ref2)
    PercentageOverlap=np.array(PercentageOverlap)
    
    if len(IdxOverlap)==0:
        IdxOverlap=np.reshape(IdxOverlap,(0,2))
        IdxOverlap=IdxOverlap.astype(int)
        AreaOverlap=np.reshape(AreaOverlap,(0))
        AreaRef1Ref2=np.reshape(AreaRef1Ref2,(0,2))
        MaxTurbinesRef1Ref2=np.reshape(MaxTurbinesRef1Ref2,(0,2))
    
    return IdxOverlap, AreaOverlap, AreaRef1Ref2, MaxTurbinesRef1Ref2, PercentageOverlap