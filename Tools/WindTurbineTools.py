#Tools to manipulate wind turbine data
#geo_env
import csv
import numpy as np
import xarray as xr
import pandas as pd
from datetime import datetime, timedelta
from GeneralGeoTools import GetDepth, GetDistanceToShore, GetTimeList, MinDistanceSetPoints
import geopandas as gpd
from shapely.geometry import Polygon, Point
from tqdm import tqdm


#Read wind turbine data
def GetTurbineData(InputDataPath, WindTurbine):
    
    PowerCurve=[]
    PowerCurve.append(['0','0'])
    
    File_Wind=open(InputDataPath+"/Wind/"+WindTurbine+'.txt', "r")
    File_Wind_csv=csv.reader(File_Wind,delimiter=';')
    
    for EachLine in File_Wind_csv:
        
        if File_Wind_csv.line_num==1:
            HubHeight=float(EachLine[1])
            
        if File_Wind_csv.line_num==2:
            Total_Efficiency=1-float(EachLine[1])
            
        if File_Wind_csv.line_num==3:
            RotorDiameter=float(EachLine[1])

        if File_Wind_csv.line_num==4:
            RatedPower=float(EachLine[1])
            
        if File_Wind_csv.line_num>=6:
            PowerCurve.append(EachLine[:])

    PowerCurve=np.asarray(PowerCurve, dtype='float32')
    
    File_Wind.close()
    
    return HubHeight, Total_Efficiency, RotorDiameter, PowerCurve, RatedPower


def FilterOnDepthShoreDistance(Depth, DistanceShore, DepthMinMax, DistanceShoreMinMax):
    IdxIn=(DistanceShore>=DistanceShoreMinMax[0]) * (DistanceShore<=DistanceShoreMinMax[1]) * (Depth>=DepthMinMax[0]) * (Depth<=DepthMinMax[1])
    
    return IdxIn

#Convert wind speed to energy in pu
def WindToEnergy(InputDataPath, WindTurbine, WindSpeedHeightsAvailable, ResolutionKm=2, SavePath=None):
       
    HubHeight, Total_Efficiency, RotorDiameter, PowerCurve, RatedPower=GetTurbineData(InputDataPath, WindTurbine)

    HeightsAvaialble=np.asarray(list(WindSpeedHeightsAvailable.keys()),dtype="int")

    #Get files we are going to use to perfrom conversion of wind speed to the hub height
    #ul/uh=(zl/zh)^alpha

    IdxUpperWindFile=np.argmin(np.abs(HeightsAvaialble-HubHeight)) #Get the closest height available
    if HeightsAvaialble[IdxUpperWindFile]-HubHeight<0:
        IdxUpperWindFile=min(IdxUpperWindFile+1,len(HeightsAvaialble)-1)

    IdxUpperWindFile=max(IdxUpperWindFile,1)
    IdxLowerWindFile=IdxUpperWindFile-1

    LWindFileName=WindSpeedHeightsAvailable[list(WindSpeedHeightsAvailable.keys())[IdxLowerWindFile]]
    HWindFileName=WindSpeedHeightsAvailable[list(WindSpeedHeightsAvailable.keys())[IdxUpperWindFile]]



    LWindData=np.load(InputDataPath+"/Wind/"+LWindFileName)["windspeed"]
    HWindData=np.load(InputDataPath+"/Wind/"+HWindFileName)["windspeed"]
    LatLong=np.load(InputDataPath+"/Wind/"+LWindFileName)["LatLong"]
       
    Depth=GetDepth(InputDataPath, LatLong)
    DistanceShore=GetDistanceToShore(InputDataPath, LatLong)
    
    #Exclude points that are too shallow or too far from shore
    IdxIn=FilterOnDepthShoreDistance(Depth, DistanceShore, DepthMinMax=(6,2500), DistanceShoreMinMax=(10,2000))
    LWindData=LWindData[:,IdxIn]
    HWindData=HWindData[:,IdxIn]
    LatLong=LatLong[IdxIn,:]
    Depth=Depth[IdxIn]
    DistanceShore=DistanceShore[IdxIn]
    
    
    # HubHeight
    LWindData[LWindData<0.01]=0.01 #To avoid log(0)
    HWindData[HWindData<0.01]=0.01 #To avoid log(0)
    alpha=np.log(HWindData/LWindData)/np.log(HeightsAvaialble[IdxUpperWindFile]/HeightsAvaialble[IdxLowerWindFile])
    alpha=np.median(alpha)
    
    if HeightsAvaialble[IdxUpperWindFile]<=HubHeight:
        h_ref=HeightsAvaialble[IdxUpperWindFile]
        WS_Ref=HWindData
    else:
        h_ref=HeightsAvaialble[IdxLowerWindFile]
        WS_Ref=LWindData
    
    WS_Hub=WS_Ref*(HubHeight/h_ref)**alpha #Wind speed at the hub height

    WindEnergy_pu=np.interp(WS_Hub,PowerCurve[:,0],PowerCurve[:,1])

    ReadMe='\
    EnergyPu: pu wind energy \n\
    1) The data is in hourly discretization starting at 1/1/2007 and going up to\
     12/31/2013 23:00'
    
    RatedPower=RatedPower
    #Attention Time is Hard Coded Here
    TimeList=GetTimeList(datetime(2007, 1, 1), datetime(2013, 12, 31, 23), TimeDeltaHours=1)
    
    if SavePath!=None:
        np.savez(SavePath, ReadMe=ReadMe, Energy_pu=WindEnergy_pu.astype(np.float16), RatedPower=RatedPower, LatLong=LatLong.astype(np.float32), Depth=Depth,\
            DistanceShore=DistanceShore.astype(np.float32),RawResource=WS_Hub.astype(np.float16),TimeList=TimeList,ResolutionKm=ResolutionKm )
        

    return WindEnergy_pu, RatedPower, LatLong, WS_Hub, Depth, DistanceShore, TimeList, ResolutionKm

def GetCostAndGenerationWindTurbine(InputDataPath, WindCostPath, WindTurbine, WindSpeedHeightsAvailable, SavePath=None):

    WindEnergy_pu, RatedPower, LatLong, WS_Hub, Depth, DistanceShore, TimeList, ResolutionKm=WindToEnergy(InputDataPath, WindTurbine, WindSpeedHeightsAvailable,ResolutionKm=2, SavePath=None)
        
    #Read EXCEL data and get NREL information
    CAPEX=pd.read_excel(WindCostPath,sheet_name="CAPEX")
    OPEX=pd.read_excel(WindCostPath,sheet_name="OPEX")

    ATB_SiteToLandfall=CAPEX["SiteToLandfall"]
    ATB_AvgDepth=CAPEX["AvgDepth"]
    ATB_TRG=CAPEX["TRG"]

    CAPEX=CAPEX[WindTurbine]
    OPEX=OPEX[WindTurbine]

    CAPEX=CAPEX*RatedPower*(10**-3) #[M$]
    OPEX =OPEX*RatedPower*(10**-3) #[$/Year]

    FCR=6.8/100 # Factor of Capital Return(WECC ATB-2023)
    IdxIn=Depth<=np.max(ATB_AvgDepth)*1.1 #Maxmimum depth which we have cost estimates
            
    #Filter for the region where we have cost estimated
    WindEnergy_pu=WindEnergy_pu[:,IdxIn]
    LatLong=LatLong[IdxIn,:]
    WS_Hub=WS_Hub[:,IdxIn]
    Depth=Depth[IdxIn] 
    DistanceShore=DistanceShore[IdxIn]


    Diff_Depth=np.reshape(Depth,(len(Depth),1))-np.reshape(ATB_AvgDepth,(1,len(ATB_AvgDepth)))
    Diff_Distance=np.reshape(DistanceShore,(len(DistanceShore),1))-np.reshape(ATB_SiteToLandfall,(1,len(ATB_SiteToLandfall)))

    EuclidianDistance=(Diff_Depth**2+Diff_Distance**2)
    Idx_NRELTurbine=np.argmin(EuclidianDistance,axis=1)

    TRG_site=ATB_TRG[Idx_NRELTurbine]
    CAPEX_site=CAPEX[Idx_NRELTurbine]
    OPEX_site=OPEX[Idx_NRELTurbine]

    AnnualizedCost=CAPEX_site*FCR + OPEX_site

    ReadMe='\
    EnergyPu: pu wind energy \n\
    1) The data is in hourly discretization starting at 1/1/2007 and going up to\
        12/31/2013 23:00\n \
    2) Cost Values in M$, energy values in MW'
    
    #Attention Time is Hard Coded Here
    TimeList=GetTimeList(datetime(2007, 1, 1), datetime(2013, 12, 31, 23), TimeDeltaHours=1)
    
    AnnualizedCost=AnnualizedCost/(10**6)
    CAPEX_site=CAPEX_site/(10**6)
    OPEX_site=OPEX_site/(10**6)
    RatedPower=RatedPower/(10**6)
    NumberOfCellsPerSite=np.ones(LatLong.shape[0]) #one site per cell

    if SavePath!=None:
        np.savez(SavePath, ReadMe=ReadMe, Energy_pu=WindEnergy_pu.astype(np.float16), RatedPower=RatedPower, LatLong=LatLong.astype(np.float32), Depth=Depth,\
            DistanceShore=DistanceShore.astype(np.float16), TRG_site=TRG_site, CAPEX_site=CAPEX_site.astype(np.float16), OPEX_site=OPEX_site.astype(np.float16),\
                AnnualizedCost=AnnualizedCost.astype(np.float16), RawResource=WS_Hub.astype(np.float16), TimeList=TimeList, NumberOfCellsPerSite=NumberOfCellsPerSite,
                ResolutionKm=ResolutionKm, ResolutionDegrees=-1)
        
    return WindEnergy_pu, RatedPower, LatLong, WS_Hub, Depth, DistanceShore, TRG_site, CAPEX_site, OPEX_site, AnnualizedCost, TimeList, NumberOfCellsPerSite


# Filter for WTK-NREL sites   
def FilterForWTKDataset(ReferenceDataPath, WTKPath, SavePath=None):
    Data=np.load(ReferenceDataPath,allow_pickle=True)

    Energy_Pu=Data["Energy_pu"]
    RawResource=Data["RawResource"]
    RatedPower=Data["RatedPower"]
    NumberOfCellsPerSite=Data["NumberOfCellsPerSite"]

    TimeList=Data["TimeList"]
    LatLong=Data["LatLong"]
    Depth=Data["Depth"]
    DistanceShore=Data["DistanceShore"]
    CAPEX_site=Data["CAPEX_site"]
    OPEX_site=Data["OPEX_site"]
    AnnualizedCost=Data["AnnualizedCost"]
    ResolutionDegrees=Data["ResolutionDegrees"]
    ResolutionKm=Data["ResolutionKm"]
    


    #def FilterUsingWTKSites(WTKPath):  
    NRELLimitSites=pd.read_csv(WTKPath)
    NREL_Lat=NRELLimitSites["latitude"].values
    NREL_Long=NRELLimitSites["longitude"].values
    fraction_of_usable_area=NRELLimitSites["fraction_of_usable_area"].values
    power_curve= NRELLimitSites["power_curve"].values


    MaxLat=np.max(LatLong[:,0])
    MinLat=np.min(LatLong[:,0])
    MaxLong=np.max(LatLong[:,1])
    MinLong=np.min(LatLong[:,1])

    IdxNRELFilter=(NREL_Lat<=MaxLat) & (NREL_Lat>=MinLat) & (NREL_Long<=MaxLong) & (NREL_Long>=MinLong)\
        & (fraction_of_usable_area==1) &(power_curve=="offshore")

    NREL_Lat=NREL_Lat[IdxNRELFilter]
    NREL_Long=NREL_Long[IdxNRELFilter]
    NREL_LatLong=np.stack((NREL_Lat, NREL_Long),axis=1)

    IdxMin, DMin=MinDistanceSetPoints(LatLong, NREL_LatLong)
    IdxMin=IdxMin[DMin<=3]# Keep all the points that distance less than 3km from a NREL site


    Energy_Pu=Energy_Pu[:,IdxMin]
    RawResource=RawResource[:,IdxMin]

    LatLong=LatLong[IdxMin,:]
    Depth=Depth[IdxMin]
    DistanceShore=DistanceShore[IdxMin]
    CAPEX_site=CAPEX_site[IdxMin]
    OPEX_site=OPEX_site[IdxMin]
    AnnualizedCost=AnnualizedCost[IdxMin]

    DataDir={"Energy_pu":Energy_Pu,
            "RawResource":RawResource,
            "TimeList":TimeList,
            "LatLong":LatLong,
            "Depth":Depth,
            "DistanceShore":DistanceShore,
            "CAPEX_site":CAPEX_site,
            "OPEX_site":OPEX_site,
            "AnnualizedCost":AnnualizedCost}
    
    if SavePath!=None:
        np.savez(SavePath,  Energy_pu=Energy_Pu.astype(np.float16), RatedPower=RatedPower, LatLong=LatLong.astype(np.float32), Depth=Depth,\
            DistanceShore=DistanceShore.astype(np.float16), CAPEX_site=CAPEX_site.astype(np.float16), OPEX_site=OPEX_site.astype(np.float16),\
                AnnualizedCost=AnnualizedCost.astype(np.float16), RawResource=RawResource.astype(np.float16), TimeList=TimeList, NumberOfCellsPerSite=NumberOfCellsPerSite,
                ResolutionDegrees=ResolutionDegrees, ResolutionKm=ResolutionKm)
    
    return DataDir

def FilterForBOEM(ReferenceDataPath, PathBOEMData,BOEM_ShpDir, SavePath=None):
    
    Data=np.load(ReferenceDataPath,allow_pickle=True)

    Energy_Pu=Data["Energy_pu"]
    RawResource=Data["RawResource"]
    RatedPower=Data["RatedPower"]
    NumberOfCellsPerSite=Data["NumberOfCellsPerSite"]

    TimeList=Data["TimeList"]
    LatLong=Data["LatLong"]
    Depth=Data["Depth"]
    DistanceShore=Data["DistanceShore"]
    CAPEX_site=Data["CAPEX_site"]
    OPEX_site=Data["OPEX_site"]
    AnnualizedCost=Data["AnnualizedCost"]
    ResolutionDegrees=Data["ResolutionDegrees"]
    ResolutionKm=Data["ResolutionKm"]

    SHP_WindLease = gpd.read_file(PathBOEMData+BOEM_ShpDir["Wind Lease Areas"])
    SHP_WindPlanning = gpd.read_file(PathBOEMData+BOEM_ShpDir["Wind Planning Areas"])

    SHP_WindLease_U=SHP_WindLease["geometry"][0]
    SHP_WindPlanning_U=SHP_WindPlanning["geometry"][0]

    for i in range(1,len(SHP_WindLease["geometry"])):
        SHP_WindLease_U=SHP_WindLease_U.union(SHP_WindLease["geometry"][i])

    for i in range(1,len(SHP_WindPlanning["geometry"])):
        SHP_WindPlanning_U=SHP_WindPlanning_U.union(SHP_WindPlanning["geometry"][i])



    SHP_LongLatPoint=[Point(x,y) for x,y in zip(LatLong[:,1],LatLong[:,0])]
    IdxIn=[]
    Code_BOEM=[]
    Distance=[]
    print("Filtering sites in lease/planning areas")
    for i in tqdm(range(len(SHP_LongLatPoint))):
        Distance.append(SHP_LongLatPoint[i].distance(SHP_WindLease_U))
        
        if SHP_LongLatPoint[i].distance(SHP_WindLease_U)<0.001:
            IdxIn.append(i)
            Code_BOEM.append("Lease")
            
            
        elif SHP_LongLatPoint[i].distance(SHP_WindPlanning_U)<0.001:
            IdxIn.append(i)
            Code_BOEM.append("Planning")
        
        else:
            Code_BOEM.append("Outside")



    Code_BOEM=np.array(Code_BOEM)

    Code_BOEM=Code_BOEM[IdxIn]
    Energy_Pu=Energy_Pu[:,IdxIn]
    RawResource=RawResource[:,IdxIn]
    NumberOfCellsPerSite=NumberOfCellsPerSite[IdxIn]
    LatLong=LatLong[IdxIn,:]
    Depth=Depth[IdxIn]
    DistanceShore=DistanceShore[IdxIn]
    CAPEX_site=CAPEX_site[IdxIn]
    OPEX_site=OPEX_site[IdxIn]
    AnnualizedCost=AnnualizedCost[IdxIn]

    DataDir={"Energy_pu":Energy_Pu,
            "RawResource":RawResource,
            "TimeList":TimeList,
            "LatLong":LatLong,
            "Depth":Depth,
            "DistanceShore":DistanceShore,
            "CAPEX_site":CAPEX_site,
            "OPEX_site":OPEX_site,
            "AnnualizedCost":AnnualizedCost,
            "Code_BOEM":Code_BOEM}

    if SavePath!=None:
        np.savez(SavePath,  Energy_pu=Energy_Pu.astype(np.float16), RatedPower=RatedPower, LatLong=LatLong.astype(np.float32), Depth=Depth,\
            DistanceShore=DistanceShore.astype(np.float16), CAPEX_site=CAPEX_site.astype(np.float16), OPEX_site=OPEX_site.astype(np.float16),\
                AnnualizedCost=AnnualizedCost.astype(np.float16), RawResource=RawResource.astype(np.float16), TimeList=TimeList, NumberOfCellsPerSite=NumberOfCellsPerSite, Code_BOEM=Code_BOEM,
                ResolutionDegrees=ResolutionDegrees, ResolutionKm=ResolutionKm)
        
    return DataDir




