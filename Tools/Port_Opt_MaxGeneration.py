#This code compute the optimal portfolio for wave wind and ocean current resources
#Considering transmission system costs, CAPEX and OPEX of each technology and its generation availability in a given region

#The objective function is the maximization of the total generation of the portfolio, costraint to limits in the portfolio LCOE, maximum
#Capacity of the transmission system, maximum number of turbines per site location, and maxmimum radious of the energy collection system.

#The model also takes into considering curtailment, and the possibility of chosing from a limited number of turbine desings.

#env Gurobi
import numpy as np
from pyomo.environ import *
import pandas as pd
from datetime import datetime, timedelta
from tqdm import tqdm
import sys

from  GetIdxInOutRadious import GetIdxOutRadious, GetIdxInRadious_Simple
from Port_Opt_Tools import GetOverlaps_Idx_Area


def PreparePotOptInputs(PathWindDesigns, PathWaveDesigns, PathKiteDesigns, PathTransmissionDesign, LCOE_RANGE=range(200,30,-2)\
    ,Max_CollectionRadious=30, MaxDesignsWind=1, MaxDesignsWave=1, MaxDesignsKite=1, MinNumWindTurb=0, MinNumWaveTurb=0, MinNumKiteTrub=0,
    WindTurbinesPerSite=4, KiteTurbinesPerSite=4, WaveTurbinesPerSite=4):

    #WindTurbinesPerSite: Number of turbines per site location based on the initial wind resolution from NREL
    #KiteTurbinesPerSite: Number of turbines per site location based on the initial kite resolution from where the data was obtained (HYCOM, MABSAB)
    #WaveTurbinesPerSite: Number of turbines per site location based on the initial wave resolution from where the data was obtained (WWIII)  
    
    #Function to prepare the inputs for the optimization
    #All portfolio data needs to be at the same time resolution and range, unless the portfolio path is empty eg. PathWaveDesigns=[]
    # LCOE_RANGE=range(200,30,-2) #Max LCOE limits investigated
    # Max_CollectionRadious=30 #Radious for the energy collection system

    WindEnergy, WindLatLong, AnnualizedCostWind, MaxNumWindPerSite, WindDesign, TimeWindData,\
    RatedPowerWindTurbine, WindResolutionDegrees, WindResolutionKm = list(), list(), list(), list(), list(), list(), list(), list(), list()
    
    KiteEnergy, KiteLatLong, AnnualizedCostKite, MaxNumKitePerSite, KiteDesign, TimeKiteData,\
    RatedPowerKiteTurbine, KiteResolutionDegrees, KiteResolutionKm = list(), list(), list(), list(), list(), list(), list(), list(), list()
    
    WaveEnergy, WaveLatLong, AnnualizedCostWave, MaxNumWavePerSite, WaveDesign, TimeWaveData,\
    RatedPowerWaveTurbine, WaveResolutionDegrees, WaveResolutionKm = list(), list(), list(), list(), list(), list(), list(), list(), list()
    
    for i in range(len(PathWindDesigns)):
        Data=np.load(PathWindDesigns[i],allow_pickle=True)
        if i==0:
            WindEnergy=Data['Energy_pu']
            WindLatLong=Data['LatLong']
            AnnualizedCostWind=Data['AnnualizedCost']
            MaxNumWindPerSite=Data["NumberOfCellsPerSite"]*WindTurbinesPerSite
            WindDesign=np.array([i]*len(Data["NumberOfCellsPerSite"]))
            TimeWindData=Data["TimeList"]
            RatedPowerWindTurbine=np.array([float(Data["RatedPower"])]*len(Data["NumberOfCellsPerSite"]))
            WindResolutionDegrees=np.array([float(Data["ResolutionDegrees"])]*len(Data["NumberOfCellsPerSite"]))
            WindResolutionKm=np.array([float(Data["ResolutionKm"])]*len(Data["NumberOfCellsPerSite"]))
            
        else:
            WindEnergy=np.concatenate((WindEnergy,Data['Energy_pu']),axis=1)
            WindLatLong=np.concatenate((WindLatLong,Data['LatLong']))
            AnnualizedCostWind=np.concatenate((AnnualizedCostWind,Data['AnnualizedCost']))
            MaxNumWindPerSite=np.concatenate((MaxNumWindPerSite,WindTurbinesPerSite*Data["NumberOfCellsPerSite"]))
            WindDesign=np.concatenate((WindDesign,[i]*len(Data["NumberOfCellsPerSite"])))
            RatedPowerWindTurbine=np.concatenate((RatedPowerWindTurbine,np.array([float(Data["RatedPower"])]*len(Data["NumberOfCellsPerSite"]))))
            WindResolutionDegrees=np.concatenate((WindResolutionDegrees, np.array([float(Data["ResolutionDegrees"])]*len(Data["NumberOfCellsPerSite"]))))
            WindResolutionKm=np.concatenate((WindResolutionKm, np.array([float(Data["ResolutionKm"])]*len(Data["NumberOfCellsPerSite"]))))
            
        Data.close()
        
    #Kite Data
    for i in range(len(PathKiteDesigns)):
        Data=np.load(PathKiteDesigns[i],allow_pickle=True)
        if i==0:
            KiteEnergy=Data['Energy_pu']
            KiteLatLong=Data['LatLong']
            AnnualizedCostKite=Data['AnnualizedCost']
            MaxNumKitePerSite=Data["NumberOfCellsPerSite"]*KiteTurbinesPerSite
            KiteDesign=np.array([i]*len(Data["NumberOfCellsPerSite"]))
            TimeKiteData=Data["TimeList"]
            RatedPowerKiteTurbine=np.array([float(Data["RatedPower"])]*len(Data["NumberOfCellsPerSite"]))
            KiteResolutionDegrees=np.array([float(Data["ResolutionDegrees"])]*len(Data["NumberOfCellsPerSite"]))
            KiteResolutionKm=np.array([float(Data["ResolutionKm"])]*len(Data["NumberOfCellsPerSite"]))
            
            
        else:
            KiteEnergy=np.concatenate((KiteEnergy,Data['Energy_pu']),axis=1)
            KiteLatLong=np.concatenate((KiteLatLong,Data['LatLong']))
            AnnualizedCostKite=np.concatenate((AnnualizedCostKite,Data['AnnualizedCost']))
            MaxNumKitePerSite=np.concatenate((MaxNumKitePerSite,KiteTurbinesPerSite*Data["NumberOfCellsPerSite"]))
            KiteDesign=np.concatenate((KiteDesign,[i]*len(Data["NumberOfCellsPerSite"])))
            RatedPowerKiteTurbine=np.concatenate((RatedPowerKiteTurbine,np.array([float(Data["RatedPower"])]*len(Data["NumberOfCellsPerSite"]))))
            KiteResolutionDegrees=np.concatenate((KiteResolutionDegrees, np.array([float(Data["ResolutionDegrees"])]*len(Data["NumberOfCellsPerSite"]))))
            KiteResolutionKm=np.concatenate((KiteResolutionKm, np.array([float(Data["ResolutionKm"])]*len(Data["NumberOfCellsPerSite"]))))
            
        Data.close()
        
    #Wave Data
    for i in range(len(PathWaveDesigns)):
        Data=np.load(PathWaveDesigns[i],allow_pickle=True)
        if i==0:
            WaveEnergy=Data['Energy_pu']
            WaveLatLong=Data['LatLong']
            AnnualizedCostWave=Data['AnnualizedCost']
            MaxNumWavePerSite=Data["NumberOfCellsPerSite"]*WaveTurbinesPerSite
            WaveDesign=np.array([i]*len(Data["NumberOfCellsPerSite"]))
            TimeWaveData=Data["TimeList"]
            RatedPowerWaveTurbine=np.array([float(Data["RatedPower"])]*len(Data["NumberOfCellsPerSite"]))  
            WaveResolutionDegrees=np.array([float(Data["ResolutionDegrees"])]*len(Data["NumberOfCellsPerSite"]))
            WaveResolutionKm=np.array([float(Data["ResolutionKm"])]*len(Data["NumberOfCellsPerSite"]))
            
        else:
            WaveEnergy=np.concatenate((WaveEnergy,Data['Energy_pu']),axis=1)
            WaveLatLong=np.concatenate((WaveLatLong,Data['LatLong']))
            AnnualizedCostWave=np.concatenate((AnnualizedCostWave,Data['AnnualizedCost']))
            MaxNumWavePerSite=np.concatenate((MaxNumWavePerSite,WaveTurbinesPerSite*Data["NumberOfCellsPerSite"]))
            WaveDesign=np.concatenate((WaveDesign,[i]*len(Data["NumberOfCellsPerSite"])))
            RatedPowerWaveTurbine=np.concatenate((RatedPowerWaveTurbine,np.array([float(Data["RatedPower"])]*len(Data["NumberOfCellsPerSite"]))))
            WaveResolutionDegrees=np.concatenate((WaveResolutionDegrees, np.array([float(Data["ResolutionDegrees"])]*len(Data["NumberOfCellsPerSite"]))))
            WaveResolutionKm=np.concatenate((WaveResolutionKm, np.array([float(Data["ResolutionKm"])]*len(Data["NumberOfCellsPerSite"]))))

    #Verify if all the data is at the same time resolution and range
    if len(PathWindDesigns)!=0 and len(PathKiteDesigns)!=0:
        
        if np.all(TimeWindData==TimeKiteData)==False:
            return print("Time resolution of the wind, and wave data is not the same")
        
    if len(PathWindDesigns)!=0 and len(PathWaveDesigns)!=0:
        if  np.all(TimeWindData==TimeWaveData)==False:
            return print("Time resolution of the wind, and wave data is not the same")

    if len(PathKiteDesigns)!=0 and len(PathWaveDesigns)!=0:
        if  np.all(TimeKiteData==TimeWaveData)==False:
            return print("Time resolution of the kite, and wave data is not the same")

    TimeList=TimeWindData

    #Transmission
    Data=np.load(PathTransmissionDesign,allow_pickle=True)["TransmissionLineParameters"].item()

    AnnualizedCostTransmission=Data['S_BestACost']
    TransLatLong=Data['TL_LatLong']
    EfficiencyTransmission=Data['S_Efficiency']
    RatedPowerMWTransmissionMW=Data['RatedPowerMW']



    PortImputDir={  #Wind data
                    "WindEnergy":WindEnergy,
                    "WindLatLong":WindLatLong,
                    "AnnualizedCostWind":AnnualizedCostWind,
                    "MaxNumWindPerSite":MaxNumWindPerSite,
                    "WindDesign":WindDesign,
                    "RatedPowerWindTurbine":RatedPowerWindTurbine,
                    "NumWindSites": len(WindLatLong),
                    "WindResolutionDegrees":WindResolutionDegrees,
                    "WindResolutionKm":WindResolutionKm,
                    
                    
                    #Kite data
                    "KiteEnergy":KiteEnergy,
                    "KiteLatLong":KiteLatLong,
                    "AnnualizedCostKite":AnnualizedCostKite,
                    "MaxNumKitePerSite":MaxNumKitePerSite,
                    "KiteDesign":KiteDesign,
                    "RatedPowerKiteTurbine":RatedPowerKiteTurbine,
                    "NumKiteSites": len(KiteLatLong),
                    "KiteResolutionDegrees":KiteResolutionDegrees,
                    "KiteResolutionKm":KiteResolutionKm,
                    
                                        
                    #Wavedata
                    "WaveEnergy":WaveEnergy,
                    "WaveLatLong":WaveLatLong,
                    "AnnualizedCostWave":AnnualizedCostWave,
                    "MaxNumWavePerSite":MaxNumWavePerSite,
                    "WaveDesign":WaveDesign,
                    "RatedPowerWaveTurbine":RatedPowerWaveTurbine,
                    "NumWaveSites": len(WaveLatLong),
                    "WaveResolutionDegrees":WaveResolutionDegrees,
                    "WaveResolutionKm":WaveResolutionKm,
                    
                    
                    "TimeList":TimeList,
                    "NumTimeSteps":len(TimeList),
                    
                    #Transmission
                    "RatedPowerMWTransmissionMW":RatedPowerMWTransmissionMW,
                    "AnnualizedCostTransmission":AnnualizedCostTransmission,
                    "TransLatLong":TransLatLong,
                    "EfficiencyTransmission":EfficiencyTransmission,
                    "NumTransSites": len(TransLatLong),
                    
                    #Optimization Params
                    "LCOE_RANGE":LCOE_RANGE,
                    "Max_CollectionRadious":Max_CollectionRadious,
                    "MaxDesignsWind":MaxDesignsWind,
                    "MaxDesignsWave":MaxDesignsWave,
                    "MaxDesignsKite":MaxDesignsKite,
                    "MinNumWindTurb":MinNumWindTurb,
                    "MinNumWaveTurb":MinNumWaveTurb,
                    "MinNumKiteTrub":MinNumKiteTrub,
                
                }
    return PortImputDir

def SolvePortOpt_MaxGen_Model(PathWindDesigns, PathWaveDesigns, PathKiteDesigns, PathTransmissionDesign, LCOE_RANGE\
    ,Max_CollectionRadious,MaxDesignsWind, MaxDesingsWave, MaxDesingsKite,MinNumWindTurb,MinNumWaveTurb,MinNumKiteTrub):


    #Create and solve the optimization problem
    InputDir=PreparePotOptInputs(PathWindDesigns, PathWaveDesigns,PathKiteDesigns, PathTransmissionDesign, LCOE_RANGE\
        ,Max_CollectionRadious,MaxDesignsWind, MaxDesingsWave, MaxDesingsKite,MinNumWindTurb,MinNumWaveTurb,MinNumKiteTrub)

    NumWindDesigns=len(PathWindDesigns)
    NumWaveDesigns=len(PathWaveDesigns)
    NumKiteDesigns=len(PathKiteDesigns)

    Model = ConcreteModel()
    BigM=1000 #Big M for the maximum total number of turbines allowed to be installed (This is used to limits deployments in the radious of the energy collection system)


    # Create Variables
    Model.Y_Wind = Var(range(InputDir["NumWindSites"]), domain=NonNegativeIntegers)# Integer variable to track the number of wind turbines used per site location
    Model.Y_Wave = Var(range(InputDir["NumWaveSites"]) , domain=NonNegativeIntegers)# Integer variable to track the number of wave turbines used per site location
    Model.Y_Kite = Var(range(InputDir["NumKiteSites"]), domain=NonNegativeIntegers)# Integer variable to track the number of kite turbines used per site location

    #Variable used to determined the turbine designs used
    Model.W_Wind = Var(range(NumWindDesigns), domain=Binary)# Binary variable to track the wind turbine design used
    Model.W_Wave = Var(range(NumWaveDesigns), domain=Binary)# Binary variable to track the wave turbine design used
    Model.W_Kite = Var(range(NumKiteDesigns), domain=Binary) # Binary variable to track the kite turbine design used

    Model.s     = Var(range(InputDir["NumTransSites"]), domain=Binary)# Binary variable to track the center of the energy collection system
    Model.Delta = Var(range(InputDir["NumTimeSteps"]),  domain=NonNegativeReals) #Curtailment amount


    # #Objective Function
    def objective_rule(Model):   
        EGWind=sum(Model.Y_Wind[i]*InputDir["WindEnergy"][:,i].mean()*InputDir["RatedPowerWindTurbine"][i]  for i in range(InputDir["NumWindSites"])) #Energy generation from wind turbines [MW Avg]
        EGWave=sum(Model.Y_Wave[i]*InputDir["WaveEnergy"][:,i].mean()*InputDir["RatedPowerWaveTurbine"][i]  for i in range(InputDir["NumWaveSites"])) #Energy generation from wave turbines [MW Avg]
        EGKite=sum(Model.Y_Kite[i]*InputDir["KiteEnergy"][:,i].mean()*InputDir["RatedPowerKiteTurbine"][i]  for i in range(InputDir["NumKiteSites"])) #Energy generation from kite turbines [MW Avg]

        TotalCurtailment=sum(Model.Delta[t] for t in range(InputDir["NumTimeSteps"]))/InputDir["NumTimeSteps"] #Average curtailment MW

        Obj=(EGWind+EGWave+EGKite-TotalCurtailment)*24*365.25 # MWh Avg per year
        return Obj

    Model.OBJ = Objective(rule = objective_rule, sense=maximize)

    #Constraints

    #Maximum number of turbines per site location wind
    def MaxTurbinesCell_Wind_rule(Model,i):
        return Model.Y_Wind[i]<=InputDir["MaxNumWindPerSite"][i]
    Model.Turbines_Cell_Wind = Constraint(range(InputDir["NumWindSites"]), rule=MaxTurbinesCell_Wind_rule)

    #Maximum number of turbines per site location kite
    def MaxTurbinesCell_Wave_rule(Model,i):
        return Model.Y_Wave[i]<=InputDir["MaxNumWavePerSite"][i]
    Model.Turbines_Cell_Wave = Constraint(range(InputDir["NumWaveSites"]), rule=MaxTurbinesCell_Wave_rule)

    #Maximum number of turbines per site location kite
    def MaxTurbinesCell_Kite_rule(Model,i):
        return Model.Y_Kite[i]<=InputDir["MaxNumKitePerSite"][i]
    Model.Turbines_Cell_Kite = Constraint(range(InputDir["NumKiteSites"]), rule=MaxTurbinesCell_Kite_rule)

    #Add new constraints to account for multiple technologiges and designs sharing the same region
    #Here pending task


    #Curtailment constraint
    def Curtailment_rule(Model,t):
        EGWind=sum(Model.Y_Wind[i]*InputDir["WindEnergy"][t,i]*InputDir["RatedPowerWindTurbine"][i]  for i in range(InputDir["NumWindSites"])) #Energy generation from wind turbines
        EGWave=sum(Model.Y_Wave[i]*InputDir["WaveEnergy"][t,i]*InputDir["RatedPowerWaveTurbine"][i]  for i in range(InputDir["NumWaveSites"])) #Energy generation from wave turbines
        EGKite=sum(Model.Y_Kite[i]*InputDir["KiteEnergy"][t,i]*InputDir["RatedPowerKiteTurbine"][i]  for i in range(InputDir["NumKiteSites"])) #Energy generation from kite turbines
        
        return -Model.Delta[t]+ EGWind+ EGWave+ EGKite <= InputDir["RatedPowerMWTransmissionMW"]

    Model.Curtailment = Constraint(range(InputDir["NumTimeSteps"]), rule=Curtailment_rule)

    #---Choose center collection system - Start
    #Slect one location for the energy collection system (One location for all technologies)
    #In the future you may be interested in having different locations for each technology or a combination of both

    Model.ChooseOneCircle= Constraint(expr=sum(Model.s[i] for i in range(InputDir["NumTransSites"]))==1)

    #Get the sites that are out of the radious of the center of the collection system
    IdxOutWind=GetIdxOutRadious(InputDir["TransLatLong"], InputDir["WindLatLong"], InputDir["Max_CollectionRadious"])
    IdxOutWave=GetIdxOutRadious(InputDir["TransLatLong"], InputDir["WaveLatLong"], InputDir["Max_CollectionRadious"])
    IdxOutKite=GetIdxOutRadious(InputDir["TransLatLong"], InputDir["KiteLatLong"], InputDir["Max_CollectionRadious"])

    def MaximumRadious(Model,i):  
        SumWind_s=sum(Model.Y_Wind[j] for j in IdxOutWind[i])
        SumWave_s=sum(Model.Y_Wave[j] for j in IdxOutWave[i])
        SumKite_s=sum(Model.Y_Kite[j] for j in IdxOutKite[i])

        return SumWind_s+SumKite_s+SumWave_s<=(1-Model.s[i])*BigM # 1000 is a big M for the maximum total number of turbines installed       

    Model.Maximum_Radious = Constraint(range(InputDir["NumTransSites"]), rule=MaximumRadious)
    #---Choose center collection system - End

    #Track the turbine designs used and limit the number of designs used (Start)
    def TrackDesignsWind_rule(Model,d):  
        IdxVarPartOfDesign=np.where(InputDir["WindDesign"]==d)[0] #Index of Variables associated with the design
        
        #The idea of this constraint is that if the design is not selected, then the sum of the variables associated with the design should be zero
        return sum(Model.Y_Wind[i] for i in IdxVarPartOfDesign)<=Model.W_Wind[d]*BigM 

    def TrackDesignsWave_rule(Model,d):  
        IdxVarPartOfDesign=np.where(InputDir["WaveDesign"]==d)[0] #Index of Variables associated with the design
        
        #The idea of this constraint is that if the design is not selected, then the sum of the variables associated with the design should be zero
        return sum(Model.Y_Wave[i] for i in IdxVarPartOfDesign)<=Model.W_Wave[d]*BigM 

    def TrackDesignsKite_rule(Model,d):  
        IdxVarPartOfDesign=np.where(InputDir["KiteDesign"]==d)[0] #Index of Variables associated with the design
        
        #The idea of this constraint is that if the design is not selected, then the sum of the variables associated with the design should be zero
        return sum(Model.Y_Kite[i] for i in IdxVarPartOfDesign)<=Model.W_Kite[d]*BigM


    Model.TrackDesigns_Wind = Constraint(range(NumWindDesigns), rule=TrackDesignsWind_rule)
    Model.TrackDesigns_Wave = Constraint(range(NumWaveDesigns), rule=TrackDesignsWave_rule)
    Model.TrackDesigns_Kite = Constraint(range(NumKiteDesigns), rule=TrackDesignsKite_rule)

    if NumWindDesigns>0:
        Model.LimitWindDesigns= Constraint(expr=sum(Model.W_Wind[d] for d in range(NumWindDesigns))==InputDir["MaxDesignsWind"])
    
    if NumWaveDesigns>0:
        Model.LimitWaveDesigns= Constraint(expr=sum(Model.W_Wave[d] for d in range(NumWaveDesigns))==InputDir["MaxDesignsWave"])
        
    if NumKiteDesigns>0:
        Model.LimitKiteDesigns= Constraint(expr=sum(Model.W_Kite[d] for d in range(NumKiteDesigns))==InputDir["MaxDesignsKite"])
    #Track the turbine designs used and limit the number of designs used (End)

    #Limit the number of turbines
    
    if NumWindDesigns>0:
        Model.SetLB_Wind= Constraint(expr=sum(Model.Y_Wind[i] for i in range(InputDir["NumWindSites"]))>=InputDir["MinNumWindTurb"])
    
    if NumWaveDesigns>0:
        Model.SetLB_Wave= Constraint(expr=sum(Model.Y_Wave[i] for i in range(InputDir["NumWaveSites"]))>=InputDir["MinNumWaveTurb"])
        
    if NumKiteDesigns>0:
        Model.SetLB_Kite= Constraint(expr=sum(Model.Y_Kite[i] for i in range(InputDir["NumKiteSites"]))>=InputDir["MinNumKiteTrub"])

    ################################### Overlap Constraints ################################### Start
    #Check for overlapping sites and constraint the number of turbines
    #overlaping with wind sites
    if NumWindDesigns>0:
        IdxOverlap_WindWind, AreaOverlap_WindWind, AreaRef1Ref2_WindWind, MaxTurbinesRef1Ref2_WindWind, PercentageOverlap_WindWind=GetOverlaps_Idx_Area(
            InputDir["WindLatLong"], InputDir["WindResolutionKm"], InputDir["WindResolutionDegrees"], InputDir["MaxNumWindPerSite"],
            InputDir["WindLatLong"], InputDir["WindResolutionKm"], InputDir["WindResolutionDegrees"], InputDir["MaxNumWindPerSite"],
            SameTech=1)
        
        IdxOvelap_WindWave, AreaOverlap_WindWave, AreaRef1Ref2_WindWave, MaxTurbinesRef1Ref2_WindWave, PercentageOverlap_WindWave=GetOverlaps_Idx_Area(
            InputDir["WindLatLong"], InputDir["WindResolutionKm"], InputDir["WindResolutionDegrees"], InputDir["MaxNumWindPerSite"],
            InputDir["WaveLatLong"], InputDir["WaveResolutionKm"], InputDir["WaveResolutionDegrees"], InputDir["MaxNumWavePerSite"],
            SameTech=0)
        
        IdxOvelap_WindKite, AreaOverlap_WindKite, AreaRef1Ref2_WindKite, MaxTurbinesRef1Ref2_WindKite, PercentageOverlap_WindKite=GetOverlaps_Idx_Area(
            InputDir["WindLatLong"], InputDir["WindResolutionKm"], InputDir["WindResolutionDegrees"], InputDir["MaxNumWindPerSite"],
            InputDir["KiteLatLong"], InputDir["KiteResolutionKm"], InputDir["KiteResolutionDegrees"], InputDir["MaxNumKitePerSite"],
            SameTech=0)
        
        #Wind sites with some overlap
        IdxOvelap_UniqueWindIdx=np.unique(np.concatenate((IdxOverlap_WindWind,IdxOvelap_WindWave,IdxOvelap_WindKite))[:,0])

        def TrackOverlaps_Wind_rule(Model,i):
            
            
            IdxWindWind_f=IdxOverlap_WindWind[IdxOverlap_WindWind[:,0]==i, 1]
            IdxWindWave_f=IdxOvelap_WindWave[IdxOvelap_WindWave[:,0]==i,   1]
            IdxWindKite_f=IdxOvelap_WindKite[IdxOvelap_WindKite[:,0]==i,   1]
            
            AreaWindWind_f=AreaRef1Ref2_WindWind[IdxOverlap_WindWind[:,0]==i]
            AreaWindWave_f=AreaRef1Ref2_WindWave[IdxOvelap_WindWave[:,0]==i]
            AreaWindKite_f=AreaRef1Ref2_WindKite[IdxOvelap_WindKite[:,0]==i]
            
            MaxTurbinesWind_f=MaxTurbinesRef1Ref2_WindWind[IdxOverlap_WindWind[:,0]==i]
            MaxTurbinesWave_f=MaxTurbinesRef1Ref2_WindWave[IdxOvelap_WindWave[:,0]==i]
            MaxTurbinesKite_f=MaxTurbinesRef1Ref2_WindKite[IdxOvelap_WindKite[:,0]==i]
            
            PercentageOverlap_WindWind_f=PercentageOverlap_WindWind[IdxOverlap_WindWind[:,0]==i]
            PercentageOverlap_WindWave_f=PercentageOverlap_WindWave[IdxOvelap_WindWave[:,0]==i]
            PercentageOverlap_WindKite_f=PercentageOverlap_WindKite[IdxOvelap_WindKite[:,0]==i]
            
            #Compute the overlaped area and estimate the equivalent number of turbines that cannot be installed on the ith location anymore
            Expression=Model.Y_Wind[i] <= InputDir["MaxNumWindPerSite"][i]\
                    -sum((AreaWindWind_f[k,1]/MaxTurbinesWind_f[k,1]*Model.Y_Wind[IdxWindWind_f[k]])*PercentageOverlap_WindWind_f[k]\
                        *MaxTurbinesWind_f[k,0]/AreaWindWind_f[k,0] for k in range(len(IdxWindWind_f)))\
                    -sum((AreaWindWave_f[k,1]/MaxTurbinesWave_f[k,1]*Model.Y_Wave[IdxWindWave_f[k]])*PercentageOverlap_WindWave_f[k]\
                        *MaxTurbinesWave_f[k,0]/AreaWindWave_f[k,0] for k in range(len(IdxWindWave_f)))\
                    -sum((AreaWindKite_f[k,1]/MaxTurbinesKite_f[k,1]*Model.Y_Kite[IdxWindKite_f[k]])*PercentageOverlap_WindKite_f[k]\
                        *MaxTurbinesKite_f[k,0]/AreaWindKite_f[k,0] for k in range(len(IdxWindKite_f)))
                    
            return Expression
            
            
        Model.OvelapWind_ALL= Constraint(list(IdxOvelap_UniqueWindIdx), rule=TrackOverlaps_Wind_rule)
        
    #overlaping with wave sites
    if NumWaveDesigns>0:
        
        IdxOvelap_WaveWind, AreaOverlap_WaveWind, AreaRef1Ref2_WaveWind, MaxTurbinesRef1Ref2_WaveWind, PercentageOverlap_WaveWind=GetOverlaps_Idx_Area(
            InputDir["WaveLatLong"], InputDir["WaveResolutionKm"], InputDir["WaveResolutionDegrees"], InputDir["MaxNumWavePerSite"],
            InputDir["WindLatLong"], InputDir["WindResolutionKm"], InputDir["WindResolutionDegrees"], InputDir["MaxNumWindPerSite"],
            SameTech=0)
        
        
        IdxOverlap_WaveWave, AreaOverlap_WaveWave, AreaRef1Ref2_WaveWave, MaxTurbinesRef1Ref2_WaveWave, PercentageOverlap_WaveWave=GetOverlaps_Idx_Area(
            InputDir["WaveLatLong"], InputDir["WaveResolutionKm"], InputDir["WaveResolutionDegrees"], InputDir["MaxNumWavePerSite"],
            InputDir["WaveLatLong"], InputDir["WaveResolutionKm"], InputDir["WaveResolutionDegrees"], InputDir["MaxNumWavePerSite"],
            SameTech=1)
        
        IdxOverlap_WaveKite, AreaOverlap_WaveKite, AreaRef1Ref2_WaveKite, MaxTurbinesRef1Ref2_WaveKite, PercentageOverlap_WaveKite=GetOverlaps_Idx_Area(
            InputDir["WaveLatLong"], InputDir["WaveResolutionKm"], InputDir["WaveResolutionDegrees"], InputDir["MaxNumWavePerSite"],
            InputDir["KiteLatLong"], InputDir["KiteResolutionKm"], InputDir["KiteResolutionDegrees"], InputDir["MaxNumKitePerSite"],
            SameTech=0)
        
        #Wave sites with some overlap
        IdxOvelap_UniqueWaveIdx=np.unique(np.concatenate((IdxOvelap_WaveWind,IdxOverlap_WaveWave,IdxOverlap_WaveKite))[:,0])
        
        def TrackOverlaps_Wave_rule(Model,i):
                
                
            IdxWaveWind_f=IdxOvelap_WaveWind[IdxOvelap_WaveWind[:,0]==i,   1]
            IdxWaveWave_f=IdxOverlap_WaveWave[IdxOverlap_WaveWave[:,0]==i, 1]
            IdxWaveKite_f=IdxOverlap_WaveKite[IdxOverlap_WaveKite[:,0]==i, 1]
            
            AreaWaveWind_f=AreaRef1Ref2_WaveWind[IdxOvelap_WaveWind[:,0]==i]
            AreaWaveWave_f=AreaRef1Ref2_WaveWave[IdxOverlap_WaveWave[:,0]==i]
            AreaWaveKite_f=AreaRef1Ref2_WaveKite[IdxOverlap_WaveKite[:,0]==i]
            
            MaxTurbinesWaveWind_f=MaxTurbinesRef1Ref2_WaveWind[IdxOvelap_WaveWind[:,0]==i]
            MaxTurbinesWaveWave_f=MaxTurbinesRef1Ref2_WaveWave[IdxOverlap_WaveWave[:,0]==i]
            MaxTurbinesWaveKite_f=MaxTurbinesRef1Ref2_WaveKite[IdxOverlap_WaveKite[:,0]==i]
            
            PercentageOverlap_WaveWind_f=PercentageOverlap_WaveWind[IdxOvelap_WaveWind[:,0]==i]
            PercentageOverlap_WaveWave_f=PercentageOverlap_WaveWave[IdxOverlap_WaveWave[:,0]==i]
            PercentageOverlap_WaveKite_f=PercentageOverlap_WaveKite[IdxOverlap_WaveKite[:,0]==i]
            
            #Compute the overlaped area and estimate the equivalent number of turbines that cannot be installed on the ith location anymore
            Expression=Model.Y_Wave[i] <= InputDir["MaxNumWavePerSite"][i]\
                    -sum((AreaWaveWind_f[k,1]/MaxTurbinesWaveWind_f[k,1]*Model.Y_Wind[IdxWaveWind_f[k]])*PercentageOverlap_WaveWind_f[k]\
                        *MaxTurbinesWaveWind_f[k,0]/AreaWaveWind_f[k,0] for k in range(len(IdxWaveWind_f)))\
                    -sum((AreaWaveWave_f[k,1]/MaxTurbinesWaveWave_f[k,1]*Model.Y_Wave[IdxWaveWave_f[k]])*PercentageOverlap_WaveWave_f[k]\
                        *MaxTurbinesWaveWave_f[k,0]/AreaWaveWave_f[k,0] for k in range(len(IdxWaveWave_f)))\
                    -sum((AreaWaveKite_f[k,1]/MaxTurbinesWaveKite_f[k,1]*Model.Y_Kite[IdxWaveKite_f[k]])*PercentageOverlap_WaveKite_f[k]\
                        *MaxTurbinesWaveKite_f[k,0]/AreaWaveKite_f[k,0] for k in range(len(IdxWaveKite_f)))
            return Expression

        Model.OverlapWave_ALL= Constraint(list(IdxOvelap_UniqueWaveIdx), rule=TrackOverlaps_Wave_rule)
        
    #overlaping with kite sites
    if NumKiteDesigns>0:
        
        IdxOvelap_KiteWind, AreaOverlap_KiteWind, AreaRef1Ref2_KiteWind, MaxTurbinesRef1Ref2_KiteWind, PercentageOverlap_KiteWind=GetOverlaps_Idx_Area(
            InputDir["KiteLatLong"], InputDir["KiteResolutionKm"], InputDir["KiteResolutionDegrees"], InputDir["MaxNumKitePerSite"],
            InputDir["WindLatLong"], InputDir["WindResolutionKm"], InputDir["WindResolutionDegrees"], InputDir["MaxNumWindPerSite"],
            SameTech=0)
        
        IdxOvelap_KiteWave, AreaOverlap_KiteWave, AreaRef1Ref2_KiteWave, MaxTurbinesRef1Ref2_KiteWave, PercentageOverlap_KiteWave=GetOverlaps_Idx_Area(
            InputDir["KiteLatLong"], InputDir["KiteResolutionKm"], InputDir["KiteResolutionDegrees"], InputDir["MaxNumKitePerSite"],
            InputDir["WaveLatLong"], InputDir["WaveResolutionKm"], InputDir["WaveResolutionDegrees"], InputDir["MaxNumWavePerSite"],
            SameTech=0)
        
        IdxOverlap_KiteKite, AreaOverlap_KiteKite, AreaRef1Ref2_KiteKite, MaxTurbinesRef1Ref2_KiteKite, PercentageOverlap_KiteKite=GetOverlaps_Idx_Area(
            InputDir["KiteLatLong"], InputDir["KiteResolutionKm"], InputDir["KiteResolutionDegrees"], InputDir["MaxNumKitePerSite"],
            InputDir["KiteLatLong"], InputDir["KiteResolutionKm"], InputDir["KiteResolutionDegrees"], InputDir["MaxNumKitePerSite"],
            SameTech=1)
        
        #Kite sites with some overlap
        IdxOvelap_UniqueKiteIdx=np.unique(np.concatenate((IdxOvelap_KiteWind,IdxOvelap_KiteWave,IdxOverlap_KiteKite))[:,0])
        
        def TrackOverlaps_Kite_rule(Model,i):
                
                
            IdxKiteWind_f=IdxOvelap_KiteWind[IdxOvelap_KiteWind[:,0]==i,   1]
            IdxKiteWave_f=IdxOvelap_KiteWave[IdxOvelap_KiteWave[:,0]==i,   1]
            IdxKiteKite_f=IdxOverlap_KiteKite[IdxOverlap_KiteKite[:,0]==i, 1]
            
            AreaKiteWind_f=AreaRef1Ref2_KiteWind[IdxOvelap_KiteWind[:,0]==i]
            AreaKiteWave_f=AreaRef1Ref2_KiteWave[IdxOvelap_KiteWave[:,0]==i]
            AreaKiteKite_f=AreaRef1Ref2_KiteKite[IdxOverlap_KiteKite[:,0]==i]
            
            MaxTurbinesKiteWind_f=MaxTurbinesRef1Ref2_KiteWind[IdxOvelap_KiteWind[:,0]==i]
            MaxTurbinesKiteWave_f=MaxTurbinesRef1Ref2_KiteWave[IdxOvelap_KiteWave[:,0]==i]
            MaxTurbinesKiteKite_f=MaxTurbinesRef1Ref2_KiteKite[IdxOverlap_KiteKite[:,0]==i]
            
            PercentageOverlap_KiteWind_f=PercentageOverlap_KiteWind[IdxOvelap_KiteWind[:,0]==i]
            PercentageOverlap_KiteWave_f=PercentageOverlap_KiteWave[IdxOvelap_KiteWave[:,0]==i]
            PercentageOverlap_KiteKite_f=PercentageOverlap_KiteKite[IdxOverlap_KiteKite[:,0]==i]
            
            #Compute the overlaped area and estimate the equivalent number of turbines that cannot be installed on the ith location anymore
            Expression=Model.Y_Kite[i] <= InputDir["MaxNumKitePerSite"][i]\
                    -sum((AreaKiteWind_f[k,1]/MaxTurbinesKiteWind_f[k,1]*Model.Y_Wind[IdxKiteWind_f[k]])*PercentageOverlap_KiteWind_f[k]\
                        *MaxTurbinesKiteWind_f[k,0]/AreaKiteWind_f[k,0] for k in range(len(IdxKiteWind_f)))\
                    -sum((AreaKiteWave_f[k,1]/MaxTurbinesKiteWave_f[k,1]*Model.Y_Wave[IdxKiteWave_f[k]])*PercentageOverlap_KiteWave_f[k]\
                        *MaxTurbinesKiteWave_f[k,0]/AreaKiteWave_f[k,0] for k in range(len(IdxKiteWave_f)))\
                    -sum((AreaKiteKite_f[k,1]/MaxTurbinesKiteKite_f[k,1]*Model.Y_Kite[IdxKiteKite_f[k]])*PercentageOverlap_KiteKite_f[k]\
                        *MaxTurbinesKiteKite_f[k,0]/AreaKiteKite_f[k,0] for k in range(len(IdxKiteKite_f)))
            return Expression
        
        Model.OverlapKite_ALL= Constraint(list(IdxOvelap_UniqueKiteIdx), rule=TrackOverlaps_Kite_rule)
    ####
    ################################### Overlap Constraints ################################### End
    
    # #LCOE Target (Attached later on the LCOE iterator)
    # def LCOETarget(Model, LCOE_Max):  
    #     EGWind=sum(Model.Y_Wind[i]*InputDir["WindEnergy"][:,i].mean()*InputDir["RatedPowerWindTurbine"][i]  for i in range(InputDir["NumWindSites"])) #Energy generation from wind turbines [MW Avg]
    #     EGWave=sum(Model.Y_Wave[i]*InputDir["WaveEnergy"][:,i].mean()*InputDir["RatedPowerWaveTurbine"][i]  for i in range(InputDir["NumWaveSites"])) #Energy generation from wave turbines [MW Avg]
    #     EGKite=sum(Model.Y_Kite[i]*InputDir["KiteEnergy"][:,i].mean()*InputDir["RatedPowerKiteTurbine"][i]  for i in range(InputDir["NumKiteSites"])) #Energy generation from kite turbines [MW Avg]

    #     TotalCurtailment=sum(Model.Delta[t] for t in range(InputDir["NumTimeSteps"]))/InputDir["NumTimeSteps"] #Average curtailment MW

    #     MWhYear=(EGWind+EGWave+EGKite-TotalCurtailment)*24*365.25 # MWh Avg per year


    #     Cost_Wind=sum(Model.Y_Wind[i]*InputDir["AnnualizedCostWind"][i]  for i in range(InputDir["NumWindSites"]))
    #     Cost_Wave=sum(Model.Y_Wave[i]*InputDir["AnnualizedCostWave"][i]  for i in range(InputDir["NumWaveSites"]))
    #     Cost_Kite=sum(Model.Y_Kite[i]*InputDir["AnnualizedCostKite"][i]  for i in range(InputDir["NumKiteSites"]))
        
        
    #     Cost_Transmission=sum(Model.s[i]*InputDir["AnnualizedCostTransmission"][i] for i in Model.SiteTrs)
        
    #     TotalCost=Cost_Wind+Cost_Wave+Cost_Kite+Cost_Transmission
        

    #     return TotalCost<=LCOE_Max*MWhYear  

    return Model, InputDir


def SolvePortOpt_MaxGen_LCOE_Iterator(PathWindDesigns, PathWaveDesigns, PathKiteDesigns, PathTransmissionDesign, LCOE_RANGE\
    ,Max_CollectionRadious,MaxDesignsWind, MaxDesingsWave, MaxDesingsKite,MinNumWindTurb,MinNumWaveTurb,MinNumKiteTrub\
    ,ReadMe,SavePath=None):

    #Create inputs and main model structure
    Model, InputDir=SolvePortOpt_MaxGen_Model(PathWindDesigns, PathWaveDesigns, PathKiteDesigns, PathTransmissionDesign, LCOE_RANGE\
        ,Max_CollectionRadious,MaxDesignsWind, MaxDesingsWave, MaxDesingsKite,MinNumWindTurb,MinNumWaveTurb,MinNumKiteTrub)

    opt = SolverFactory('gurobi', solver_io="python")
    opt.options['mipgap'] = 0.05

    #LCOE Target
    def LCOETarget_rule(Model, LCOE_Max):  
        EGWind=sum(Model.Y_Wind[i]*InputDir["WindEnergy"][:,i].mean()*InputDir["RatedPowerWindTurbine"][i]  for i in range(InputDir["NumWindSites"])) #Energy generation from wind turbines [MW Avg]
        EGWave=sum(Model.Y_Wave[i]*InputDir["WaveEnergy"][:,i].mean()*InputDir["RatedPowerWaveTurbine"][i]  for i in range(InputDir["NumWaveSites"])) #Energy generation from wave turbines [MW Avg]
        EGKite=sum(Model.Y_Kite[i]*InputDir["KiteEnergy"][:,i].mean()*InputDir["RatedPowerKiteTurbine"][i]  for i in range(InputDir["NumKiteSites"])) #Energy generation from kite turbines [MW Avg]

        TotalCurtailment=sum(Model.Delta[t] for t in range(InputDir["NumTimeSteps"]))/InputDir["NumTimeSteps"] #Average curtailment MW

        MWhYear=(EGWind+EGWave+EGKite-TotalCurtailment)*24*365.25 # MWh Avg per year


        Cost_Wind=sum(Model.Y_Wind[i]*InputDir["AnnualizedCostWind"][i]  for i in range(InputDir["NumWindSites"]))
        Cost_Wave=sum(Model.Y_Wave[i]*InputDir["AnnualizedCostWave"][i]  for i in range(InputDir["NumWaveSites"]))
        Cost_Kite=sum(Model.Y_Kite[i]*InputDir["AnnualizedCostKite"][i]  for i in range(InputDir["NumKiteSites"]))
        
        
        Cost_Transmission=sum(Model.s[i]*InputDir["AnnualizedCostTransmission"][i] for i in range(InputDir["NumTransSites"]))
        
        TotalCost=Cost_Wind+Cost_Wave+Cost_Kite+Cost_Transmission #M$
        TotalCost=TotalCost*10**6 #USD


        return TotalCost<=LCOE_Max*MWhYear  

    SaveFeasibility, Save_LCOETarget, Save_LCOE_Achieved, SaveTotalMWAvg = list(), list(), list(), list()
    Save_Y_Wind, Save_Y_Wave, Save_Y_Kite, Save_W_Wind, Save_W_Wave, Save_W_Kite, Save_s, Save_Delta = list(), list(), list(), list(), list(), list(), list(), list()

    LowestLCOE=10**10
    for LCOETarget in tqdm(InputDir["LCOE_RANGE"]):
        
        #Skip based on the algorithm progress, avoid repeating the same LCOE*
        if LCOETarget<LowestLCOE:    
            Bypass=0
            
            #Upperbound For the LCOE Activate Constraint
            LCOETarget_rule_tmp=LCOETarget_rule(Model,LCOETarget)
            Model.LCOE_Target = Constraint(rule=LCOETarget_rule_tmp)
            print("Running Model With LCOE= %.2f" % LCOETarget)
            
            try:
                results=opt.solve(Model, tee=False)
            except:
                Bypass=1
                Model.del_component(Model.LCOE_Target)  
        
            if Bypass==0:
                if (results.solver.status == SolverStatus.ok) and (results.solver.termination_condition == TerminationCondition.optimal):
                    SaveFeasibility.append(1)
                    Save_LCOETarget.append(LCOETarget)
                    
                    Optimal_Y_Wind=np.array([Model.Y_Wind[i].value for i in range(InputDir["NumWindSites"])])
                    Optimal_Y_Wave=np.array([Model.Y_Wave[i].value for i in range(InputDir["NumWaveSites"])])
                    Optimal_Y_Kite=np.array([Model.Y_Kite[i].value for i in range(InputDir["NumKiteSites"])])
                    Optimal_W_Wind=np.array([Model.W_Wind[i].value for i in range(len(Model.W_Wind))])
                    Optimal_W_Wave=np.array([Model.W_Wave[i].value for i in range(len(Model.W_Wave))])
                    Optimal_W_Kite=np.array([Model.W_Kite[i].value for i in range(len(Model.W_Kite))])
                    Optimal_s=np.array([Model.s[i].value for i in range(InputDir["NumTransSites"])])
                    Optimal_Delta=np.array([Model.Delta[i].value for i in range(InputDir["NumTimeSteps"])])
                    
                    Save_Y_Wind.append(Optimal_Y_Wind)
                    Save_Y_Wave.append(Optimal_Y_Wave)
                    Save_Y_Kite.append(Optimal_Y_Kite)
                    Save_W_Wind.append(Optimal_W_Wind)
                    Save_W_Wave.append(Optimal_W_Wave)
                    Save_W_Kite.append(Optimal_W_Kite)
                    Save_s.append(Optimal_s)
                    Save_Delta.append(Optimal_Delta)
                    

                    #Current LCOE
                    EGWind=sum(Optimal_Y_Wind[i]*InputDir["WindEnergy"][:,i].mean()*InputDir["RatedPowerWindTurbine"][i]  for i in range(InputDir["NumWindSites"])) #Energy generation from wind turbines [MW Avg]
                    EGWave=sum(Optimal_Y_Wave[i]*InputDir["WaveEnergy"][:,i].mean()*InputDir["RatedPowerWaveTurbine"][i]  for i in range(InputDir["NumWaveSites"])) #Energy generation from wave turbines [MW Avg]
                    EGKite=sum(Optimal_Y_Kite[i]*InputDir["KiteEnergy"][:,i].mean()*InputDir["RatedPowerKiteTurbine"][i]  for i in range(InputDir["NumKiteSites"])) #Energy generation from kite turbines [MW Avg]

                    TotalCurtailment=sum(Optimal_Delta[t] for t in range(InputDir["NumTimeSteps"]))/InputDir["NumTimeSteps"] #Average curtailment MW

                    MWhYear=(EGWind+EGWave+EGKite-TotalCurtailment)*24*365.25 # MWh Avg per year


                    Cost_Wind=sum(Optimal_Y_Wind[i]*InputDir["AnnualizedCostWind"][i]  for i in range(InputDir["NumWindSites"]))
                    Cost_Wave=sum(Optimal_Y_Wave[i]*InputDir["AnnualizedCostWave"][i]  for i in range(InputDir["NumWaveSites"]))
                    Cost_Kite=sum(Optimal_Y_Kite[i]*InputDir["AnnualizedCostKite"][i]  for i in range(InputDir["NumKiteSites"]))
                    
                    
                    Cost_Transmission=sum(Optimal_s[i]*InputDir["AnnualizedCostTransmission"][i] for i in range(InputDir["NumTransSites"]))
                    
                    TotalCost=Cost_Wind+Cost_Wave+Cost_Kite+Cost_Transmission #M$
                    TotalCost=TotalCost*10**6 #USD
                
                    CurrentLCOE=TotalCost/MWhYear
                    LowestLCOE=CurrentLCOE
                    
                    Save_LCOE_Achieved.append(CurrentLCOE)
                    SaveTotalMWAvg.append(MWhYear)
                    
                    
                    print("LCOE OPT: %.2f,\n MW Wind: %.2f,\nMW Wave: %.2f,\nMW Kite: %.2f,\nMW Curtailment: %.2f,\nMW Total: %.2f\n" % (CurrentLCOE,EGWind,EGWave,EGKite,TotalCurtailment,EGWind+EGWave+EGKite-TotalCurtailment))
                    
                    #Delete constraint for its modification in the next step of the for loop
                    Model.del_component(Model.LCOE_Target)

                else:# Something else is wrong
                    Model.del_component(Model.LCOE_Target)
                    SaveFeasibility.append(0)
                    Save_LCOETarget.append(None)
                    Save_LCOE_Achieved.append(None)
                    SaveTotalMWAvg.append(None)   
                    
                    Save_Y_Wind.append(None)
                    Save_Y_Wave.append(None)
                    Save_Y_Kite.append(None)
                    Save_W_Wind.append(None)
                    Save_W_Wave.append(None)
                    Save_W_Kite.append(None)
                    Save_s.append(None)
                    Save_Delta.append(None)
                    break

    #Save Results
    if SavePath!=None:
        np.savez(SavePath, 
                ReadMe=ReadMe,
                #Model Inputs
                PathWindDesigns=PathWindDesigns,
                PathWaveDesigns=PathWaveDesigns,
                PathKiteDesigns=PathKiteDesigns,
                PathTransmissionDesign=PathTransmissionDesign,
                LCOE_RANGE=LCOE_RANGE,
                Max_CollectionRadious=Max_CollectionRadious,
                MaxDesignsWind=MaxDesignsWind,
                MaxDesingsWave=MaxDesingsWave,
                MaxDesingsKite=MaxDesingsKite,
                MinNumWindTurb=MinNumWindTurb,
                MinNumWaveTurb=MinNumWaveTurb,
                MinNumKiteTrub=MinNumKiteTrub,

                #Model Outputs
                SaveFeasibility=SaveFeasibility,
                Save_LCOETarget=Save_LCOETarget,
                Save_LCOE_Achieved=Save_LCOE_Achieved,
                SaveTotalMWAvg=SaveTotalMWAvg,
                
                Save_Y_Wind=Save_Y_Wind,
                Save_Y_Wave=Save_Y_Wave,
                Save_Y_Kite=Save_Y_Kite,
                Save_W_Wind=Save_W_Wind,
                Save_W_Wave=Save_W_Wave,
                Save_W_Kite=Save_W_Kite,
                Save_s=Save_s,
                Save_Delta=Save_Delta,
                
                )
