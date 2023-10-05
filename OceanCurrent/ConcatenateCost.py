#Concatenate costs in the .npz file containing information for ocean current energy
import numpy as np
import pandas as pd
import geoplot as gplt
import geopandas as gpd
import csv
import matplotlib.colors as clrs
import matplotlib.pyplot as plt

def ComputeAnnualCost(DistanceShoreToPlatform, TotalMorringChainLength, LengthCable34kV=None):
    
    #Default Values
    NumTurbines=50
    RatedPower=4#Rated power for one device [MW]
    #DistanceShoreToPlatform=150# [km]
    #TotalMorringChainLength=NumTurbines*1000# [km]
    
    if LengthCable34kV==None:
        LengthCable34kV=NumTurbines*1.5 #[km]
    
    FCR=11.3/100 #Fixed Charge Rate
    
    ##---------------------#Transmission Cost------------------------------------------------------------##
    #For a TL with Capacity for 50 Devices
    if DistanceShoreToPlatform<66:
        CAPEX_Transmission=19.31 + 1.28*DistanceShoreToPlatform
        
    else:
        CAPEX_Transmission=69.98 + 0.50*DistanceShoreToPlatform  
    
    OPEX_Transmission=2.5/100*CAPEX_Transmission
    
    
    ##---------------------#Ocean Current CAPEX----------------------------------------------------------##
    
    #Development costs
    DevelopmentCost=-13.5643 + 10.37*np.log(NumTurbines*RatedPower + 3.7)# Regression from SANDIA report 
    
    #Infrastructure
    Infrastructure= 0.4245*LengthCable34kV+ 30*np.ceil(NumTurbines/40)
    
    #Mooring/Foundation
    Mooring_Foundation= 1.5317*NumTurbines+ 0.344*TotalMorringChainLength + 0.4034
    
    #Device Structural Components
    StructuralComponentsCost=1.23*(-3.879*10**(-3)*NumTurbines**2 + 5.0876*NumTurbines + 5.613)
    
    #Power Take Off
    PowerTakeOff=-18.787*10**(-3)*NumTurbines**2 + 12.3961*NumTurbines + 16.0626 
    
    ##### Subsystem Integration and Profit Margin
    Integration_ProfitMargin=0.1*(StructuralComponentsCost + PowerTakeOff)
    
    #Instalation
    Instalation=1.271375*NumTurbines + 6.9445*10**(-3)*NumTurbines*DistanceShoreToPlatform \
    + 18.2521*10**(-3)*DistanceShoreToPlatform + 0.11357*LengthCable34kV + 8.73053
    
    #Contingency
    Contingency=0.1*(Infrastructure+Mooring_Foundation\
                +StructuralComponentsCost+PowerTakeOff+Integration_ProfitMargin+Instalation)
    
    CAPEX_OC=DevelopmentCost+Infrastructure+Mooring_Foundation+StructuralComponentsCost+PowerTakeOff\
    +Integration_ProfitMargin+Instalation+Contingency
    
    ##---------------------#Ocean Current OPEX----------------------------------------------------------##
    
    #Insurance
    Insurance=0.01*CAPEX_OC
    
    #Environmental and Regulatory Compliance
    Environmental_RegulatoryCompliance=2.036 
    
    #Offshore Operations
    MarineOperations=5.128*10**(-6)*NumTurbines*(17925 + 150*DistanceShoreToPlatform)
    
    #Shore Operations
    ShoresideOperations=15.555*10**(-3)*NumTurbines + 0.2205 
    
    #Replacement Parts
    ReplacementParts= 0.86/100*PowerTakeOff
    
    #Consumables
    Consumables=17.494 *10**(-3)*NumTurbines
    
    OPEX_OC=Insurance+Environmental_RegulatoryCompliance+ MarineOperations +ShoresideOperations\
        + ReplacementParts+ Consumables
       
    ##---------------------#Annualized Cost----------------------------------------------------------##      
    AnnCost=(CAPEX_Transmission+ CAPEX_OC)*FCR + (OPEX_Transmission+ OPEX_OC)
    Capex=CAPEX_Transmission+ CAPEX_OC
    Opex=(OPEX_Transmission+ OPEX_OC)
    return AnnCost, Capex, Opex


AnnualizedCostOcean=[]
Capex=[]
Opex=[]

Ocean_Data=np.load('OceanCurrentEnergyRM4.npz',allow_pickle=True)
DepthSites=Ocean_Data['DepthSites']
ShoreDistance=Ocean_Data['ShoreDistance']
NumOfSites=len(DepthSites)

for i in range(NumOfSites):
    TotalMorringChainLength=DepthSites[i]*50/1000#[km]
    DistanceShoreToPlatform=ShoreDistance[i]
    
    UnitCost,CAPEX, OPEX=ComputeAnnualCost(DistanceShoreToPlatform, TotalMorringChainLength)
    UnitCost=UnitCost/50
    AnnualizedCostOcean.append(UnitCost)
    Capex.append(CAPEX/50)
    Opex.append(OPEX/50)
    
AnnualizedCostOcean=np.array(AnnualizedCostOcean)
    
ReadMe='\
CurrentEnergy_pu: Ocean current pu energy\m\
CurrentSpeed: Ocean current speed\n\
LatLong: Latitude,Logitude data\n\
DepthSites: Depth in each site location\n\
ShoreDistance: Smallest distance from site location to shore\n\
RatedPower: Rated power of the turbine model studied\n\
OceanDateTime: Time of the energy estimation\n\
AnnualizedCostOcean: Annualized cost per turbine unit in each site location'

np.savez('OceanCurrentEnergyRM4.npz',ReadMe=ReadMe,CurrentSpeed=Ocean_Data['CurrentSpeed'],\
         CurrentEnergy_pu=Ocean_Data['CurrentEnergy_pu'], RatedPower=Ocean_Data['RatedPower'],DepthSites=Ocean_Data['DepthSites'],\
         LatLong=Ocean_Data['LatLong'],OceanDateTime=Ocean_Data['OceanDateTime'],\
         ShoreDistance=Ocean_Data['ShoreDistance'],AnnualizedCostOcean=AnnualizedCostOcean,\
            Capex=Capex, Opex=Opex)             