#Compute costs for wave energy
import numpy as np
def ComputeAnnualCost(DistanceShoreToPlatform, TotalMorringChainLength, LengthCable34kV=None):
    
    #TotalMorringChainLength: Not used in the current code, but kept here in case we decide to do
    #a more detailed modeling for the wave technology
    
    #Default Values
    NumTurbines=100
    RatedPower=1.5#Rated power for one device [MW]

    
    if LengthCable34kV==None:
        LengthCable34kV=NumTurbines*1# (No platform was considered)[km]
    
    FCR=11.3/100 #Fixed Charge Rate
      
    ##---------------------#Wave energy CAPEX----------------------------------------------------------##    
    #Manufacturing cost
    CAPEX_OC = 1269.6 + 1.003*DistanceShoreToPlatform + 0.186*LengthCable34kV

    OPEX_OC=30.4868 + 0.025*DistanceShoreToPlatform + 4.5*(10**-3)*LengthCable34kV
    ##---------------------#Annualized Cost----------------------------------------------------------##      
    AnnCost=(CAPEX_OC)*FCR + (OPEX_OC)
    
    return AnnCost, CAPEX_OC, OPEX_OC


AnnualizedCostWave=[]
Capex=[]
Opex=[]
Wave_Data=np.load('WaveEnergy_Pelamis_2009_2013.npz',allow_pickle=True)
DepthSites=Wave_Data['DepthSites']
ShoreDistance=Wave_Data['ShoreDistance']

for i in range(len(DepthSites)):
    TotalMorringChainLength=DepthSites[i]*100/1000#[km]
    DistanceShoreToPlatform=ShoreDistance[i]
    
    UnitCost, CAPEX_OC, OPEX_OC=ComputeAnnualCost(DistanceShoreToPlatform, TotalMorringChainLength)
    UnitCost=UnitCost/100
    AnnualizedCostWave.append(UnitCost)
    Capex.append(CAPEX_OC/100)
    Opex.append(OPEX_OC/100)
    
    
AnnualizedCostWave=np.array(AnnualizedCostWave)
    
np.savez('WaveEnergy_Pelamis_2009_2013.npz',
         Energy_pu=Wave_Data['Energy_pu'], DepthSites=Wave_Data['DepthSites'],RatedPower=Wave_Data['RatedPower'],\
         LatLong=Wave_Data['LatLong'],\
         ShoreDistance=Wave_Data['ShoreDistance'],AnnualizedCostWave=AnnualizedCostWave,\
             Capex=Capex, Opex=Opex)                  
