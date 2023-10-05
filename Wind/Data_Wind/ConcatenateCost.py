#Compute costs for wind energy
import numpy as np
#import matplotlib.colors as clrs
#import matplotlib.pyplot as plt

WindEnergyFile='WindEnergyNREL_100m_Haliade150_6MW.npz'
WindEnergyData=np.load(WindEnergyFile)

Depth_Data=WindEnergyData['Depth']
DistanceToShore_Data=WindEnergyData['DistanceToShore']
#plt.scatter(NREL_WaterDepth,NREL_ShoreDistance,c='b',  s=5)
#plt.scatter(Depth_Data,DistanceToShore_Data,c='r' , s=5)

WaterDepth=100 #[m]
DistanceToShore=30 #[km]
FCR=8.8/100 # Factor of Capital Return(NREL: https://www.nrel.gov/docs/fy18osti/72167.pdf)


#def ComputeAnnualCost(DistanceShoreToPlatform, TotalMorringChainLength,LengthCable34kV=None):
RatedPowerOfTurbine=6 #[MW]
NREL_WaterDepth=np.array([18,22,24,29,31,144,159,157,148,107,375,467,663,432,468]) #[m]
NREL_ShoreDistance=np.array([27,29,33,57,65,38,45,46,64,101,116,116,166,147,133]) #[km]

NREL_CAPEX=np.array([3361,3475,3660,4001,4633,4602,4661,4710,4936,5209,5320,5331,5785,6145,6599]) #[$/kW]
NREL_OPEX=np.array([105,106,109,112,107,77,78,78,84,91,98,97,92,87,83]) #[$/(kW*Year)]

NREL_CAPEX=NREL_CAPEX*RatedPowerOfTurbine*(10**-3) #[M$]
NREL_OPEX=NREL_OPEX*RatedPowerOfTurbine*(10**-3) #[$/Year]


Diff_Depth=np.reshape(Depth_Data,(len(Depth_Data),1))-np.reshape(NREL_WaterDepth,(1,len(NREL_WaterDepth)))
Diff_Distance=np.reshape(DistanceToShore_Data,(len(DistanceToShore_Data),1))-np.reshape(NREL_ShoreDistance,(1,len(NREL_ShoreDistance)))

EuclidianDistance=(Diff_Depth**2+Diff_Distance**2)
Idx_NRELTurbine=np.argmin(EuclidianDistance,axis=1)

AnnualizedCostWind=NREL_CAPEX[Idx_NRELTurbine]*FCR + NREL_OPEX[Idx_NRELTurbine]

np.savez(WindEnergyFile,ReadMe=WindEnergyData['ReadMe'],WindEnergy=WindEnergyData['WindEnergy'],\
         RatedPower=WindEnergyData['RatedPower'],LatLong=WindEnergyData['LatLong'],\
         Depth=WindEnergyData['Depth'],DistanceToShore=WindEnergyData['DistanceToShore'],AnnualizedCostWind=AnnualizedCostWind,\
             Capex=NREL_CAPEX[Idx_NRELTurbine], Opex=NREL_OPEX[Idx_NRELTurbine])