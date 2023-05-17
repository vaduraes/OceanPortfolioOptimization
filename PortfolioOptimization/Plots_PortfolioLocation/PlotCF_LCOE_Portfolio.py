#Plot Deployments

import numpy as np
import datetime as dt
import numpy as np
import geopandas as gpd
import datetime as dt
import csv
import matplotlib.colors as clrs
import matplotlib.pyplot as plt

def DistanceBetweenLatLong(LatLong1, LatLong2):
    LatLong1=LatLong1*2*np.pi/360
    LatLong2=LatLong2*2*np.pi/360
    
    dLat=np.reshape(LatLong1[:,0],(len(LatLong1[:,0]),1))-np.reshape(LatLong2[:,0],(1,len(LatLong2[:,0])))
    dLong=np.reshape(LatLong1[:,1],(len(LatLong1[:,1]),1))-np.reshape(LatLong2[:,1],(1,len(LatLong2[:,1])))
    
    P1=np.reshape(np.cos(LatLong1[:,0]),(LatLong1.shape[0],1))
    P2=np.reshape(np.cos(LatLong2[:,0]),(1,LatLong2.shape[0]))
    
    a=np.power(np.sin(dLat/2),2) + (P1*P2)*np.power(np.sin(dLong/2),2)
    c=2*np.arcsin(np.minimum(1,np.sqrt(a)))
    Distance=6367*c #[km]
    
    return Distance

def FixTimeResolutionCF(WindE_File, WaveE_File, OceanE_File):
    #The wind data starts at 1/1/2007 and ends at 12/31/2013
    WindEnergy=WindE_File['WindEnergy']
    WindLatLong_CF=WindE_File['LatLong']
    
    #The wave data starts at 1/1/2009 and ends at 12/31/2013
    WaveEnergy=WaveE_File['Energy_pu']
    WaveLatLong_CF=WaveE_File['LatLong']
    
    #The wave data starts at 1/1/2009 and ends at 11/30/2014
    OceanEnergy=OceanE_File['CurrentEnergy_pu']
    OceanLatLong_CF=OceanE_File['LatLong']
    
#-------------------------- WIND-----------------#
    SIdxWind=(dt.datetime(2009,1,1,0)-dt.datetime(2007,1,1,0)).days*24# Start index
    EIdxWind=(dt.datetime(2014,1,1,0)-dt.datetime(2007,1,1,0)).days*24# Last index
    
    Num3HourIntervals=int((EIdxWind-SIdxWind)/3)-1
    CFWind=np.zeros((WindEnergy.shape[0],Num3HourIntervals),dtype=float)
    
    for H3_Steps in range(Num3HourIntervals):
        
        StartIdx=SIdxWind + 3*(H3_Steps)
        
        CFWind[:,H3_Steps]=(WindEnergy[:,StartIdx])
        
#-------------------------- WAVE-----------------#
    SIdxWave=(dt.datetime(2009,1,1,0)-dt.datetime(2009,1,1,0)).days*8
    EIdxWave=((dt.datetime(2014,1,1,0)-dt.datetime(2009,1,1,0)).days)*8
    CFWave=WaveEnergy[:,SIdxWave:EIdxWave]

#-------------------------- OCEAN CURRENT-----------------#   
    SIdxCurrent=(dt.datetime(2009,1,1,0)-dt.datetime(2009,1,1,0)).days*8
    EIdxCurrent=((dt.datetime(2014,1,1,0)-dt.datetime(2009,1,1,0)).days)*8-1
    CFOcean=OceanEnergy[:,SIdxCurrent:EIdxCurrent]
    
    CFWind=np.average(CFWind, axis=1)
    CFWave=np.average(CFWave, axis=1)
    CFOcean=np.average(CFOcean, axis=1)

    return CFWind, CFWave, CFOcean, WindLatLong_CF, WaveLatLong_CF, OceanLatLong_CF

def GetSystemData(WindE_File, WaveE_File, OceanE_File):
    
    CFWind, CFWave, CFOcean, WindLatLong_CF, WaveLatLong_CF, OceanLatLong_CF=FixTimeResolutionCF(WindE_File, WaveE_File, OceanE_File)
    
    #[M$/Turbine][M$/Turbine]
    AnnualizedCostWind=WindE_File['AnnualizedCostWind']
    AnnualizedCostWave=WaveE_File['AnnualizedCostWave']
    AnnualizedCostOcean=OceanE_File['AnnualizedCostOcean']

    #[$/Mwh]
    LCOE_Wind=AnnualizedCostWind*10**6/  (CFWind*6*8766)
    LCOE_Wave=AnnualizedCostWave*10**6/  (CFWave*6*8766)
    LCOE_Ocean=AnnualizedCostOcean*10**6/(CFOcean*6*8766)
    
    return LCOE_Wind, LCOE_Wave, LCOE_Ocean, CFWind, CFWave, CFOcean, WindLatLong_CF, WaveLatLong_CF, OceanLatLong_CF
    

#########
def ComputeEquivalentCF (LatLong_Results, Y_NumTurbines, LatLong_Costs, CF):
    
    Distance=DistanceBetweenLatLong(LatLong_Costs, LatLong_Results)
    Idx_ClosestSiteCost=np.argmin(Distance,axis=0)
    Y_NumTurbines=np.reshape(Y_NumTurbines,(Y_NumTurbines.shape[0],1))
    
    CF=np.sum(CF[Idx_ClosestSiteCost,:]*Y_NumTurbines,axis=0)/sum(Y_NumTurbines)
    return CF


def TEMOA_INPUTS(LOCE_target, PathCapex_Opex_And_Simulations,PathPortOPT_Result, StartDate):
    WindData=np.load(PathCapex_Opex_And_Simulations+"WindEnergyNREL_100m_Haliade150_6MW.npz")
    WaveData=np.load(PathCapex_Opex_And_Simulations+"WaveEnergy_Pelamis_2009_2013.npz")
    OceanData=np.load(PathCapex_Opex_And_Simulations+"OceanCurrentEnergyRM4.npz")
    
    
    RatedPowerWind=WindData["RatedPower"]*10**-9 #[GW]
    RatedPowerWave=WaveData["RatedPower"]*10**-9
    RatedPowerOcean=OceanData["RatedPower"]*10**-9
    
    WindC_latlong=WindData["WindEnergy"]
    WaveC_latlong=WaveData["LatLong"]
    OceanC_latlong=OceanData["LatLong"]
    
    

    
LOCE_target="L" #H, M, L -> high, medium, low (LCOE)

PathCapex_Opex_And_Simulations="../"
PathPortOPT_Result="PortfolioOptimizationWindWaveOcean(50_0_0).npz"

WindData=np.load(PathCapex_Opex_And_Simulations+"WindEnergyNREL_100m_Haliade150_6MW.npz")
WaveData=np.load(PathCapex_Opex_And_Simulations+"WaveEnergy_Pelamis_2009_2013.npz")
OceanData=np.load(PathCapex_Opex_And_Simulations+"OceanCurrentEnergyRM4.npz")

    
#Capacity Factor
PortOPT_R=np.load(PathCapex_Opex_And_Simulations+PathPortOPT_Result,allow_pickle=True)
TotalNumTurbines=PortOPT_R["TotalNumTurbines"]
LCOEs=PortOPT_R["MINLPSolutionsLCOE"]
LCOEs=LCOEs[LCOEs!=None]
    
if LOCE_target=="L":
    CaseIDX=len(LCOEs)-1
    
if LOCE_target=="H":
    CaseIDX=0
    
if LOCE_target=="M":
    Avg=(LCOEs[0]+LCOEs[-1])/2
    CaseIDX=np.argmin(np.abs(LCOEs-Avg))

LCOE=LCOEs[CaseIDX]

WindR_latlong =PortOPT_R["MINLPSolutionsLatLongWind"]
WaveR_latlong =PortOPT_R["MINLPSolutionsLatLongWave"]
OceanR_latlong=PortOPT_R["MINLPSolutionsLatLongOcean"]


LCOE_Wind, LCOE_Wave, LCOE_Ocean, CFWind, CFWave, CFOcean, WindLatLong_CF, WaveLatLong_CF, OceanLatLong_CF=GetSystemData(WindData, WaveData, OceanData)



ShapeFileCoast="../GEO_data/ne_10m_coastline.shp"
ShapeFileStates="../GEO_data/ne_10m_admin_1_states_provinces_lines.shp"

min_longitude=-78.5
max_longitude=-74

min_latitude=33.5
max_latitude=37

xlim =[min_longitude,max_longitude]
ylim=[min_latitude, max_latitude]

df = gpd.read_file(ShapeFileCoast)
df1 = gpd.read_file(ShapeFileStates)

fig, ax = plt.subplots(figsize  = None)

df.plot(color='gray',linewidth=0.5,ax=ax)
df1.plot(color='gray',linewidth=0.5,ax=ax)

plt.scatter(WindLatLong_CF[:,1],WindLatLong_CF[:,0], c=LCOE_Wind, s=0.2, cmap='jet')
plt.scatter(WaveLatLong_CF[:,1],WaveLatLong_CF[:,0], c=LCOE_Wave, s=0.4, cmap='jet')
plt.scatter(OceanLatLong_CF[:,1],OceanLatLong_CF[:,0], c=LCOE_Ocean, s=0.2, cmap='jet')
plt.clim(np.floor(min(LCOE_Wind)),500)

clb = plt.colorbar()
clb.ax.set_title('LCOE\n $/MWh')

ax.set_xlim(xlim)
ax.set_ylim(ylim)
plt.xlabel("Longitude")
plt.ylabel("Latitude")

plt.title("Portfolio")

plt.savefig('LCOE.png',dpi=700)

############
ig, ax = plt.subplots(figsize  = None)

df.plot(color='gray',linewidth=0.5,ax=ax)
df1.plot(color='gray',linewidth=0.5,ax=ax)

Kite=np.load("../tempKite.npz")
LCOE_Kite=Kite["LCOE"]
LCOE_Lat=Kite["LAT_LCOE"]
LCOE_Long=Kite["LAT_LONG"]
# WindLatLong_CF=WindLatLong_CF[CFWind>0.25]
# WaveLatLong_CF=WaveLatLong_CF[CFWave>0.25]
# OceanLatLong_CF=OceanLatLong_CF[CFOcean>0.25]

# CFOcean=CFOcean[CFOcean>0.25]
# CFWave=CFWave[CFWave>0.25]
# CFOcean=CFOcean[CFOcean>0.25]
Filter=(LCOE_Kite<300) *(LCOE_Kite>10)


# plt.scatter(WaveLatLong_CF[:,1],WaveLatLong_CF[:,0], c=CFWave, s=1, cmap='jet')
# plt.scatter(OceanLatLong_CF[:,1],OceanLatLong_CF[:,0], c=CFOcean, s=0.2, cmap='jet')
# #plt.clim(np.floor(min(LCOE_Wind)),500)
X_LAT=np.zeros((int(len(LCOE_Lat)*len(LCOE_Long))),dtype=float)
Y_LONG=np.zeros((int(len(LCOE_Lat)*len(LCOE_Long))),dtype=float)
count=0
for I_lat in range(len(LCOE_Lat)):
    for I_long in range(len(LCOE_Long)):

      Y_LONG[count]=LCOE_Lat[I_lat]
      X_LAT[count]=LCOE_Long[I_long]
      count=count+1
          
          

# plt.scatter(X_LAT[np.reshape(Filter,-1)],Y_LONG[np.reshape(Filter,-1)], c=LCOE_Kite[Filter], s=0.2, cmap='jet')
# clb = plt.colorbar()
# clb.ax.set_title('Marginal LCOE[$/MWh]')


WindR_latlong_T=WindR_latlong[13]
P1=plt.scatter(WindR_latlong_T[:,1],WindR_latlong_T[:,0], s=15, marker="s",facecolor="none", color='k',)
    
if TotalNumTurbines[1]!=0:
    WaveR_latlong_T=WaveR_latlong[CaseIDX]
    P2=plt.scatter(WaveR_latlong_T[:,1],WaveR_latlong_T[:,0], s=15, marker="^",facecolor="none", color='k')
    

OceanR_latlong_T=OceanR_latlong[35]
P3=plt.scatter(OceanR_latlong_T[:,1]-0.25,OceanR_latlong_T[:,0]+0.5, s=15, marker="x", color='k')

# plt.legend((P1, P2, P3),
#            ('Wind', 'Wave', 'Ocean C.'),
#            scatterpoints=1,
#            loc='lower right',
#            ncol=1,
#            fontsize=8,
#            frameon=False)

plt.legend((P1, P3),
            ('Wind', 'Kite'),
            scatterpoints=1,
            loc='upper right',
            ncol=1,
            fontsize=8,
            frameon=False)

ax.set_xlim(xlim)
ax.set_ylim(ylim)
plt.xticks(np.arange(min(xlim), max(xlim)+.5, 0.5))
plt.yticks(np.arange(min(ylim), max(ylim)+.5, 0.5))
plt.xlabel("Longitude")
plt.ylabel("Latitude")

plt.savefig('CF.png',dpi=700)

