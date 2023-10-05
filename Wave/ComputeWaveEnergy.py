#Compute wave energy

import numpy as np
from ReadTurbineData import ReadTurbineData
import geopandas as gpd
import matplotlib.pyplot as plt

#Get turbine data
TurbineFile='RM3.csv'
Turbine=ReadTurbineData(TurbineFile)

RatedPower=Turbine["RatedPower"]
E_Mec2El=Turbine["E_Mec2El"]
E_Av=Turbine["E_Av"]
E_Tr=Turbine["E_Tr"]
Te_Bins=Turbine["Te_Bins"]
Hs_Bins=Turbine["Hs_Bins"]
MP_Matrix=Turbine["MP_Matrix"]

#Get wave data
WaveData='WaveHsTp_WWIII.npz'
WaveFile=np.load(WaveData)
HsNC=WaveFile["HsNC"]
TpNC=WaveFile["TpNC"]

LatLong=WaveFile["LatLong"]

HsShape=HsNC.shape

#Create an one dimensional vector with the Hs data. This procedure facilitates significantly the future steps
D1_HsNC=np.reshape(HsNC,(np.size(HsNC),1))

#Index of the Hs set value closest to the observed Hs value
IdxHs=np.reshape(np.argmin(np.abs(D1_HsNC-Hs_Bins.T),axis=1), HsShape)

#Create an one dimensional vector with the Tp data.
D1_TpNC=np.reshape(TpNC,(np.size(TpNC),1))

#Index of the Tp set value closest to the observed Tp value
IdxTp=np.reshape(np.argmin(np.abs(D1_TpNC-Te_Bins.T),axis=1),HsShape)

EnergyProduction=np.minimum(MP_Matrix[IdxHs,IdxTp]*E_Mec2El,RatedPower)

WakeEffect=0.95
EnergyPu=EnergyProduction/RatedPower*WakeEffect

EnergyPu=EnergyPu*E_Av*E_Tr

    
    #### Plot
ShapeFileCoast="./GEO_data/ne_10m_coastline.shp"
ShapeFileStates="./GEO_data/ne_10m_admin_1_states_provinces_lines.shp"

DepthSites=WaveFile["Depth"]
ShoreDistance=WaveFile["ShoreDistance"]
AvgWaveEnergy=np.average(EnergyPu,axis=1)

min_longitude=-78.7
max_longitude=-74.5

min_latitude=33.7
max_latitude=36.7

xlim =[min_longitude,max_longitude]
ylim=[min_latitude, max_latitude]

df = gpd.read_file(ShapeFileCoast)
df1 = gpd.read_file(ShapeFileStates)

#### CF
fig, ax = plt.subplots(figsize  = None)

df.plot(color='black',linewidth=1,ax=ax)
df1.plot(color='black',linewidth=1,ax=ax)

plt.scatter(LatLong[:,1],LatLong[:,0],c=AvgWaveEnergy, s=8, cmap='jet')

clb = plt.colorbar()
clb.ax.set_title('CF')

ax.set_xlim(xlim)
ax.set_ylim(ylim)
plt.xlabel("Longitude")
plt.ylabel("Latitude")

plt.savefig('CF'+TurbineFile+".png",dpi=700)
