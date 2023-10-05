#Plot geo graphs for the wave CF
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt

ShapeFileCoast="./GEO_data/ne_10m_coastline.shp"
ShapeFileStates="./GEO_data/ne_10m_admin_1_states_provinces_lines.shp"

HourlyData=np.load('WaveEnergy_Pelamis_2009_2013.npz')

HourWaveEnergy=HourlyData["Energy_pu"]
DepthSites=HourlyData["DepthSites"]
ShoreDistance=HourlyData["ShoreDistance"]
AvgWaveEnergy=np.average(HourWaveEnergy,axis=1)

LatLong=HourlyData["LatLong"]

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

plt.title("Wave Energy Sites 2009-2013\n(WWIII)")

plt.savefig('CF_WaveNC_2009_2013.png',dpi=700)

Record=np.load("WaveHsTp_WWIII.npz")    
LatLong=Record["LatLong"]
ShoreDistance=Record["ShoreDistance"]
DepthSites=Record["Depth"]

#### Depth
fig, ax = plt.subplots(figsize  = None)

df.plot(color='black',linewidth=1,ax=ax)
df1.plot(color='black',linewidth=1,ax=ax)

plt.scatter(LatLong[:,1],LatLong[:,0],c=DepthSites, s=8, cmap='jet')

clb = plt.colorbar()
clb.ax.set_title('[m]')

ax.set_xlim(xlim)
ax.set_ylim(ylim)
plt.xlabel("Longitude")
plt.ylabel("Latitude")


plt.savefig('Depth_WaveNC_2009_2013.png',dpi=700)

#### Distance
fig, ax = plt.subplots(figsize  = None)

df.plot(color='black',linewidth=1,ax=ax)
df1.plot(color='black',linewidth=1,ax=ax)

plt.scatter(LatLong[:,1],LatLong[:,0],c=ShoreDistance, s=8, cmap='jet')

clb = plt.colorbar()
clb.ax.set_title('[km]')

ax.set_xlim(xlim)
ax.set_ylim(ylim)
plt.xlabel("Longitude")
plt.ylabel("Latitude")


plt.savefig('Distance_WaveNC_2009_2013.png',dpi=700)