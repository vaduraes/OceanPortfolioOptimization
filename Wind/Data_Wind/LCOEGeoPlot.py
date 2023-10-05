#Geo plot of LCOE for wind
import numpy as np
import datetime as dt
import geopandas as gpd
import csv
import matplotlib.colors as clrs
import matplotlib.pyplot as plt

ShapeFileCoast="./GEO_data/ne_10m_coastline.shp"
ShapeFileStates="./GEO_data/ne_10m_admin_1_states_provinces_lines.shp"

HourlyData=np.load('WindEnergyNREL_100m_Haliade150_6MW.npz', allow_pickle=True)

#Indexes from 2009-2013
FromDate=dt.datetime(2007,1,1,0)
FirstIDX=(dt.datetime(2009,1,1,0)-FromDate).days*24
LastIDX=(dt.datetime(2014,1,1,0)-FromDate).days*24

HourlyCurrentEnergy=HourlyData["WindEnergy"][:,FirstIDX: LastIDX]
AvgCurrentEnergy=np.average(HourlyCurrentEnergy,axis=1)

LatLong=HourlyData["LatLong"]
DepthSites=HourlyData["Depth"]
ShoreDistance=HourlyData["DistanceToShore"]


AnnualizedCostWind=HourlyData["AnnualizedCostWind"]
RatedPower=HourlyData["RatedPower"]/10**6

LCOE=[]
for i in range(len(AvgCurrentEnergy)):
    LCOE_i=AnnualizedCostWind[i]/(AvgCurrentEnergy[i]*RatedPower*8766)*10**6
    LCOE.append(LCOE_i)

LCOE=np.array(LCOE)


min_longitude=-78.7
max_longitude=-74.5

min_latitude=33.7
max_latitude=36.7

xlim =[min_longitude,max_longitude]
ylim=[min_latitude, max_latitude]

df = gpd.read_file(ShapeFileCoast)
df1 = gpd.read_file(ShapeFileStates)

fig, ax = plt.subplots(figsize  = None)

df.plot(color='black',linewidth=1,ax=ax)
df1.plot(color='black',linewidth=1,ax=ax)

cmap = clrs.LinearSegmentedColormap.from_list("", ["red","yellow","green","blue"])

plt.scatter(LatLong[:,1],LatLong[:,0], c=LCOE, s=0.2, cmap='jet')
plt.clim(np.floor(min(LCOE)),140)

clb = plt.colorbar()
clb.ax.set_title('LCOE\n $/MWh')


ax.set_xlim(xlim)
ax.set_ylim(ylim)
plt.xlabel("Longitude")
plt.ylabel("Latitude")

plt.title("Wind Energy LCOE")

plt.savefig('LCOE_WindNC_2009_2013.png',dpi=700)


#### Depth
fig, ax = plt.subplots(figsize  = None)

df.plot(color='black',linewidth=1,ax=ax)
df1.plot(color='black',linewidth=1,ax=ax)

plt.scatter(LatLong[:,1],LatLong[:,0],c=DepthSites, s=0.2, cmap='jet',vmax=45)

clb = plt.colorbar()
clb.ax.set_title('[m]')


ax.set_xlim(xlim)
ax.set_ylim(ylim)
plt.xlabel("Longitude")
plt.ylabel("Latitude")


plt.savefig('Depth_WindNC_2009_2013.png',dpi=700)

#### Distance
fig, ax = plt.subplots(figsize  = None)

df.plot(color='black',linewidth=1,ax=ax)
df1.plot(color='black',linewidth=1,ax=ax)

plt.scatter(LatLong[:,1],LatLong[:,0],c=ShoreDistance, s=0.2, cmap='jet')

clb = plt.colorbar()
clb.ax.set_title('[km]')

ax.set_xlim(xlim)
ax.set_ylim(ylim)
plt.xlabel("Longitude")
plt.ylabel("Latitude")


plt.savefig('Distance_WindNC_2009_2013.png',dpi=700)
