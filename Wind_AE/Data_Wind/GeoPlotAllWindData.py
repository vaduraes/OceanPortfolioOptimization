#Plot geo graph for the wind speed- do not limit for the feasible site locations (run all NC region)

import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt

ShapeFileCoast="./GEO_data/ne_10m_coastline.shp"
ShapeFileStates="./GEO_data/ne_10m_admin_1_states_provinces_lines.shp"

Data=np.load('windspeed_100m.npz')

windspeed=Data["windspeed"]
windspeed=windspeed.T

windspeed=windspeed[:,17544:61320] #Indexes for 2009-2013
AvgWindSpeed=np.average(windspeed,axis=1)

LatLong=Data["LatLong"]

min_longitude=-78.7
max_longitude=-74.5

min_latitude=34
max_latitude=36.7

xlim =[min_longitude,max_longitude]
ylim=[min_latitude, max_latitude]

df = gpd.read_file(ShapeFileCoast)
df1 = gpd.read_file(ShapeFileStates)

fig, ax = plt.subplots(figsize  = None)

df.plot(color='black',linewidth=1,ax=ax)
df1.plot(color='black',linewidth=1,ax=ax)

plt.scatter(LatLong[:,1],LatLong[:,0],c=AvgWindSpeed, s=0.1, cmap='jet')

clb = plt.colorbar()
clb.ax.set_title('m/s')

ax.set_xlim(xlim)
ax.set_ylim(ylim)
plt.xlabel("Longitude")
plt.ylabel("Latitude")

plt.title("Wind Speed 2009-2013\n(NREL)")

plt.savefig('WindSpeed_NC_2009_2013.png',dpi=700)