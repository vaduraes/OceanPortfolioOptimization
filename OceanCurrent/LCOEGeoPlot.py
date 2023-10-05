#Geo plot of the LCOE for ocean current energy
import numpy as np
import geopandas as gpd
import datetime as dt
import csv
import matplotlib.colors as clrs
import matplotlib.pyplot as plt

#Distance between two lat long points
def DistanceToShore (CoastLine, LatLong1): #Compute distance to shore in km of a lat long point
   
    CoastLine=CoastLine*2*np.pi/360
    LatLong1=LatLong1*2*np.pi/360
    
    LatLong1=np.reshape(LatLong1,(1,2))
    dLat=LatLong1[:,0]-CoastLine[:,0]
    dLong=LatLong1[:,1]-CoastLine[:,1]
    
    a=np.power(np.sin(dLat/2),2)+np.cos(CoastLine[:,0])*np.cos(LatLong1[:,0])*np.power(np.sin(dLong/2),2)
    c=2*np.arcsin(np.minimum(1,np.sqrt(a)))
    d=6367*c
    
    Distance=np.min(d) #Minimum distance to shore km
    
    return Distance

def CurveDistanceFromShore(CoastLine, TargetDistance):
    LatLongForDistance=[]

    for Lat in np.arange(33,36.2,0.005): 
        MinD=1000
        StoreLong=0
        
        for Long in np.arange(-74,-79,-0.005):
            D=DistanceToShore(CoastLine, np.array([Lat, Long]))      
            
            if abs(D-TargetDistance)<MinD:
                MinD=abs(D-TargetDistance)
                StoreLong=Long
                
            if abs(D-TargetDistance)<0.6 or D<0.6:
                break
            
        LatLongForDistance.append([Lat, StoreLong])

    return np.array(LatLongForDistance)

#Get coastline data
CoastLine=[]

File_CoastLine=open('./GEO_data/Coastline_NC.csv', "r")
File_CoastLine_csv=csv.reader(File_CoastLine,delimiter=',')

for EachLine in File_CoastLine_csv:

    if File_CoastLine_csv.line_num > 1:
        CoastLine.append([float(EachLine[1]), float(EachLine[0])] ) #LatLong


CoastLine=np.array(CoastLine)

ShapeFileCoast="./GEO_data/ne_10m_coastline.shp"
ShapeFileStates="./GEO_data/ne_10m_admin_1_states_provinces_lines.shp"

HourlyData=np.load('OceanCurrentEnergyRM4.npz', allow_pickle=True)

LatLong=HourlyData["LatLong"]
HourlyCurrentEnergy=HourlyData["CurrentEnergy_pu"]
AnnualizedCostOcean=HourlyData["AnnualizedCostOcean"]
OceanDateTime=HourlyData['OceanDateTime']

#Last index for the last data of 2013
Idx_2013=np.where(OceanDateTime==dt.datetime(2013,12,31,21))[0][0]

AvgCurrentEnergy=np.average(HourlyCurrentEnergy[:,0:Idx_2013],axis=1)

LCOE=[]
for i in range(len(AvgCurrentEnergy)):
    if AvgCurrentEnergy[i]!=0:
        LCOE_i=AnnualizedCostOcean[i]/(AvgCurrentEnergy[i]*4*8766)*10**6#[$]
        LCOE.append(LCOE_i)
        
    else:
        LCOE.append(1000000)#No ocean current in the region

LCOE=np.array(LCOE)

LatLong66km=CurveDistanceFromShore(CoastLine,  66)
LatLong100km=CurveDistanceFromShore(CoastLine, 100)
LatLong150km=CurveDistanceFromShore(CoastLine, 150)

min_longitude=-77
max_longitude=-74

min_latitude=33
max_latitude=36.2

xlim =[min_longitude,max_longitude]
ylim=[min_latitude, max_latitude]

df = gpd.read_file(ShapeFileCoast)
df1 = gpd.read_file(ShapeFileStates)

fig, ax = plt.subplots(figsize  = None)

df.plot(color='black',linewidth=1,ax=ax)
df1.plot(color='black',linewidth=1,ax=ax)

cmap = clrs.LinearSegmentedColormap.from_list("", ["red","yellow","green","blue"])

plt.plot(LatLong66km[:,1],LatLong66km[:,0], linestyle='--', color='black', linewidth=0.8)
plt.plot(LatLong100km[:,1],LatLong100km[:,0], linestyle='--', color='black', linewidth=0.8)
plt.plot(LatLong150km[:,1],LatLong150km[:,0], linestyle='--', color='black', linewidth=0.8)

plt.scatter(LatLong[:,1],LatLong[:,0], c=LCOE, s=0.2, cmap='jet')
plt.clim(np.floor(min(LCOE)),500)

clb = plt.colorbar()
clb.ax.set_title('LCOE\n $/MWh')

plt.text(-76.95,34.05,"66km",fontsize=8)

plt.text(-75,34.5,"100km", fontsize=8,  bbox={'facecolor':'white', 'edgecolor':'none', 'pad':1})

plt.text(-75,34,"150km", fontsize=8,  bbox={'facecolor':'white', 'edgecolor':'none', 'pad':1})

ax.set_xlim(xlim)
ax.set_ylim(ylim)
plt.xlabel("Longitude")
plt.ylabel("Latitude")

plt.title("Ocean Current LCOE")

plt.savefig('LCOE_CurrentNC_2009_2013.png',dpi=700)

