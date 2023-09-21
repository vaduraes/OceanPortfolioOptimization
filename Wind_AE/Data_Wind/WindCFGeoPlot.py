#Geo plot of CF for wind energy
import numpy as np
import pandas as pd
import geoplot as gplt
import geopandas as gpd
import matplotlib.colors as clrs
import matplotlib.pyplot as plt


def GeoPlot(Variable, LatLong, SaveName, CF=False, WS=False, NameCB="None"):
    ShapeFileCoast="./GEO_data/ne_10m_coastline.shp"
    ShapeFileStates="./GEO_data/ne_10m_admin_1_states_provinces_lines.shp"

    min_longitude=-78.7
    max_longitude=-74.5

    min_latitude=33.0
    max_latitude=37.0

    xlim =[min_longitude,max_longitude]
    ylim=[min_latitude, max_latitude]

    df = gpd.read_file(ShapeFileCoast)
    df1 = gpd.read_file(ShapeFileStates)

    fig, ax = plt.subplots(figsize  = None)

    df.plot(color='black',linewidth=1,ax=ax)
    df1.plot(color='black',linewidth=1,ax=ax)

    
    if CF==True:
        plt.title("CF Wind Energy Sites\n2007-2014")
    if WS==True:
        plt.title("Wind Speed\n2007-2014")


    plt.scatter(LatLong[:,1],LatLong[:,0],c=Variable, s=0.2, cmap='jet')
    clb = plt.colorbar()

    if CF==True:
        clb.ax.set_title('CF')
        plt.title("CF Wind Energy Sites\n2007-2014")
    if WS==True:
        clb.ax.set_title('[m/s]')
        plt.title("Wind Speed\n2007-2014")
    if NameCB!="None":
        clb.ax.set_title(NameCB)
        

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")

    plt.savefig("./Figures/"+SaveName+'_2007_2014.png',dpi=700, bbox_inches = 'tight')


