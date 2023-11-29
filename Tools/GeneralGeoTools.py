# env geo_env
# Tools to compute distance, and depth for ocean portfolios
from itertools import product
from numpy.random import randn
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import os
import csv
import xarray as xr
import pandas as pd
from tqdm import tqdm
import matplotlib.colors as mcolors
import sys
import geopandas as gpd

def GetDepthData(InputDataPath):   

    Depth_NETCDF = xr.open_dataset(InputDataPath+"./Depths.nc")
    
    Lat=Depth_NETCDF.lat.data
    Long=Depth_NETCDF.lon.data
    RefDepth=np.reshape(Depth_NETCDF.elevation.data,-1)
    RefDepth_LatLong = np.array(list(product(Lat, Long)))
    
    return RefDepth, RefDepth_LatLong    

def GetCoastLine_LatLong(InputDataPath):
    CoastLine=[]
    File_CoastLine=open(InputDataPath+'./Coastline/Coastline_NC.csv', "r")
    File_CoastLine_csv=csv.reader(File_CoastLine,delimiter=',')

    for EachLine in File_CoastLine_csv:

        if File_CoastLine_csv.line_num > 1:
            CoastLine.append([float(EachLine[1]), float(EachLine[0])] ) #LatLong

    CoastLine=np.array(CoastLine)

    return CoastLine

def GetDistanceToShore(InputDataPath, LatLong):
    #InputDataPath: Path with input data
    #LatLong: np.array with #Point,Lat[:,0]Long[:,1]     
    
    CoastLine=GetCoastLine_LatLong(InputDataPath)
    CoastLine=CoastLine*2*np.pi/360
    LatLong=LatLong*2*np.pi/360
    
    
    MaxSize=20# Maximum number of points to calculate distance at once (reduce memory)
    NIterations=int(np.ceil(LatLong.shape[0]/MaxSize))
    
    Distance=[]
    
    print("\n Calculating distance to shore for viable site locations")
    for i in tqdm(range(NIterations)):
        LatLong_i=LatLong[i*MaxSize:(i+1)*MaxSize,:]

        if LatLong_i.ndim==2:
            LatLong_i=np.reshape(LatLong_i,(np.shape(LatLong_i)[0],2))
        else:
            LatLong_i=np.reshape(LatLong_i,(1,2))   
            
        Lat=np.repeat(LatLong_i[:,0][:,np.newaxis], CoastLine[:,0].shape[0], axis=1)
        Long=np.repeat(LatLong_i[:,1][:,np.newaxis], CoastLine[:,1].shape[0], axis=1)
        CoastLat=np.repeat(np.reshape(CoastLine[:,0],(1,CoastLine.shape[0])), LatLong_i[:,0].shape[0], axis=0)
        CoastLong=np.repeat(np.reshape(CoastLine[:,1],(1,CoastLine.shape[0])), LatLong_i[:,1].shape[0], axis=0)
        
        dLat=Lat-CoastLat
        dLong=Long-CoastLong
            
            
        a=np.power(np.sin(dLat/2),2)+np.multiply(np.multiply(np.cos(CoastLat),np.cos(Lat)),np.power(np.sin(dLong/2),2))
        c=2*np.arcsin(np.minimum(1,np.sqrt(a)))
        d=6367*c
        
        Distance.append(np.min(d,axis=1)) #Minimum distance to shore km
    
    Distance=np.concatenate(Distance)
    
    #Distance: Minimum distance to shore for each point in LatLong [km]
    return Distance


def GetDepth(InputDataPath, LatLong):
    #Estimate depth in an array of latlong points
    
    #InputDataPath: Path with input data
    #LatLong: np.array with #Point,Lat[:,0]Long[:,1]     
    
    MaxSize=50# Maximum number of points to calculate distance at once (reduce memory)
    NIterations=int(np.ceil(LatLong.shape[0]/MaxSize))
    RefDepth, RefDepth_LatLong =GetDepthData(InputDataPath)
    
    #Filter for RefPoints in the range of LatLong (Reduce memory)
    IdxIn=np.where((RefDepth_LatLong[:,0]>=np.min(LatLong[:,0])-0.1) & (RefDepth_LatLong[:,0]<=np.max(LatLong[:,0])+0.1) & (RefDepth_LatLong[:,1]>=np.min(LatLong[:,1])-0.1) & (RefDepth_LatLong[:,1]<=np.max(LatLong[:,1])+0.1))[0]
    RefDepth=RefDepth[IdxIn]
    RefDepth_LatLong=RefDepth_LatLong[IdxIn,:]
    
    Depth=[]
    print("\n Calculating water depth at each point")
    for i in tqdm(range(NIterations)):
        LatLong_i=LatLong[i*MaxSize:(i+1)*MaxSize,:]

        if LatLong_i.ndim==2:
            LatLong_i=np.reshape(LatLong_i,(np.shape(LatLong_i)[0],2))
        else:
            LatLong_i=np.reshape(LatLong_i,(1,2))   
            
        Lat=np.repeat(LatLong_i[:,0][:,np.newaxis], RefDepth_LatLong[:,0].shape[0], axis=1)
        Long=np.repeat(LatLong_i[:,1][:,np.newaxis], RefDepth_LatLong[:,1].shape[0], axis=1)
        RefDephLat=np.repeat(np.reshape(RefDepth_LatLong[:,0],(1,RefDepth_LatLong.shape[0])), LatLong_i[:,0].shape[0], axis=0)
        RefDephLong=np.repeat(np.reshape(RefDepth_LatLong[:,1],(1,RefDepth_LatLong.shape[0])), LatLong_i[:,1].shape[0], axis=0)
        
        dLat=Lat-RefDephLat
        dLong=Long-RefDephLong
        
        #Complete version
        # a=np.power(np.sin(dLat/2),2)+np.multiply(np.multiply(np.cos(RefDephLat),np.cos(Lat)),np.power(np.sin(dLong/2),2))
        # c=2*np.arcsin(np.minimum(1,np.sqrt(a)))
        # d=6367*c
        
        #Simple Fast version
        d=np.sqrt(np.power(dLat,2)+np.power(dLong,2))
        
        Depth.append(RefDepth[np.argmin(d,axis=1)]) #Get depth for each point in LatLong [m]
    
    Depth=-np.concatenate(Depth)
    
    return Depth


def PlotGeneralGeoData(LatLong, Y_Variable, GeoDataPath, ColorBarTitle=None, Title=None, SavePath=None, s=6, LatMaxMin=(33.3, 37.2), LongMaxMin=(-78.7, -74.3)):

    ShapeFileCoast=GeoDataPath+"ne_10m_coastline.shp"
    ShapeFileStates=GeoDataPath+"ne_10m_admin_1_states_provinces_lines.shp"

    df = gpd.read_file(ShapeFileCoast)
    df1 = gpd.read_file(ShapeFileStates)


    fig, ax = plt.subplots(figsize  = None)
    df.plot(color='black',linewidth=1,ax=ax)
    df1.plot(color='black',linewidth=1,ax=ax)

    plt.scatter(LatLong[:,1],LatLong[:,0],c=Y_Variable, s=s, cmap="jet")

    clb = plt.colorbar()
    
    ax.set_xlim(LongMaxMin) 
    ax.set_ylim(LatMaxMin)
    plt.xlabel("Longitude", fontsize=12)
    plt.ylabel("Latitude", fontsize=12)
   
    if ColorBarTitle!=None:
        clb.ax.set_title(ColorBarTitle, fontsize=12) 

    if Title!=None:
        plt.title(Title, fontsize=12)

    if SavePath!=None:
        plt.savefig(SavePath, bbox_inches='tight', dpi=700)
        
        
    return plt.show()


def PlotGeneralGeoData_Class(LatLong, Y_Variable, GeoDataPath, Title=None, SavePath=None, s=6, LatMaxMin=(33.3, 37.2), LongMaxMin=(-78.7, -74.3)):

    ShapeFileCoast=GeoDataPath+"ne_10m_coastline.shp"
    ShapeFileStates=GeoDataPath+"ne_10m_admin_1_states_provinces_lines.shp"

    df = gpd.read_file(ShapeFileCoast)
    df1 = gpd.read_file(ShapeFileStates)


    fig, ax = plt.subplots(figsize  = None)
    df.plot(color='black',linewidth=1,ax=ax)
    df1.plot(color='black',linewidth=1,ax=ax)


    for label in np.unique(Y_Variable):
        IdIn=Y_Variable==label
        plt.scatter(LatLong[IdIn,1],LatLong[IdIn,0], s=s, label=label)


    ax.set_xlim(LongMaxMin) 
    ax.set_ylim(LatMaxMin)
    plt.xlabel("Longitude", fontsize=12)
    plt.ylabel("Latitude", fontsize=12)
    plt.legend(frameon=False)
   

    if Title!=None:
        plt.title(Title, fontsize=12)

    if SavePath!=None:
        plt.savefig(SavePath, bbox_inches='tight', dpi=700)
        
        
    return plt.show()