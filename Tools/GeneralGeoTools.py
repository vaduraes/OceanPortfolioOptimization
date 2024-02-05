# env geo_env
# Tools to compute distance, and depth for ocean portfolios
from itertools import product
from numpy.random import randn
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

from tqdm import tqdm
import os
import csv
import xarray as xr
import pandas as pd
from tqdm import tqdm
import matplotlib.colors as mcolors
import sys
import geopandas as gpd
import re
from datetime import datetime, timedelta



def GetTimeList(start_date, end_date, TimeDeltaHours=1):
    
    # start_date = datetime(2007, 1, 1)
    # end_date = datetime(2013, 12, 31, 23)

    TimeList = []

    current_date = start_date
    while current_date <= end_date:
        TimeList.append(current_date)
        current_date += timedelta(hours=TimeDeltaHours)
        
    return TimeList

def GetDepthData(InputDataPath):   

    Depth_NETCDF = xr.open_dataset(InputDataPath+"./Depths.nc")
    
    Lat=Depth_NETCDF.lat.data
    Long=Depth_NETCDF.lon.data
    RefDepth=np.reshape(Depth_NETCDF.elevation.data,-1)
    RefDepth_LatLong = np.array(list(product(Lat, Long)))
    
    return RefDepth, RefDepth_LatLong    

def GetCoastLine_LatLong(InputDataPath):
    CoastLine=[]
    File_CoastLine=open(InputDataPath+'./Coastline/Coastline_Full_NC.csv', "r")
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

    MaxSize=10# Maximum number of points to calculate distance at once (reduce memory)
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
                
        IdxIn=np.where((RefDepth_LatLong[:,0]>=np.min(LatLong_i[:,0])-0.1) * (RefDepth_LatLong[:,0]<=np.max(LatLong_i[:,0])+0.1) * (RefDepth_LatLong[:,1]>=np.min(LatLong_i[:,1])-0.1) * (RefDepth_LatLong[:,1]<=np.max(LatLong_i[:,1])+0.1))[0]
        RefDepth_i=RefDepth[IdxIn]
        RefDepth_LatLong_i=RefDepth_LatLong[IdxIn,:]  
        
        
        Lat=np.repeat(LatLong_i[:,0][:,np.newaxis], RefDepth_LatLong_i[:,0].shape[0], axis=1)
        Long=np.repeat(LatLong_i[:,1][:,np.newaxis], RefDepth_LatLong_i[:,1].shape[0], axis=1)
        RefDephLat=np.repeat(np.reshape(RefDepth_LatLong_i[:,0],(1,RefDepth_LatLong_i.shape[0])), LatLong_i[:,0].shape[0], axis=0)
        RefDephLong=np.repeat(np.reshape(RefDepth_LatLong_i[:,1],(1,RefDepth_LatLong_i.shape[0])), LatLong_i[:,1].shape[0], axis=0)
        
        dLat=Lat-RefDephLat
        dLong=Long-RefDephLong
        
        #Complete version
        # a=np.power(np.sin(dLat/2),2)+np.multiply(np.multiply(np.cos(RefDephLat),np.cos(Lat)),np.power(np.sin(dLong/2),2))
        # c=2*np.arcsin(np.minimum(1,np.sqrt(a)))
        # d=6367*c
        
        #Simple Fast version
        d=np.sqrt(np.power(dLat,2)+np.power(dLong,2))
        
        Depth.append(RefDepth_i[np.argmin(d,axis=1)]) #Get depth for each point in LatLong [m]

    Depth=-np.concatenate(Depth)
    
    return Depth

def PlotGeneralGeoData(LatLong, Y_Variable, GeoDataPath="./InputData/CoastLine/", ColorBarTitle=None, Title=None, SavePath=None, s=6, LatMaxMin=(33.3, 37.2), LongMaxMin=(-78.7, -74.3)):

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



def custom_sort(value):
    # Split the string into letters and numbers
    letters, numbers = re.match(r'([a-zA-Z]+)([0-9]+)', value).groups()
    return (letters, int(numbers))  # Convert numbers to integers for sorting



def PlotGeneralGeoData_Class(LatLong, Y_Variable, GeoDataPath, Title=None, SavePath=None, s=6, LatMaxMin=(33.3, 37.2), LongMaxMin=(-78.7, -74.3)):

    ShapeFileCoast=GeoDataPath+"ne_10m_coastline.shp"
    ShapeFileStates=GeoDataPath+"ne_10m_admin_1_states_provinces_lines.shp"

    df = gpd.read_file(ShapeFileCoast)
    df1 = gpd.read_file(ShapeFileStates)


    fig, ax = plt.subplots(figsize  = None)
    df.plot(color='black',linewidth=1,ax=ax)
    df1.plot(color='black',linewidth=1,ax=ax)

    # Sort the legends using the custom sorting function
    if any(char.isdigit() for char in np.unique(Y_Variable)[0]):#there is a number in the legend
        sorted_legends = sorted(np.unique(Y_Variable), key=custom_sort)

        for label in sorted_legends:
            IdIn=Y_Variable==label
            plt.scatter(LatLong[IdIn,1],LatLong[IdIn,0], s=s, label=label)
    else:
        for label in np.unique(Y_Variable):
            IdIn=Y_Variable==label
            plt.scatter(LatLong[IdIn,1],LatLong[IdIn,0], s=s, label=label)  


    ax.set_xlim(LongMaxMin) 
    ax.set_ylim(LatMaxMin)
    plt.xlabel("Longitude", fontsize=12)
    plt.ylabel("Latitude", fontsize=12)
    
    #legend outsize figure
    plt.legend(frameon=False, bbox_to_anchor=(1.04, 0.5), loc="center left")
   

    if Title!=None:
        plt.title(Title, fontsize=12)

    if SavePath!=None:
        plt.savefig(SavePath, bbox_inches='tight', dpi=700)
        
        
    return plt.show()


#Change space and time resolution of a energy dataset
def ChangeTimeSpaceResolution(ReferenceDataPath, CurrentTimeResolution, NewTimeResolution, StepsPerDegree, StartDateTime, EndDateTime, NewSavePath=None):
    Data=np.load(ReferenceDataPath,allow_pickle=True)

    Energy_Pu=Data["Energy_pu"]
    RawResource=Data["RawResource"]
    RatedPower=Data["RatedPower"]
    TimeList=Data["TimeList"]
    LatLong=Data["LatLong"]
    Depth=Data["Depth"]
    DistanceShore=Data["DistanceShore"]
    CAPEX_site=Data["CAPEX_site"]
    OPEX_site=Data["OPEX_site"]
    AnnualizedCost=Data["AnnualizedCost"]

    #######--------- Change Time Resolution and start end dates
    NewTimeList=GetTimeList(StartDateTime, EndDateTime, TimeDeltaHours=NewTimeResolution)
    NewTimeList=np.array(NewTimeList) #Reporting the start of the time period only
    TimeList=np.array(TimeList) 

    NewEnergy_Pu_time=np.zeros((len(NewTimeList),Energy_Pu.shape[1]))
    NewRawResource_time=np.zeros((len(NewTimeList),Energy_Pu.shape[1]))

    for i in range(len(NewTimeList)):
        IdxIn=(TimeList>=NewTimeList[i])*(TimeList<NewTimeList[i]+timedelta(hours=NewTimeResolution))
        NewEnergy_Pu_time[i,:]=np.average(Energy_Pu[IdxIn,:],axis=0)
        NewRawResource_time[i,:]=np.average(RawResource[IdxIn,:],axis=0)

    #######--------- Change Spatial Resolution
    #Create Grid of Latitudes and Longitudes   
    ResolutionDegrees=1/StepsPerDegree
    
    lat_min=np.min(LatLong[:,0])-1/StepsPerDegree
    lat_max=np.max(LatLong[:,0])+1/StepsPerDegree
    lon_min=np.min(LatLong[:,1])-1/StepsPerDegree
    lon_max=np.max(LatLong[:,1])+1/StepsPerDegree

    xlim =[lon_min,lon_max]
    ylim =[lat_min, lat_max]

    lat_range = np.linspace(lat_min, lat_max,int((lat_max-lat_min)*StepsPerDegree))
    lon_range = np.linspace(lon_min, lon_max, int((lon_max-lon_min)*StepsPerDegree))
            
    NewEnergy_Pu_TimeSpace=[]
    NewRawResource_TimeSpace=[]
    NewLatLong=[]
    NewDepth=[]
    NewDistanceShore=[]
    NewCAPEX_site=[]
    NewOPEX_site=[]
    NewAnnualizedCost=[]
    NumberOfCellsPerSite=[]

    for lat in lat_range:
        for long in lon_range:
            IdxIn=(LatLong[:,0]>=lat)*(LatLong[:,0]<lat+1/StepsPerDegree)*(LatLong[:,1]>=long)*(LatLong[:,1]<long+1/StepsPerDegree)
            
            if np.sum(IdxIn)>0:
                NewEnergy_Pu_TimeSpace.append(np.average(NewEnergy_Pu_time[:,IdxIn],axis=1))
                NewRawResource_TimeSpace.append(np.average(NewRawResource_time[:,IdxIn],axis=1))
                #NewLatLong.append([np.average(LatLong[IdxIn,0]),np.average(LatLong[IdxIn,1])])
                NewLatLong.append([lat+1/StepsPerDegree*0.5,long + 1/StepsPerDegree*0.5])#Center of the cell
                NewDepth.append(np.average(Depth[IdxIn]))
                NewDistanceShore.append(np.average(DistanceShore[IdxIn]))
                NewCAPEX_site.append(np.average(CAPEX_site[IdxIn]))
                NewOPEX_site.append(np.average(OPEX_site[IdxIn]))
                NewAnnualizedCost.append(np.average(AnnualizedCost[IdxIn]))
                NumberOfCellsPerSite.append(np.sum(IdxIn))
            
    NewEnergy_Pu_TimeSpace=np.array(NewEnergy_Pu_TimeSpace).T
    NewRawResource_TimeSpace=np.array(NewRawResource_TimeSpace).T
    NewLatLong=np.array(NewLatLong)
    NewDepth=np.array(NewDepth)
    NewDistanceShore=np.array(NewDistanceShore)
    NewCAPEX_site=np.array(NewCAPEX_site)
    NewOPEX_site=np.array(NewOPEX_site)
    NewAnnualizedCost=np.array(NewAnnualizedCost)
    NumberOfCellsPerSite=np.array(NumberOfCellsPerSite)

    if NewSavePath!=None:
        np.savez(NewSavePath, Energy_pu=NewEnergy_Pu_TimeSpace, RawResource=NewRawResource_TimeSpace, TimeList=NewTimeList, LatLong=NewLatLong, Depth=NewDepth,\
            DistanceShore=NewDistanceShore, CAPEX_site=NewCAPEX_site, OPEX_site=NewOPEX_site, AnnualizedCost=NewAnnualizedCost, NumberOfCellsPerSite=NumberOfCellsPerSite,RatedPower=RatedPower, 
            ResolutionDegrees=ResolutionDegrees, ResolutionKm=-1)

    return NewEnergy_Pu_TimeSpace, NewTimeList, NewRawResource_TimeSpace, NewLatLong, NewDepth, NewDistanceShore, NewCAPEX_site, NewOPEX_site, NewAnnualizedCost, NumberOfCellsPerSite, RatedPower
  

#Compute minimum distance from set of points
#Index from SetPoints1 that is closest to SetPoints2 element by element
#Usually SetPoints1 is the set of points in the grid and SetPoints2 is the set of points in the portfolio
def MinDistanceSetPoints (SetPoints1, SetPoints2): 
    SetPoints1=SetPoints1*2*np.pi/360
    SetPoints2=SetPoints2*2*np.pi/360

    DMin=[]
    IdxMin=[]

    for i in range(SetPoints2.shape[0]):
        LatLong1=SetPoints2[i,:]
        LatLong1=np.reshape(LatLong1,(1,2))
        dLat=LatLong1[:,0]-SetPoints1[:,0]
        dLong=LatLong1[:,1]-SetPoints1[:,1]

        a=np.power(np.sin(dLat/2),2)+np.cos(SetPoints1[:,0])*np.cos(LatLong1[:,0])*np.power(np.sin(dLong/2),2)
        c=2*np.arcsin(np.minimum(1,np.sqrt(a)))
        d=6367*c

        Distance=np.min(d) #Minimum distance km
        DMin.append(Distance)
        Idx=np.argmin(d)
        IdxMin.append(Idx)
        
    IdxMin=np.array(IdxMin)
    DMin=np.array(DMin)
    
    
    return IdxMin, DMin

def PlotsWithBOEM(GeoDataPath, PathBOEMData, BOEM_ShpDir, LatLong1, Y_Variable1, LatLong2=None, Y_Variable2=None, ColorBarTitle1=None, ColorBarTitle2=None, Title=None, SavePath=None, s=6, LatMaxMin=(32,40), LongMaxMin=(-83,-72)):

    ShapeFileCoast=GeoDataPath+"ne_10m_coastline.shp"
    ShapeFileStates=GeoDataPath+"ne_10m_admin_1_states_provinces_lines.shp"

        
    df = gpd.read_file(ShapeFileCoast)
    df1 = gpd.read_file(ShapeFileStates)

    fig, ax = plt.subplots()
    fig.set_size_inches(18.5, 10.5)

    df.plot(color='black',linewidth=1,ax=ax)
    df1.plot(color='black',linewidth=1,ax=ax)


    kings_county_map = gpd.read_file(PathBOEMData+BOEM_ShpDir["Wind Lease Areas"])
    kings_county_map.plot(ax=ax, color='green', edgecolor='none', alpha=0.3,legend=True)

    kings_county_map = gpd.read_file(PathBOEMData+BOEM_ShpDir["Wind Planning Areas"])
    kings_county_map.plot(ax=ax, color='red', edgecolor='none', alpha=0.3)

    ax.arrow( -73.7, 37.3, 0.5,-0.5, fc="k", ec="k",head_width=0.05, head_length=0.1)
    plt.text(-73.7+0.5-0.4, 37.3-0.5-0.25, 'Wind Planning Area', fontsize = 10)

    ax.arrow( -73.5, 39.4, 0.5,-0.5, fc="k", ec="k",head_width=0.05, head_length=0.1)
    plt.text(-73.5+0.5-0.6, 39.4-0.5-0.25, 'Wind Lease Area', fontsize = 10)


    plt.scatter(LatLong1[:,1], LatLong1[:,0], c=Y_Variable1, s=s,cmap="jet", alpha=0.3)
    clb=plt.colorbar()
    clb.ax.set_title(ColorBarTitle1)
    
    try:
        plt.scatter(LatLong2[:,1], LatLong2[:,0], c=Y_Variable2, s=s,cmap="jet", alpha=0.3)
        clb=plt.colorbar()
        clb.ax.set_title(ColorBarTitle2)
    except:
        pass

    ax.set_xlim(LongMaxMin) 
    ax.set_ylim(LatMaxMin)
    plt.xlabel("Longitude", fontsize=12)
    plt.ylabel("Latitude", fontsize=12)
    plt.title(Title)
    
    if SavePath!=None:
        plt.savefig(SavePath, bbox_inches='tight', dpi=700)
        
        
def PlotEfficientFrontier(SolutionPaths, Legend, Title, linestyle=None,ColorList=None, SavePath=None):
    
    if linestyle==None:
        linestyle=['-','-','-','-', '--','--','--','--',':',':',':',':','-.','-.','-.','-.']
        ColorList=['k','b','g','r', 'k','b','g','r', 'k','b','g','r', 'k','b','g','r'] 
              
    for i, SolutionPath in enumerate(SolutionPaths):
        Solution=np.load(SolutionPath,allow_pickle=True)
        LCOE=Solution["Save_LCOE_Achieved"]
        MWAvg=Solution["SaveTotalMWAvg"]
        
        LCOE=LCOE[MWAvg!=None]
        
        #Temporary solution for some simulations that I forgot to save the MWAvg and saved the MWH
        if MWAvg[MWAvg!=None][0]>5000: #5GW
            MWAvg=MWAvg[MWAvg!=None]/(24*365.25)
        else:
            MWAvg=MWAvg[MWAvg!=None]
        
        plt.plot(MWAvg,LCOE, label=Legend[i],linestyle=linestyle[i],color=ColorList[i])

    plt.xlabel("MW Avg - Delivered to Shore")
    plt.ylabel("LCOE [$/MWh]")
    plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0., frameon=False)
    
    if SavePath!=None:
        plt.savefig(SavePath, bbox_inches='tight', dpi=700)

def PlotPowerTechDistribution(SolutionPath, SavePath=None):

    Data=np.load(SolutionPath,allow_pickle=True)
    Save_TotalMWAvgWind=Data["Save_TotalMWAvgWind"]
    Save_TotalMWAvgWave=Data["Save_TotalMWAvgWave"]
    Save_TotalMWAvgKite=Data["Save_TotalMWAvgKite"]
    Save_totalMWAvgCurtailment=Data["Save_totalMWAvgCurtailment"]
    SaveTotalMWAvg=Data["SaveTotalMWAvg"]
    Save_LCOE_Achieved=Data["Save_LCOE_Achieved"]


    # Colors for each energy source
    colors = ['skyblue', 'seagreen', 'coral', 'purple']

    fig, ax = plt.subplots()

    # Stacking wind, wave, and ocean current energy
    ax.bar(Save_LCOE_Achieved, Save_TotalMWAvgWind, width=3, color=colors[0], label='Wind Avg Gen')
    ax.bar(Save_LCOE_Achieved, Save_TotalMWAvgWave, width=3, color=colors[1], label='Wave Avg Gen', bottom=Save_TotalMWAvgWind)
    ax.bar(Save_LCOE_Achieved, Save_TotalMWAvgKite, width=3, color=colors[2], label='Kite Avg Gen', bottom=Save_TotalMWAvgWind + Save_TotalMWAvgWave)

    # Adding curtailment below the x-axis
    ax.bar(Save_LCOE_Achieved, -Save_totalMWAvgCurtailment, width=3, color=colors[3], label='Curtailment')

    ax.plot(Save_LCOE_Achieved, SaveTotalMWAvg, color='black', label='Avg Gen to Shore')

    # Setting labels and title
    ax.set_xlabel('LCOE ($/MWh)')
    ax.set_ylabel('Avg Power (MW Avg)')
    ax.set_title('Energy Generation and Curtailment by Source')
    plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0., frameon=False)

    plt.grid(linestyle='--')

    if SavePath!=None:
        plt.savefig(SavePath, bbox_inches='tight', dpi=700)
        
#Plot the turbine locations for a given solution
def PlotTurbineLocations(SolutionPath, Legend, PathDataUnderLayer, LegendUnderLayerColorbar, StateCountoursPath, UnderLayerVariable=["Energy_pu", "RawResource","RawResource"],
    LCOE_Target=-1, LatMaxMin=(33.3, 37.2), LongMaxMin=(-78.7, -74.3), SavePath=None):

    #Legend=[] #one list per solution, each list contains the legend for each design of each technology 3d list
    #UB_LCOE: Plot will use the LCOE closes to this value to plot the turbine locations

    ShapeFileCoast=StateCountoursPath+"ne_10m_coastline.shp"
    ShapeFileStates=StateCountoursPath+"ne_10m_admin_1_states_provinces_lines.shp"

        
    df = gpd.read_file(ShapeFileCoast)
    df1 = gpd.read_file(ShapeFileStates)

    fig, ax = plt.subplots()
    fig.set_size_inches(18.5, 10.5)
    

    df.plot(color='black',linewidth=1,ax=ax)
    df1.plot(color='black',linewidth=1,ax=ax)

    LatLongList_All_Wind=[]
    LatLongList_All_Wave=[]
    LatLongList_All_OC=[]

    ColorList=["black","red","blue","green","orange","purple"]

    LegendIn=[]
    for i, path in enumerate(SolutionPath):
        data=np.load(path,allow_pickle=True)
        SolutionLCOE=data["Save_LCOE_Achieved"][:-1].astype(float)
        SolutionIDX=np.abs(LCOE_Target - SolutionLCOE).argmin()
        
        
        if len(data["PathWindDesigns"])>0:
            RangeOfRotations=[0, 180, 90, 45, -45, -90]
            
            Y_Wind=data["Save_Y_Wind"][SolutionIDX]
            
            LatLongWind=np.empty((0, 2))
            IDX_Designs=[]
            IDX_Designs.append(0)
            for PathWindDesigns in data["PathWindDesigns"]:

                WindData=np.load(PathWindDesigns,allow_pickle=True)
                LatLongWind=np.concatenate((LatLongWind,WindData["LatLong"]))
                IDX_Designs.append(IDX_Designs[-1]+WindData["LatLong"].shape[0])
            
            
            for j in range(len(IDX_Designs)-1):
                LatLong_tmp=LatLongWind[IDX_Designs[0+j]:IDX_Designs[1+j],:][Y_Wind[IDX_Designs[0+j]:IDX_Designs[1+j]]>0,:]
                #if len(LatLong_tmp)>0:
                plt.scatter(LatLong_tmp[:,1], LatLong_tmp[:,0], marker=(3, 0, RangeOfRotations[i+j]), s=100, edgecolors=ColorList[i+j], facecolors='none',alpha=1, label=Legend[i][0][j])


        if len(data["PathWaveDesigns"])>0:
            RangeOfRotations=[45, 0, 30, 60]
            
            Y_Wave=data["Save_Y_Wave"][SolutionIDX]

            LatLongWave=np.empty((0, 2))
            IDX_Designs=[]
            IDX_Designs.append(0)
            for PathWaveDesigns in data["PathWaveDesigns"]:

                WaveData=np.load(PathWaveDesigns,allow_pickle=True)
                LatLongWave=np.concatenate((LatLongWave,WaveData["LatLong"]))
                IDX_Designs.append(IDX_Designs[-1]+WaveData["LatLong"].shape[0])
            
            for j in range(len(IDX_Designs)-1):
                LatLong_tmp=LatLongWave[IDX_Designs[0+j]:IDX_Designs[1+j],:][Y_Wave[IDX_Designs[0+j]:IDX_Designs[1+j]]>0,:]
                #if len(LatLong_tmp)>0:
                plt.scatter(LatLong_tmp[:,1], LatLong_tmp[:,0], marker=(4, 0, RangeOfRotations[i+j]), s=100, edgecolors=ColorList[i+j], facecolors='none',alpha=1, label=Legend[i][1][j]) 

                

        if len(data["PathKiteDesigns"])>0:
            RangeOfRotations=["1", "2", "3", "4","+"]
            
            Y_Kite=data["Save_Y_Kite"][SolutionIDX]
            LatLongKite=np.empty((0, 2))
            IDX_Designs=[]
            IDX_Designs.append(0)
            for PathKiteDesigns in data["PathKiteDesigns"]:

                KiteData=np.load(PathKiteDesigns,allow_pickle=True)
                LatLongKite=np.concatenate((LatLongKite,KiteData["LatLong"]))
                IDX_Designs.append(IDX_Designs[-1]+KiteData["LatLong"].shape[0])
            
            for j in range(len(IDX_Designs)-1):
                
                LatLong_tmp=LatLongKite[IDX_Designs[0+j]:IDX_Designs[1+j],:][Y_Kite[IDX_Designs[0+j]:IDX_Designs[1+j]]>0,:]
                # if len(LatLong_tmp)>0:
                plt.scatter(LatLong_tmp[:,1], LatLong_tmp[:,0], marker=RangeOfRotations[j+i], c=ColorList[i+j], s=100,  alpha=0.9, label=Legend[i][2][j]) 

            
    ax.set_xlim(LongMaxMin) 
    ax.set_ylim(LatMaxMin)
    plt.xlabel("Longitude", fontsize=12)
    plt.ylabel("Latitude", fontsize=12)

    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), fancybox=True, shadow=True, ncol=5, fontsize=12,frameon=False)

    divider = make_axes_locatable(ax)

    if len(PathDataUnderLayer)!=0:
        for i, PathData in enumerate(PathDataUnderLayer):
            data=np.load(PathData,allow_pickle=True)
            if UnderLayerVariable[i]=="Energy_pu":
                RefData=data["Energy_pu"].mean(axis=0)
            else:
                RefData=data[UnderLayerVariable[i]]
                if len(np.shape(RefData))==2:
                    RefData=RefData.mean(axis=0)
                
            LatLong=data["LatLong"]
            
            if data["ResolutionKm"]!=-1:
                s=1 #2km
            else:
                s=10

            norm = mcolors.Normalize(vmin=min(RefData), vmax=max(RefData))

            
            if i==0:
                sm = cm.ScalarMappable(norm=norm, cmap='jet')
                colors = plt.cm.jet(norm(RefData))
                ax.scatter(LatLong[:,1], LatLong[:,0], c=sm.to_rgba(RefData), s=s, alpha=0.3)
            
            else:
                sm = cm.ScalarMappable(norm=norm, cmap='jet')
                colors = plt.cm.jet(norm(RefData))
                ax.scatter(LatLong[:,1], LatLong[:,0], edgecolor=sm.to_rgba(RefData), facecolors='none', marker='s', s=s+5,alpha=0.7)


            clb=plt.colorbar(sm, norm=norm, cax=divider.append_axes("right", size="5%", pad=0.05 + i*0.6))
            clb.ax.set_title(LegendUnderLayerColorbar[i])
                
                
    if SavePath!=None:
        plt.savefig(SavePath, bbox_inches='tight', dpi=700)
        
