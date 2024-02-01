#py37 env
import numpy as np
from datetime import datetime, timedelta
import sys
import numpy as np
import datetime as dt
from dateutil.relativedelta import relativedelta
import os
import urllib.request
from multiprocessing import Pool
from tqdm import tqdm  
import datetime as dt
import matlab.engine
import multiprocessing
from joblib import Parallel, delayed
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
from scipy.io import savemat
import random

sys.path.append('./Tools')
from OC_nc2npz import ReadHYCOMData
from GeneralGeoTools import GetDepth
sys.path.append('./Tools/KitesMatlab')

def WriteMatlabInputs2Kites (HycomPath, StartDTime, EndDTime, TimeDiscretization, DepthDataPath="./InputData", SaveMatPath=None):
    Data=np.load(HycomPath,allow_pickle=True)
    OCSpeed=Data["OCSpeed"]
    LatLong=Data["LatLong"]
    depth=Data["depth"]
    TimeList=Data["TimeList"]
    IdxTimeIn=(TimeList>=StartDTime) * (TimeList<=EndDTime)
    OCSpeed=OCSpeed[IdxTimeIn,:,:]
    TimeList=TimeList[IdxTimeIn]


    NewTimeList=[]
    #Get Max Depth
    Dmax=GetDepth(DepthDataPath, LatLong)
    IdxIn=(Dmax>10)* (Dmax<=2500) #only use data deeper than 10m and less than 2500m
    # IdxInDepth= (depth>10)*(depth<=2500) #Dont need data shallower than 10m and deeper than 2500m
    IdxInDepth= (depth<=2500) #Dont need data deeper than 2500m
    depth=depth[IdxInDepth]
    Dmax=Dmax[IdxIn]
    OCSpeed=OCSpeed[:,:,IdxIn]
    OCSpeed=OCSpeed[:,IdxInDepth,:]

    LatLong=LatLong[IdxIn,:]

    #change time discretization
    print("Changing time discretization")
    BaseTimeDiscretization=3 #3hours at each time discretization from hycom
    TimeStep=int(TimeDiscretization/BaseTimeDiscretization)

    UpdatedOCSpeed=np.ones((OCSpeed.shape[2],OCSpeed.shape[1],int(OCSpeed.shape[0]/TimeStep)))*-1 #Site, Depth, Time
    for i in tqdm(range(0,OCSpeed.shape[0],TimeStep)):
        for d in range(OCSpeed.shape[1]):
            for s in range(OCSpeed.shape[2]):
                
                ValidIdxs=(OCSpeed[i:i+TimeStep,d,s]!=-1)
                if np.sum(ValidIdxs)>0:
                    UpdatedOCSpeed[s, d, int(i/TimeStep)] = np.mean(OCSpeed[i:i+TimeStep,d,s][ValidIdxs])
            

        NewTimeList.append(TimeList[i])
    NewTimeList=np.array(NewTimeList)

    IdxFilterVeryLowMedians=np.max(np.median(UpdatedOCSpeed,axis=2),axis=1)>=0.15 #Maximum median across all depths for each site location. Filter out sites with very low median currents
    UpdatedOCSpeed=UpdatedOCSpeed[IdxFilterVeryLowMedians,:,:]
    LatLong=LatLong[IdxFilterVeryLowMedians,:]
    Dmax=Dmax[IdxFilterVeryLowMedians]


    #Fill Missing Data With Previous Data
    CountFillTimes=0
    print("Filling missing data")
    for s in tqdm(range(UpdatedOCSpeed.shape[0])):
        for t in range(1,UpdatedOCSpeed.shape[2],1):
            if UpdatedOCSpeed[s,0,t]==-1:
                UpdatedOCSpeed[s,:,t]=UpdatedOCSpeed[s,:,t-1]
                CountFillTimes+=1
    print("Percentage of missing data filled: ", CountFillTimes/(UpdatedOCSpeed.shape[0]*UpdatedOCSpeed.shape[2])*100)
    
    #Plot a sample Time Series
    print("Plotting a sample time series")  
    plt.figure()
    RandSite=random.randint(0,UpdatedOCSpeed.shape[0]-1)
    plt.plot(NewTimeList,UpdatedOCSpeed[RandSite,0,:])
    plt.xticks(rotation=90)
    plt.title("Sample Time Series at Lat: "+str(LatLong[RandSite,0])+" Long: "+str(LatLong[RandSite,1]))
    
    plt.xlabel("Time")
    plt.ylabel("Current Speed (m/s)")
    
    if SaveMatPath is not None:
        Data={"LatLong":LatLong, "TimeList":[NewTimeList[i].strftime("%Y%m%d") for i in range(len(NewTimeList))], "OCSpeed":UpdatedOCSpeed, "depth":depth, "Dmax":Dmax}
        savemat(SaveMatPath+"OCSpeedHycom_"+StartDTime.strftime("%Y%m%d")+"_"+EndDTime.strftime("%Y%m%d")+".mat", Data)
        
    else:
        return LatLong, NewTimeList, UpdatedOCSpeed, depth, Dmax
    
    
def TimeSeriesGeneration_Kite (PathKiteParams, i_vd, i_cs, SavePowerTimeSeriesPath, FullMatlabHycomDataPath, NumEnvs=5, MatlabKitePath='./Tools/KitesMatlab/'):
    BaseCSpeed=[0.5, 1 , 1.5 , 2.0 , 2.5]
    VerticalDepth=[50, 100, 150, 200]

    dataframe = pd.read_excel(PathKiteParams)
    Power=dataframe.iloc[1:5,1:].values
    StructuralMass=dataframe.iloc[1+7:5+7,1:].values
    Span=dataframe.iloc[1+14:5+14,1:].values
    AspectRatio=dataframe.iloc[1+21:5+21,1:].values
    Length=dataframe.iloc[1+28:5+28,1:].values
    Diameter=dataframe.iloc[1+35:5+35,1:].values
    LCOE=dataframe.iloc[1+42:5+42,1:].values
    AnnualizedCost=LCOE*Power*365*24 #Annualized Cost (in $/year)


    Envs=[matlab.engine.start_matlab() for i in range(NumEnvs)] #Create matlab envs

    for i in range(NumEnvs):
        Envs[i].eval("addpath('{}');".format(MatlabKitePath) ,nargout=0)
        Envs[i].eval("addpath('{}');".format(MatlabKitePath+"functions/") ,nargout=0)
        Envs[i].eval("addpath('{}');".format(MatlabKitePath+"scripts/") ,nargout=0)
        Envs[i].eval("addpath('{}');".format(MatlabKitePath+"functions/SWDT/") ,nargout=0)
        # Envs[i].eval("addpath('{}');".format(SaveMatPath) ,nargout=0)

        Envs[i].eval("load('{}');".format(FullMatlabHycomDataPath),nargout=0) #Load Ocean Current Data


    Envs[0].eval('[Size_s, SizeDepth, SizeTime]=size(OCSpeed);' ,nargout=0) #Get Shape of input data
    Size_s_matlab =int(Envs[0].workspace["Size_s"])
    Size_depth_matlab =int(Envs[0].workspace["SizeDepth"])
    Size_time_matlab =int(Envs[0].workspace["SizeTime"])
    LatLong=np.array(Envs[0].workspace["LatLong"])
    dmax=np.array(Envs[0].workspace["Dmax"][0])[0,:]

    #Create list of commands
    print("Running Design: Vertical Depth Index: %d --- Base Current Speed Index: %d" % (i_vd, i_cs))
    #Create list of commands
    cmd1=[]
    cmd2=[]        
    for i in range(Size_s_matlab):
        if dmax[i]>30: #only consider sites with depth above 30m
            
            Envs[0].eval('MaxMeanDepthCurrentAtSite=max((mean(OCSpeed({},:,:),3)),[],"all");'.format(i+1) ,nargout=0) #(mean for each depth, then max across depths)
            MaxMeanDepthCurrentAtSite =Envs[0].workspace["MaxMeanDepthCurrentAtSite"]
            if MaxMeanDepthCurrentAtSite>0.5:


                
                RatedPower_tmp=Power[i_vd, i_cs]    #Rated Power (in Kw)
                Span_tmp = Span[i_vd, i_cs]         #wingspan (in m) 
                AR_tmp = AspectRatio[i_vd, i_cs]    #wing aspect ratio 
                Len_tmp = Length[i_vd, i_cs]        #fuselage length (in m)
                Dia_tmp = Diameter[i_vd, i_cs]      #fuselage diameter (in m)


                cmd1_tmp = str('iSite ='+str(i+1)+';')
                cmd2_tmp = str('uGeo =['+str(Span_tmp)+','+str(AR_tmp)+','+str(Len_tmp)+','+str(Dia_tmp)+','+str(RatedPower_tmp)+'];')
                cmd1.append(cmd1_tmp)
                cmd2.append(cmd2_tmp)

    #Run all site locations for a given design
    G_count=0
    NumRuns=0
    PowerTimeSeries=[]
    MatlabSiteIdx=[]
    while G_count!=len(cmd1):
        G_count=G_count+NumRuns

        print("%.2f %% complete" % (G_count/len(cmd1)*100))

        NumRuns=np.min([len(Envs),len(cmd1)-G_count]) #Number of runs to do in parallel

        #Insert Data
        for i in range(0,NumRuns,1):
            Envs[i].eval(cmd1[G_count+i],nargout=0) #Inputs
            Envs[i].eval(cmd2[G_count+i],nargout=0) #Inputs

        #Separate loop to run in parallel (only for running in parallel)
        for i in range(0,NumRuns,1):
            Envs[i].powerFunc_PythonHycom(nargout=0,background=True) #Solve problem in background

        for i in range(0,NumRuns,1):   
            Jopt_vec =np.array(Envs[i].workspace["Jopt_vec"][0])[0,:]
            PowerTimeSeries.append(Jopt_vec)
            MatlabSiteIdx.append(int(Envs[i].workspace["iSite"]))
            
    PowerTimeSeries=np.array(PowerTimeSeries)
    MatlabSiteIdx=np.array(MatlabSiteIdx)

    np.savez(SavePowerTimeSeriesPath+'PowerTimeSeriesKite_VD'+str(VerticalDepth[i_vd])+'_BCS'+str(BaseCSpeed[i_cs])+'.npz', PowerTimeSeries=PowerTimeSeries, MatlabSiteIdx=MatlabSiteIdx,
                LatLong=LatLong, dmax=dmax, AnnualizedCost=AnnualizedCost[i_vd, i_cs], StructuralMass=StructuralMass[i_vd, i_cs], RatedPower=Power[i_vd, i_cs],
                Span=Span[i_vd, i_cs], AspectRatio=AspectRatio[i_vd, i_cs], Length=Length[i_vd, i_cs], Diameter=Diameter[i_vd, i_cs])

    #Close environments
    for i in range(NumEnvs):
        Envs[i].quit()