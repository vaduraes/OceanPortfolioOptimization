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
from functools import partial


LatMinMax=(33.3, 37.2)
LongMinMax=(-78.7, -74.3)
    
def ListHYCOM_MultiProcess(StartDTime, EndDTime, LatMinMax=LatMinMax,  LongMinMax=LongMinMax):
    TimeSteps=15
    List_DownloadLinks=[]
    ListSavePaths=[]
    ListStartDate=[]
    ListEndDate=[]

    SDate=StartDTime
    EDate=SDate
    while EDate<=EndDTime:
        
        Year=SDate.year
        EDate=SDate+timedelta(hours=3*TimeSteps)
        
        if EDate.year!=Year:
            EDate=datetime(Year, 12, 31, 21, 0, 0)
        
        TimeStart_Str=SDate.strftime("%Y-%m-%dT%H")
        TimeEnd_Str=(EDate).strftime("%Y-%m-%dT%H")
        

        SavePath="./InputData/OceanCurrent/Hycom/"
        
        LinkDownload_p1="https://ncss.hycom.org/thredds/ncss/GLBv0.08/expt_53.X/data/"+str(Year)\
            +"?var=water_u&var=water_v&north="+str(LatMinMax[1])+"&west="+str(LongMinMax[0])+"&east="+str(LongMinMax[1])+"&south="+str(LatMinMax[0])\
            +"&horizStride=1&time_start="+TimeStart_Str+"%3A00%3A00Z&time_end="+TimeEnd_Str+"%3A00%3A00Z&timeStride=1&vertCoord=&addLatLon=true&accept=netcdf"

        SaveFileName=TimeStart_Str+"_"+TimeEnd_Str+".nc"
        # urllib.request.urlretrieve(LinkDownload_p1, SavePath+SaveFileName)
        
        #use list for multiprocessing on download
        if os.path.isfile(SavePath+SaveFileName)==False:
            List_DownloadLinks.append(LinkDownload_p1)
            ListSavePaths.append(SavePath+SaveFileName)
            ListStartDate.append(SDate)
            ListEndDate.append(EDate)
        
        SDate=EDate+timedelta(hours=3)
        # print("Data from "+TimeStart_Str+" to "+TimeEnd_Str+" downloaded")
        # print("Percentage Completed: "+str(round((EDate-StartDTime).total_seconds()/(EndDTime-StartDTime).total_seconds()*100, 2))+"%")
    return List_DownloadLinks, ListSavePaths, ListStartDate, ListEndDate
        
def DonwloadHYCOM(i, List_DownloadLinks,ListSavePaths):


    try:
        print("Downloading File: "+ListSavePaths[i],flush=True)
        urllib.request.urlretrieve(List_DownloadLinks[i], ListSavePaths[i])
        
    except:
        pass
    
def DownloadHYCOM_MultiProcess (StartDTime, EndDTime, PoolSize=10, LatMinMax=LatMinMax,  LongMinMax=LongMinMax):
    List_DownloadLinks, ListSavePaths, ListStartDate, ListEndDate=ListHYCOM_MultiProcess(StartDTime, EndDTime, LatMinMax=LatMinMax,  LongMinMax=LongMinMax)
    
    Iterator=partial(DonwloadHYCOM,List_DownloadLinks=List_DownloadLinks,ListSavePaths=ListSavePaths)
    
    with Pool(PoolSize) as p:
        ListOfGrib2Results = list(tqdm(p.imap(Iterator, range(len(List_DownloadLinks))), total=len(List_DownloadLinks)))