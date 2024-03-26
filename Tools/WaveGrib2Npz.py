#Run on py37 (pygrib works on windows but it is tricky to install)

import numpy as np
from datetime import datetime, timedelta
import sys
import numpy as np
import datetime as dt
from dateutil.relativedelta import relativedelta
import xarray as xr
import pygrib
from tqdm import tqdm
from multiprocessing import Pool
from functools import partial
from GeneralGeoTools import GetDepth, GetDistanceToShore

#Process a single data for the WWIII Grib2 data
def ProcessSpecificDateGRIB2WWIII(Date, ValidLatLongIdx):
    NumSites=np.sum(ValidLatLongIdx)
    DateStr=Date.strftime('%Y%m')
    print("The File Date %s is Being Processed"%DateStr)
    
    file_hs="./InputData/Wave/WWIII/multi_1.at_4m.hs."+DateStr+".grb2"
    grb_hs = pygrib.open(file_hs)
    NumTimeIntervals=len(grb_hs.select())
    Hs_Month=np.zeros((NumSites ,NumTimeIntervals),dtype=float)
    
    for TimeIdx in range (1,NumTimeIntervals+1,1):
        Temp_hs=np.reshape(grb_hs[TimeIdx]['values'].data,-1)
        Hs_Month[:,TimeIdx-1]=Temp_hs[ValidLatLongIdx]

        
    grb_hs.close()   
    
    
    file_tp="./InputData/Wave/WWIII/multi_1.at_4m.tp."+DateStr+".grb2"
    grb_tp = pygrib.open(file_tp)
    Tp_Month=np.zeros((NumSites ,NumTimeIntervals),dtype=float)
    
    for TimeIdx in range (1,NumTimeIntervals+1,1):
        Temp_tp=np.reshape(grb_tp[TimeIdx]['values'].data,-1)
        Tp_Month[:,TimeIdx-1]=Temp_tp[ValidLatLongIdx]
    
    grb_tp.close()   
    
    StartDTime=datetime(Date.year, Date.month, Date.day, 0, 0, 0)
    DateTimeList=[]
    for i in range(NumTimeIntervals):
       DateTimeList.append(StartDTime+timedelta(hours=3*i))
    
    
    return Hs_Month, Tp_Month, DateTimeList

 
def ConvertGrib2NPZ_WWIII(FromDate, ToDate, SavePath=None, LatMinMax=(33, 37), LongMinMax=(-81,-73), DepthMinMax=(1,1000), InputDataPath="./InputData"):
   Date=FromDate
   Date_i=FromDate
   ListofDates=[]
   while Date_i<=ToDate:
      ListofDates.append(Date_i)
      Date_i=Date_i+relativedelta(months=+1)
      
   LongMinMax=np.array(LongMinMax)
   LongMinMax=LongMinMax+360



   StrDate=Date.strftime("%Y%m")
   year=Date.strftime("%Y")
   month=Date.strftime("%m")

   filenameTp1 = 'multi_1.at_4m.tp.'+StrDate+'.grb2'
   filenameHs1 = 'multi_1.at_4m.hs.'+StrDate+'.grb2'

               
   # Specify the path to your GRIB2 file
   grib2_file_pathTp1 = './InputData/Wave/WWIII/'+filenameTp1
   grib2_file_pathHs1 = './InputData/Wave/WWIII/'+filenameHs1

   # Read GRIB2 file
   grb = pygrib.open(grib2_file_pathHs1)
   Grb_Lat=grb[1]['latitudes']
   Grb_Long=grb[1]['longitudes']
   grb.close()

   #1st filter for valid lat long values
   ValidLatIdx=(Grb_Lat<=LatMinMax[1])*(Grb_Lat>=LatMinMax[0])
   ValidLongIdx=(Grb_Long<=LongMinMax[1])*(Grb_Long>=LongMinMax[0])
   ValidLatLongIdx=(ValidLatIdx*ValidLongIdx)

   LatLong=[Grb_Lat[ValidLatLongIdx],Grb_Long[ValidLatLongIdx]]
   LatLong=np.array(LatLong).T
   NumSites=LatLong.shape[0]

   Iterator=partial(ProcessSpecificDateGRIB2WWIII,ValidLatLongIdx=ValidLatLongIdx)

   #Run in parallel
   # if __name__ == '__main__':
   try:
      with Pool(24) as p:
         ListOfGrib2Results = list(tqdm(p.imap(Iterator, ListofDates), total=len(ListofDates)))
         #Hs_Month, Tp_Month, Date=ProcessSpecificDateGRIB2WWIII(Date,ValidLatLongIdx)
   except:
      with Pool(4) as p:
         ListOfGrib2Results = list(tqdm(p.imap(Iterator, ListofDates), total=len(ListofDates)))
         #Hs_Month, Tp_Month, Date=ProcessSpecificDateGRIB2WWIII(Date,ValidLatLongIdx)

   #Concatenate for Hs and Tp
   DateTimeList=[]
   for i in range(len(ListOfGrib2Results)):
      DateTimeList=DateTimeList+ListOfGrib2Results[i][2]
   Hs=np.concatenate([i[0] for i in ListOfGrib2Results],axis=1) #Location, time
   Tp=np.concatenate([i[1] for i in ListOfGrib2Results],axis=1)

   LatLong[:,1]=LatLong[:,1]-360



   #Remove duplicate dates. End of one month and start of another month may have the same date
   index_dict = {}
   unique_indices = []

   # Iterate over the list using enumerate
   for index, value in enumerate(DateTimeList):
      
      if value not in index_dict:
         index_dict[value] = index
         unique_indices.append(index)

   DateTimeList=np.array(DateTimeList)[unique_indices]

   Hs=Hs[:, unique_indices]
   Tp=Tp[:, unique_indices]

   
   DistanceShore=GetDistanceToShore(InputDataPath, LatLong)
   Depth=GetDepth(InputDataPath, LatLong)
   
   IdxIn=(Depth>=DepthMinMax[0]) * (Depth<=DepthMinMax[1])
   Depth=Depth[IdxIn]
   DistanceShore=DistanceShore[IdxIn]
   LatLong=LatLong[IdxIn,:]
   Hs=Hs[IdxIn,:]
   Tp=Tp[IdxIn,:]
   
   #exclude sites with more than 0.1%missing values
   SecondIdxFilter1=(np.sum(Hs>999,axis=1)<=Hs.shape[0]*0.1/100) * (np.sum(Tp>999,axis=1)<=Tp.shape[0]*0.1/100)
   Hs=Hs[SecondIdxFilter1,:]
   Tp=Tp[SecondIdxFilter1,:]
   LatLong=LatLong[SecondIdxFilter1,:]
   Depth=Depth[SecondIdxFilter1]
   DistanceShore=DistanceShore[SecondIdxFilter1]
   
   #Fill the missing values with the value in the previous 3h
   #We used this simple procedure because there are very few missing values
   for i in range(Tp.shape[0]):
      for j in range(Tp.shape[1]):
         if Tp[i,j]==9999:
               Tp[i,j]=Tp[i,j-1]
               
         if Hs[i,j]==9999:
               Hs[i,j]=Hs[i,j-1]

   if SavePath!=None:
      np.savez(SavePath, LatLong=LatLong, Hs=Hs, Tp=Tp, DateTimeList=DateTimeList, DistanceShore=DistanceShore, Depth=Depth)
   else:
      return LatLong, Hs, Tp, DateTimeList, DistanceShore, Depth