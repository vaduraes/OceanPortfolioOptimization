#Read Grib2 File from WWIII (Run on LINUX or MAC because of pygrib)
#Still needs to treat missing values

import numpy as np
import pygrib
import datetime as dt
from datetime import datetime
from datetime import timedelta
from dateutil.relativedelta import relativedelta


file="WWW3Data/multi_1.at_4m.hs.200502.grb2"
grb = pygrib.open(file)
Grb_Lat=grb[1]['latitudes']
Grb_Long=grb[1]['longitudes']
grb.close()

LatMax=36.7
LatMin=33.7

LongMax=-74.5+360
LongMin=-78.7+360

#1st filter for valid lat long values
ValidLatIdx=(Grb_Lat<=LatMax)*(Grb_Lat>=LatMin)
ValidLongIdx=(Grb_Long<=LongMax)*(Grb_Long>=LongMin)
ValidIdx=(ValidLatIdx*ValidLongIdx)

LatLong=[Grb_Lat[ValidIdx],Grb_Long[ValidIdx]]
LatLong=np.array(LatLong).T
NumSites=LatLong.shape[0]


FromDate=dt.date(2009, 1, 1)
ToDate=dt.date(2013, 12, 31)


Date=FromDate - relativedelta(months=1)


while Date!=(ToDate+timedelta(days=1)):

    Date = Date + relativedelta(months=1)
    DateStr=Date.strftime('%Y%m')
    print("The File Date %s is Being Processed"%DateStr)
    
    file_hs="WWW3Data/multi_1.at_4m.hs."+DateStr+".grb2"
    grb_hs = pygrib.open(file_hs)
    NumTimeIntervals=len(grb_hs.select())
    Hs_Month=np.zeros((NumSites ,NumTimeIntervals),dtype=float)
    
    for TimeIdx in range (1,NumTimeIntervals+1,1):
        Temp_hs=np.reshape(grb_hs[TimeIdx]['values'].data,-1)
        Hs_Month[:,TimeIdx-1]=Temp_hs[ValidIdx]

        
    grb_hs.close()   
    
    
    file_tp="WWW3Data/multi_1.at_4m.tp."+DateStr+".grb2"
    grb_tp = pygrib.open(file_tp)
    Tp_Month=np.zeros((NumSites ,NumTimeIntervals),dtype=float)
    
    for TimeIdx in range (1,NumTimeIntervals+1,1):
        Temp_tp=np.reshape(grb_tp[TimeIdx]['values'].data,-1)
        Tp_Month[:,TimeIdx-1]=Temp_tp[ValidIdx]
    
    grb_tp.close()   
     
    if Date==FromDate:
        Hs=Hs_Month
        Tp=Tp_Month
        
    else:
        Hs=np.concatenate((Hs, Hs_Month),axis=1)
        Tp=np.concatenate((Tp, Tp_Month),axis=1)
    
SecondIdxFilter1=(np.sum((Hs==9999),axis=1)!=Hs.shape[1])

Hs=Hs[SecondIdxFilter1,:]
Tp=Tp[SecondIdxFilter1,:]
LatLong=LatLong[SecondIdxFilter1,:]

np.savez("WWIII_RawData.npz", LatLong=LatLong, Hs=Hs, Tp=Tp, FromDate=FromDate, ToDate=ToDate)
