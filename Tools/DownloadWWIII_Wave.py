# #Download significant wave height and wave period from Wave watch III model
from dateutil.relativedelta import relativedelta
import datetime as dt
import urllib.request
import time
import os


def DonwloadWave_WWIII(FromDate, ToDate, SavePath):
    NumberOfFiles=(ToDate-FromDate).days+1

    date=FromDate
    attempts=0
    while date!=(ToDate + relativedelta(months=1)):
        try:
            attempts=0
            StrDate=date.strftime("%Y%m")
            year=date.strftime("%Y")
            month=date.strftime("%m")
        
            path = "https://polar.ncep.noaa.gov/waves/hindcasts/multi_1/" + year + month + '/gribs/'

            filenameTp1 = 'multi_1.at_4m.tp.'+StrDate+'.grb2'
            filenameTp2 = 'multi_1.at_4m.tp.'+StrDate+'.grb2.MD5'
            
            filenameHs1 = 'multi_1.at_4m.hs.'+StrDate+'.grb2'
            filenameHs2 = 'multi_1.at_4m.hs.'+StrDate+'.grb2.MD5'
            
            #Afeter 07/2013 some files are missing (skip this all together)
            # if not os.path.exists(SavePath+filenameTp2):
            #     urllib.request.urlretrieve(path + filenameTp2, SavePath+filenameTp2) #very small file
            # if not os.path.exists(SavePath+filenameHs2):
            #     urllib.request.urlretrieve(path + filenameHs2, SavePath+filenameHs2)
                
            # Only download if the file does not exist
            if not os.path.exists(SavePath+filenameTp1):
                urllib.request.urlretrieve(path + filenameTp1, SavePath+filenameTp1)
            else:
                print('---- File %s already exists, skipping'%(filenameTp1))
            
            if not os.path.exists(SavePath+filenameHs1): 
                urllib.request.urlretrieve(path + filenameHs1, SavePath+filenameHs1)
            else:
                print('---- File %s already exists, skipping'%(filenameHs1))
            
            print('---- Year %s /Month %s Done,    Percentage Complete %s'%(year,(month),((date-FromDate).days*100/NumberOfFiles)))
            date=date + relativedelta(months=1)

        
        except:
            attempts=attempts+1
            print('---- %s or %s Failed'%(filenameTp1, filenameHs1))
            print('---- Waiting 3 min to try again, after 3 attempts this file will be skipped')
            
            time.sleep(3*60)
            if attempts==3:
                date=date + relativedelta(months=1)
                print('---- Year (%s) /Month (%s) Skipped'%(year,month))
                attempts=0