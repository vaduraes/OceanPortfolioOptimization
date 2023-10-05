#Download significant wave height and wave period from Wave watch III model

from dateutil.relativedelta import relativedelta
import datetime as dt
import ftplib

FromDate=dt.date(2005, 2, 1)
ToDate=dt.date(2019, 1, 1)
NumberOfFiles=(ToDate-FromDate).days+1

ftp = ftplib.FTP("polar.ncep.noaa.gov") 
ftp.login("", "") 

date=FromDate
while date!=(ToDate + relativedelta(months=1)):
    
    StrDate=date.strftime("%Y%m")
    year=date.strftime("%Y")
    month=date.strftime("%m")
    
    
    path = '/pub/history/waves/multi_1/'+year+month+'/gribs/'
    ftp.cwd(path)
    
    filenameTp1 = 'multi_1.at_4m.tp.'+StrDate+'.grb2'
    filenameTp2 = 'multi_1.at_4m.tp.'+StrDate+'.grb2.MD5'
    
    filenameHs1 = 'multi_1.at_4m.hs.'+StrDate+'.grb2'
    filenameHs2 = 'multi_1.at_4m.hs.'+StrDate+'.grb2.MD5'

    
    LocalFileTp1= open("E:\OceanProject\WaveData\WWW3Data\\"+filenameTp1, 'wb')
    #LocalFileTp2= open("E:\OceanProject\WaveData\WWW3Data\\"+filenameTp2, 'wb')
    
    LocalFileHs1= open("E:\OceanProject\WaveData\WWW3Data\\"+filenameHs1, 'wb')
    #LocalFileHs2= open("E:\OceanProject\WaveData\WWW3Data\\"+filenameHs2, 'wb')
    
    ftp.retrbinary('RETR '+filenameTp1, LocalFileTp1.write)
    #ftp.retrbinary('RETR '+filenameTp2, LocalFileTp2.write)
    
    ftp.retrbinary('RETR '+filenameHs1, LocalFileHs1.write)
    #ftp.retrbinary('RETR '+filenameHs2, LocalFileHs2.write)
    
    LocalFileTp1.close()
    #LocalFileTp2.close()
    LocalFileHs1.close()
    #LocalFileHs2.close()
    
    print('---- Month %s Done'%(month))
    date=date + relativedelta(months=1)
    
ftp.quit()

