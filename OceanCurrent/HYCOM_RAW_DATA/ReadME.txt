1)To fill missing values of HYCOM and save the results in .npz format run:

PreprocessHYCOM.py

Inputs:
#.nc files with the hycom data. It is possible to donwload all hycom though their FTP, however this data is very large.
here we used their API to get the data. It is important to mention that their API has a limit in the size of the file you can request, what justifies the 
several .nc files (P1, ..., P13) we currently have in the HYCOM_RAW_DATA folder 

Outputs:
#.npz current speed, time and latlong for all hycom data (CurrentSpeedHYCOM2009_2014.npz)