1) To download wave data from WWIII model FTP run:

DownloadDataFTPWWIII.py

Inputs: 
#FromDate / ToDate> Range of days you want to download the data 

Outputs:
#A series of .grb2 files that we will process later

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

2) To convert the .grb2 files to .npz run:

Grib_to_NPZ_WWIII.py

Inputs:
#FromDate / ToDate> Range of days you want save as NPZ file
#Series of .grb2 files you downloaded in the previous step

Outputs:
#.npz file with latitude langitude for a series of site locations as well as its significant wave height and period
#in this output the data is not yet "cleaned", containing a few missing values

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

3) To generate a clean .npz file (no missing values)  containing wave data run.

CleanDataWWII_NPZ.py

Inputs:
#WWIII_Data: .npz file obtained from step (2)

Outputs:
#.npz file containing information of significant wave height and period, as well as latitude and longitude for all grid cells with its respective depth and 
minimum distance from the shore (WaveHsTp_WWIII.npz)


--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

4) To convert significant have height and periods to energy values run:

ComputeWaveEnergy.py

Inputs:
#TurbineFile: .csv file containing data of the turbine model used
#WaveData: .npz file containing data of historical wave height and wave period

Outputs:
#.npz file containing the energy generated for each grid cell in the entire NC coast where the turbine can be deployed.

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

5) To generate geo plots with the CF distribution of the wave resource run:

WaveCFGeoPlot.py

Inputs:
#files in GEO_data describing US contours
#HourlyData: Data of CF for each site location of the wave energy resource

Outputs:
#Plot of the spacial distribution of the wave resource with its CF

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

3) To concatenate annualized cost information in the .npz file constructed in (4) run:

ConcatenateCost.py

Input:
#.npz file of step (1)

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

7) To generate an geo plot of the LCOE of the ocean current resource in the NC region run:

LCOEGeoPlot.py

Input:
#.npz file of step (1)

#Output: "LCOE_WaveNC_2009_2013.png"

