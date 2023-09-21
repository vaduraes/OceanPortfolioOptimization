1) To download NREL wind data (speed m3/s) run:

WindData_hsdsConfiguration.py
DownloadWindDataNREL.py

Inputs: 
#LatNC> Minimum and maximum values for Latitude 
#LonNC> Minimum and maximum values for Longitude
#Data_Name: Name of the file you want to download from the NREL database (available 10,40,60,80,100,120,140,160,200m. Wind speed and direction) 

Outputs:
#Wind speed and Lat Long for each wind speed measurements from 2007 to 2013 in hourly discretization (saved in .npz)

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

2)To reduce the size of the wind speed data run:

GeneralWindSpeed_To_NREL_Sites.py

Inputs:
#HubHeight> Hub height for the turbine studied
#NREL_File> NREL file with the "feasible" site locations for wind turbines in the NC
#File_CoastLine> .csv file describing the coastline of the NC
#Depth_NETCDF> .nc file describing the ocean depth in the NC region

Outputs:
#Wind speed and Lat Long for each wind speed measurements from 2007 to 2013 in hourly discretization just for the sites in the "NREL_File"
#Minimum distance from coastline and depth for each site location

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

3) To get wind energy run:

ComputeWindEnergy.py

Inputs:
#WindTurbine: Wind turbine file. A txt file with the information needed for the wind turbine
#WindSpeedFile: Adequate winds speed file for the turbine in hands (in terms of hub height)

Outputs:
#.npz file with energy information (in pu units) for each feasible wind site location in the NC region. The information is in daily discretization

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

4) To generate a Geo plot for the capacity factor of wind resource run:

WindCFGeoPlot.py

Inputs:
#ShapeFileCoast: Shapefile with information of the USA coast
#ShapeFileStates: Shapefile with information of the USA states
#HourlyData: Energy data for wind 

Outputs:
#.png geo plot showing the CF for the wind resource



-------------Comments/References
Limit sites based on NREL analysis: https://www.nrel.gov/gis/wind-supply-curves.html
Limit sites based on NREL: https://data.nrel.gov/submissions/54