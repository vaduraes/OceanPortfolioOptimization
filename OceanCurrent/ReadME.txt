1)To compute current energy in each site location merging the data from HYCOM and MABSAB run:

Merge_Hycom_Mabsab.py

Inputs:
#TurbineFile: File containing information of the turbine model used
#MABSAB_Dictionary: File with MBSAB data
#HYCOM_Dictionary: File with HYCOM data
#GEO_data: directory with the NC contours 

Outputs:
#.npz file with CF and current speed in each feasible site location together with, LatLong, Depth,
minimum distance from shore for each site location, time variables, and rated power of the turbine model 
used



--------------------------------------------------------

2) to generate geo plots of the ocean current CF and speed run:

OceanCurrentGeoPlot.py

Inputs:
#HourlyData: .npz file with CF and ocean current Speed as well as LatLong for each site location
#GEO_data: directory with the NC contours

--------------------------------------------------------

3) To concatenate annualized cost information in the .npz file constructed in (1) run:

ConcatenateCost.py

Input:
#.npz file of step (1)

--------------------------------------------------------

4) To generate an geo plot of the LCOE of the ocean current resource in the NC region run:

LCOEGeoPlot.py

Input:
#.npz file of step (1)

#Output: "LCOE_CurrentNC_2009_2013.png"
