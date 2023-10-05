1)To convert the .mat files to a single .npz file with all current speeds run:

Create_NPZ_From_MAT.py

Inputs:
#YearsAnalysed> Array with the set of years you want to get the data
#FileSaveName> Name of the output saved files (.npz format)
#TurbineFile> .txt file with the turbine information (e.g. RM4Sandia.txt)
#CurrentRawData> Raw data with current speed (e.g. MABSAB_2014.mat)
#DepthDomain> Depth domain for the speed data (.mat)

Outputs:
#.npz file with ocean current speeds in each feasible site location in a specific year
