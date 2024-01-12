This code contains functions to compute optimal offshore energy porfolios {Wind-ATB22 devices}, {Wave- }, and {Ocean Current- Kite Energy Converters} 


#######-------Main
1) Transmission.ipynb: Interactive code to determine transmission system cost/configuration for all locations of interest
2) DownloadWindDataNREL.ipynb: Donwload wind energy data from NREL
3) GetWindEnergyCostGeneration.ipynb: Get costs and energy generation for wind energy on a given region
4) PortfolioOptimization.ipynb: Solve portfolio optimization problems using multiprocessing
5) PlotPortfolioResults.ipynb: Generate plots from the portfolio results
6) DownloadWaveDataWWIII.ipynb: Code to download Hs and Tp (wave data) from wave watch III model
7) Wave_Grib2NPZ.ipynb: Convert GRIB2 to .npz file for the Wave Watch 3 data

#######-------./Tools/
1) TransmissionTools.py: Tools to compute transmission cost and circuit parameters
2) GeneralGeoToos.py: Toods to compute distance and depth for offshore portfolios, and perfrom general plots
3) NREL_hsdsConfiguration.py: Configure NREL API for wind energy download (register Key)
4) WindTurbineTools.py: Tools to compute wind energy and costs
5) Port_Opt_MaxGeneration: Model to solve the portfolio optimization for the max generation constrained to LCOE (can run wind, wave and kite)
6) DownloadWWIII_Wave.py: Code to download Hs and Tp (wave data) from wave watch III model
7) WWaveGrib2Npz.py: Convert GRIB2 to .npz file for the Wave Watch 3 data

#######-------./InputData/
1) CoastLine: Division for the NC coast and states
2) Transmission: Costs and parameters for AC and DC cables. On Refecens/ you will fid the references used to estimate these values
3) Wind: Turbine Parameters, Turbine costs, wind speed data downloaded from NREL toolkit
4) Wave: Wave Turbine data and Hs and Tp data from WWIIII

#######-------./OutputData/
1) Transmission/Transmission_...: Transmission annualized costs and efficiency for different rate transmission capacities
2) Wind/GenCost_...: Generation and Annualized Costs for different turbine designs:





#####to do things
Depth curves on BOEM plots?

