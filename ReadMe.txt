This code contains functions to compute optimal offshore energy porfolios {Wind-ATB22 devices}, {Wave- }, and {Ocean Current- Kite Energy Converters} 



#######-------Main
1) Transmission.ipynb: Interactive code to determine transmission system cost/configuration for all locations of interest
2) DownloadWindDataNREL.ipynb: Donwload wind energy data from NREL
3) GetWindEnergyCostGeneration.ipynb: Get costs and energy generation for wind energy on a given region

#######-------./Tools/
1) TransmissionTools.py: Tools to compute transmission cost and circuit parameters
2) GeneralGeoToos.py: Toods to compute distance and depth for offshore portfolios, and perfrom general plots
3) NREL_hsdsConfiguration.py: Configure NREL API for wind energy download (register Key)
4) WindTurbineTools.py: Tools to compute wind energy and costs


#######-------./InputData/
1) CoastLine: Division for the NC coast and states
2) Transmission: Costs and parameters for AC and DC cables. On Refecens/ you will fid the references used to estimate these values
3) Wind: Turbine Parameters, Turbine costs, wind speed data downloaded from NREL toolkit

#######-------./OutputData/
1) Transmission/Transmission_...: Transmission annualized costs and efficiency for different rate transmission capacities
2) Wind/GenCost_...: Generation and Annualized Costs for different turbine designs:





#####to do things
Depth curves on BOEM plots?

