#Compute energy generated based on turbine type and wind speed data
import csv
import numpy as np

WindTurbine="Haliade150_6MW" #Name of the wind turbine file

Range_of_WindSpeeds=np.array([40,60,80,100,120,140,160,200])#m³/s
GridScale=2#km


#Read wind turbine data
def Condition_on_TurbineData(WindTurbine, Range_of_WindSpeeds):
    
    PowerCurve=[]
    PowerCurve.append(['0','0'])
    
    File_Wind=open(WindTurbine+'.txt', "r")
    File_Wind_csv=csv.reader(File_Wind,delimiter=';')
    
    for EachLine in File_Wind_csv:
        
        if File_Wind_csv.line_num==1:
            HubHeight=float(EachLine[1])
            
        if File_Wind_csv.line_num==2:
            Total_Efficiency=1-float(EachLine[1])
            
        if File_Wind_csv.line_num==3:
            RotorDiameter=float(EachLine[1])

        if File_Wind_csv.line_num==4:
            RatedPower=float(EachLine[1])
            
        if File_Wind_csv.line_num>=6:
            PowerCurve.append(EachLine[:])

    PowerCurve=np.asarray(PowerCurve, dtype='float32')
    
    
    #Get hub height as the closest value between the wind speed data we have and the value in the data sheet
    Index_Hub_WindSpeed=len(Range_of_WindSpeeds)-np.argmin(np.flip(np.abs(Range_of_WindSpeeds-HubHeight)))-1
    AdjustedHubHeight=Range_of_WindSpeeds[Index_Hub_WindSpeed]
    
    File_Wind.close()
    
    return AdjustedHubHeight, Total_Efficiency, RotorDiameter, PowerCurve, RatedPower


#Convert wind speed to energy in pu
def WindToEnergy(AdjustedHubHeight, PowerCurve, WindTurbine, Total_Efficiency):
       
    Modified_PowerCurve=np.zeros((80),dtype=float)
    Modified_PowerCurve[0:PowerCurve.shape[0]]=PowerCurve[:,1]
    
    WindSpeedFile="WindSpeedNREL_"+str(AdjustedHubHeight)+"m.npz"
    WindData=np.load(WindSpeedFile) 
    
    WindSpeed=WindData["windspeed"]
    
    WindSpeed=np.round(WindSpeed)
    
    IntValuesWindSpeed=WindSpeed.astype(int)
    
    WindEnergy=Modified_PowerCurve[IntValuesWindSpeed]*Total_Efficiency
    WindEnergy=WindEnergy.astype('float32')
    

    ReadMe='\
    EnergyPu: pu wind energy (Base 6MW)-1 unit per site\n\
    1) The data is in hourly discretization starting at 1/1/2007 and going up to\
     12/31/2013 23:00'
    
    np.savez("WindEnergyNREL_"+str(AdjustedHubHeight)+"m_"+WindTurbine+'.npz',ReadMe=ReadMe,WindEnergy=WindEnergy, RatedPower=RatedPower, LatLong=WindData["LatLong"],Depth=WindData["Depth"], DistanceToShore=WindData["DistanceToShore"])
    

AdjustedHubHeight, Total_Efficiency, RotorDiameter, PowerCurve,RatedPower=Condition_on_TurbineData(WindTurbine, Range_of_WindSpeeds)
WindToEnergy(AdjustedHubHeight, PowerCurve, WindTurbine, Total_Efficiency)


#Using the wind speed data from NREL each point can be considered as a 2km x 2km grid cell
#With this information and the hub height it is possible to determine the maximum number of units per grid cell