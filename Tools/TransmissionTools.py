#env geo_env
import numpy as np
import os
import xarray as xr
import pandas as pd
from tqdm import tqdm
from GeneralGeoTools import GetDistanceToShore, GetDepth


#Tools to compute relevant transmission data

def TL_PowerParam_AC(AC_CableData, idxCable, D_SL, RatedPower_Generation, AddExtraCable=0):
    #Determined power parameters for a transmission line and admitance (AC Configuration)
    
    #Params:
    #D_SL: Distance from the shore to the wind farm [km]
    #AC_CableData: Excel table with AC cable data
    #idxCable: Index of the cable you are investigating
    #AddExtraCable: Force the addition of and extra cable to the transmission line
    
    LineLength=D_SL*1.2# Line length [km], considered 120% of shore distance
    V=AC_CableData['Voltage [KV]'][idxCable]
    R=AC_CableData['Resistance [Ohm/Km]'][idxCable]
    C=AC_CableData['Capacitance [MicroF/Km]'][idxCable]
    L=AC_CableData['Inductance [mF/Km]'][idxCable]
    Max_MVA=AC_CableData["MVA Capacity"][idxCable]
    
    Cost_PM=AC_CableData['Cable Cost [$/m]'][idxCable]
    
    Gamma=np.sqrt((R + (2*np.pi*60*L*10**(-3))*1.j)*(2*np.pi*60*C*10**(-6)*1.j)) #Propagation Constant [km^(-1)]
    Y=2*np.pi*60*C*10**(-6)*LineLength*np.tanh(Gamma*LineLength/2)/(Gamma*LineLength/2)*1.j #Admitance [1/Ω]
    Q=Y.imag*V**2 #Reactive power [MVar] for a single line 
      
    #P_LT: Max MW possible to transmit
    if (Max_MVA**2-(Q/2)**2)<=0:
        P_LT=-1
        S_LT=-1
        NumConductors=0 #the cable cannot be used to trasmit power for this system
        #It would only carry reactive power
        
    else:
        P_LT=np.sqrt(Max_MVA**2-(Q/2)**2) #Per conductor max [MW]
        if  RatedPower_Generation>P_LT:
            NumConductors=np.ceil(RatedPower_Generation/P_LT)
            
        if  RatedPower_Generation<=P_LT:
            NumConductors=1

        NumConductors=NumConductors+AddExtraCable
        
        S_LT=Max_MVA*NumConductors # Total power capacity of the line [MVA]
        P_LT=P_LT*NumConductors    # Total power capacity of the line [MW]
    
    Total_Q=Q*NumConductors #Total reactive power requirements [MVar]
    
    #P_LT: Max active power possible to transmit (Sum of all conductors) [MW]
    #S_LT: Max power capacity of the line (Sum of all conductors) [MVA]
    #Total_Q: Total reactive power requirements for the transmission line (Sum of all conductors) [MVar]
    #Y: Admitance of one line [1/Ω]
    
    #When P_LT=S_LT=-1 and NumConductors=0 the cable cannot be used to trasmit power for this system
    return P_LT, S_LT, Total_Q, NumConductors, Y
    
    
def TL_Efficiency_AC(AC_CableData, idxCable, NumConductors, D_SL, RatedPower_Generation, CF=0.5):
    #Determined the efficiency of a transmission line (AC Configuration)
    
    #Params:
    #AC_CableData: Excel table with AC cable data
    #idxCable: Index of the cable you are investigating, from AC_CableData
    #NumConductors: Number of conductors in the transmission line
    #D_SL: Distance from the shore to the wind farm [km]
    #RatedPower_Generation: Rated power of the wind farm [MW]  
    
    V=AC_CableData['Voltage [KV]'][idxCable]
    R=AC_CableData['Resistance [Ohm/Km]'][idxCable]
    C=AC_CableData['Capacitance [MicroF/Km]'][idxCable]
    L=AC_CableData['Inductance [mF/Km]'][idxCable]
    Max_MVA=AC_CableData["MVA Capacity"][idxCable]
    
    RatedPower_Generation=RatedPower_Generation*CF
    EffTransformer=EffTransformer=0.997 
    LT_Length=D_SL*1.2# Line length [km], considered 120% of shore distance

    TRL_OFF_AC=(1-EffTransformer)*RatedPower_Generation
    LC_AC=NumConductors*R*LT_Length*((RatedPower_Generation*EffTransformer/(V*NumConductors))**2)#Cable energy losses [MW]

    TRL_ON_AC=(RatedPower_Generation*EffTransformer-LC_AC)*(1-EffTransformer)# Onshore terminal losses [MW]
    
    LossesAC=TRL_OFF_AC+LC_AC+TRL_ON_AC
    EfficiencyAC=1-LossesAC/(RatedPower_Generation)
    
    #EfficiencyAC: Efficiency of the transmission line [EnergyOut/EnergyIn]
    #LossesAC: Total losses of the transmission line [MW]
    return EfficiencyAC, LossesAC


def TL_AnnualizedCost_AC(AC_CableData, idxCable, D_SL, RatedPower_Generation, AddExtraCable=0):
    #Determined the annualized cost of a transmission line (AC Configuration)
    
    #Params:
    #AC_CableData: Excel table with AC cable data 
    #idxCable: Index of the cable you are investigating, from AC_CableData
    #D_SL: Distance from the shore to the wind farm [km]
    #RatedPower_Generation: Rated power of the wind farm [MW]
    
    
    LineLength=D_SL*1.2# Line length [km], considered 120% of shore distance
    
    
    P_LT, S_LT, Total_Q, NumConductors, Y=TL_PowerParam_AC(AC_CableData, idxCable, D_SL, RatedPower_Generation, AddExtraCable=AddExtraCable)
    
    Cost_PC=AC_CableData['Cable Cost [$/m]'][idxCable]*10**-3 #Cost per conductor [M$/km]
    
    if NumConductors!=0: 
        
        EfficiencyAC, LossesAC=TL_Efficiency_AC(AC_CableData, idxCable, NumConductors, D_SL, RatedPower_Generation, CF=0.5)       
        
        FCR=0.113 #Factor of capital retrun
        OPPC_AC=6.55 + 0.0472*S_LT #Platform and plant cost HVAC [M$]
        OPC_AC=0.03434*S_LT**(0.7513)# Onshore plant cost [M$]
        QC_AC=0.0262*Total_Q# Cost of reactive power compensation [M$]

        CC_AC=Cost_PC*NumConductors*LineLength + 0.221*D_SL + 4.245*10**-3*S_LT + 0.629# Cable cost + cable landing + cable instalation [M$]
        CAPEX_AC=OPPC_AC+OPC_AC+QC_AC+CC_AC# Capital expenditures [M$]


        OPEX_AC=0.025*CAPEX_AC# Operational Expenditures [M$/year]

        AnnualizedCost_AC=FCR*CAPEX_AC + OPEX_AC# Anualized transmission costs [M$/year]
        
    else:#Specific circuit not feasible
        AnnualizedCost_AC=-1
        EfficiencyAC=-1
    
    MaxActivePower=P_LT #MW
    return AnnualizedCost_AC, NumConductors, EfficiencyAC, MaxActivePower

def TL_Efficiency_DC(DC_CableData, idxCable, NumConductors, D_SL, RatedPower_Generation, CF=0.5):
    #Compute Efficiency of the DC transmission line
    
    #Params:
    #DC_CableData: Excel table with DC cable data
    #idxCable: Index of the cable you are investigating, from DC_CableData
    #NumConductors: Number of conductors in the transmission line
    #D_SL: Distance from the shore to the wind farm [km]
    #RatedPower_Generation: Rated power of the wind farm [MW]  
    

    EffConverter=0.982 #Efficiency of the transformer
    LineLength=D_SL*1.2# Line length [km], considered 120% of shore distance
    S_LT=RatedPower_Generation*CF
    
    V=DC_CableData['Voltage [KV]'][idxCable]
    R=DC_CableData['Resistance [Ohm/Km]'][idxCable]
    Max_MVA=DC_CableData["MVA Capacity"][idxCable]
    CableCostDC=DC_CableData['Cable Cost [$/m]'][idxCable]
    
    TRL_OFF_DC=(1-EffConverter)*S_LT#Offshore terminal losses [MW]
    LC_DC=2*R*LineLength*((S_LT*EffConverter/(2*V*NumConductors))**2)#Cable energy losses [MW]
    TRL_ON_DC=(S_LT*EffConverter-LC_DC)*(1-EffConverter) # Onshore terminal losses [MW]
    LossesDC=TRL_OFF_DC+LC_DC+TRL_ON_DC
    EfficiencyDC=1-LossesDC/(S_LT)
    
    #EfficiencyAC: Efficiency of the transmission line [EnergyOut/EnergyIn]
    MaxActivePower=Max_MVA*NumConductors #Max input active power on the transmission system [MW]
    return EfficiencyDC, MaxActivePower

   
def TL_AnnualizedCost_DC(DC_CableData, idxCable, D_SL, RatedPower_Generation, AddExtraCable=0):
    #Determined the annualized cost of a transmission line (DC Configuration)
    
    #Params:
    #AC_CableData: Excel table with DC cable data 
    #idxCable: Index of the cable you are investigating, from DC_CableData
    #D_SL: Distance from the shore to the wind farm [km]
    #RatedPower_Generation: Rated power of the wind farm [MW]
    
    
    EffConverter=0.982 #Efficiency of the transformer
    LineLength=D_SL*1.2# Line length [km], considered 120% of shore distance
    FCR=0.113 #Factor of capital retrun
    
    S_LT=RatedPower_Generation
    V=DC_CableData['Voltage [KV]'][idxCable]
    R=DC_CableData['Resistance [Ohm/Km]'][idxCable]
    Max_MVA=DC_CableData["MVA Capacity"][idxCable]
    CableCostDC=DC_CableData['Cable Cost [$/m]'][idxCable]*10**-3
    
    NumConductors=np.ceil(RatedPower_Generation/Max_MVA)+AddExtraCable
    
    OPPC_DC=32.75 + 0.07205*S_LT #Platform and plant cost HVDC [M$]
    OPC_DC=0.1067*S_LT # Onshore plant cost [M$]
    
    
    CC_DC=CableCostDC*NumConductors*LineLength + 0.221*D_SL + 4.245*10**-3*S_LT + 0.629# Cable cost + cable landing + cable instalation [M$]
    CAPEX_DC=OPPC_DC+OPC_DC+CC_DC # Capital expenditures [M$]
    
    OPEX_DC=0.025*CAPEX_DC# Operational Expenditures [M$/year]
    
    AnnualizedCost_DC=FCR*CAPEX_DC + OPEX_DC# Anualized transmission costs [M$/year]
    
    EfficiencyDC, MaxActivePower=TL_Efficiency_DC(DC_CableData, idxCable, NumConductors, D_SL, RatedPower_Generation, CF=0.5)
    
    return AnnualizedCost_DC, NumConductors, EfficiencyDC, MaxActivePower



def GetBestTransmission(InputDataPath, AC_DataPath, DC_DataPath, RatedPower_Generation, LatMaxMin=(33.5, 37), LongMaxMin=(-78.5, -74.5), StepsPerDegree=10, SavePath=None, FixToCopper=False):
    #StepsPerDegree: Number of possible locations per degree of latitude and longitude
    #LatMaxMin: Maximum and Minimum Latitude
    #LongMaxMin: Maximum and Minimum Longitude
    
    
    AC_CableData = pd.read_excel(AC_DataPath)
    DC_CableData = pd.read_excel(DC_DataPath)
    
    if FixToCopper:
        #Fix to copper
        AC_CableData=AC_CableData.loc[AC_CableData["Observations"]=="Copper",:]
    
    #Create Grid of Latitudes and Longitudes for the center of energy colection of the transmission system   
    lat_min=LatMaxMin[0]-0.1
    lat_max=LatMaxMin[1]+0.1
    lon_min=LongMaxMin[0]-0.1
    lon_max=LongMaxMin[1]+0.1
    
    xlim =[lon_min,lon_max]
    ylim =[lat_min, lat_max]

    lat_range = np.linspace(lat_min, lat_max,int((lat_max-lat_min)*StepsPerDegree))
    lon_range = np.linspace(lon_min, lon_max, int((lon_max-lon_min)*StepsPerDegree))
            
    TL_LatLong = np.array(np.meshgrid(lat_range,lon_range)).T.reshape(-1, 2)
    

    #Get Minimum Distance to Shore for all points in the grid
    ShoreDistance=GetDistanceToShore(InputDataPath, TL_LatLong)
    Depth=GetDepth(InputDataPath, TL_LatLong)

    IdxFeasibleSites=(Depth>5)*(ShoreDistance>10)*(Depth<=2500)# Keep only sites with depth greater than 9 m and distance to shore greater than 2 km

    #Filter sites
    TL_LatLong=TL_LatLong[IdxFeasibleSites,:]
    TL_Depth=Depth[IdxFeasibleSites]
    TL_ShoreDistance=ShoreDistance[IdxFeasibleSites]
    
    #Get the best transmission line for each site
    S_BestCable=[]
    S_BestAnnualizedCost=[]
    S_Efficiency=[]
    S_Mode=[]
    S_NumConductors=[]
    S_MaxCableCapacity=[] #Maximum Active power for the transmission >= Design power

    for i in tqdm(range(len(TL_ShoreDistance))):
        D_SL=TL_ShoreDistance[i]
        
        MinEstimateCostWithLoss=np.inf
        BestCable=-1
        Best_Efficiency=-1
        MinCost=-1
        Mode=-1
        Tmp_NumConductors=-1
        Tmp_MaxCableCapacity=-1
        
        #AC Transmission
        #Run with minimum number of conductors
        for idxCable in range(len(AC_CableData)):
            
            AnnualizedCost_AC, NumConductors, EfficiencyAC, MaxCableCapacity = TL_AnnualizedCost_AC(AC_CableData, idxCable, D_SL, RatedPower_Generation)

            
            if AnnualizedCost_AC!=-1:
                
                #Find the cost of the transmission line taking into account the cost of enegy losses priced at 83$/MWh
                CostWithLoses=AnnualizedCost_AC*10**6 + 83*RatedPower_Generation*0.5*8760*(1-EfficiencyAC) #83$ from https://atb.nrel.gov/electricity/2022/offshore_wind
            
                if MinEstimateCostWithLoss>CostWithLoses :
                    MinEstimateCostWithLoss=CostWithLoses
                    
                    MinCost=AnnualizedCost_AC
                    BestCable=idxCable
                    Best_Efficiency=EfficiencyAC
                    Mode="HVAC"
                    Tmp_NumConductors=NumConductors
                    Tmp_MaxCableCapacity=MaxCableCapacity
                    
        #Run with minimum number of conductors + 1
        for idxCable in range(len(AC_CableData)):
            
            AnnualizedCost_AC, NumConductors, EfficiencyAC, MaxCableCapacity=TL_AnnualizedCost_AC(AC_CableData, idxCable, D_SL, RatedPower_Generation, AddExtraCable=1)
            
            if AnnualizedCost_AC!=-1:
                
                #Find the cost of the transmission line taking into account the cost of enegy losses priced at 83$/MWh
                CostWithLoses=AnnualizedCost_AC*10**6 + 83*RatedPower_Generation*0.5*8760*(1-EfficiencyAC) #83$ from https://atb.nrel.gov/electricity/2022/offshore_wind
            
                if MinEstimateCostWithLoss>CostWithLoses :
                    MinEstimateCostWithLoss=CostWithLoses
                    
                    MinCost=AnnualizedCost_AC
                    BestCable=idxCable
                    Best_Efficiency=EfficiencyAC
                    Mode="HVAC"
                    Tmp_NumConductors=NumConductors
                    Tmp_MaxCableCapacity=MaxCableCapacity
                    
                    
        #DC Transmission
        #Run with minimum number of conductors
        for idxCable in range(len(DC_CableData)):
            AnnualizedCost_DC, NumConductors, EfficiencyDC, MaxCableCapacity=TL_AnnualizedCost_DC(DC_CableData, idxCable, D_SL, RatedPower_Generation)
            

            if AnnualizedCost_DC!=-1:
                CostWithLoses=AnnualizedCost_DC*10**6 + 83*RatedPower_Generation*0.5*8760*(1-EfficiencyDC) #83$ from https://atb.nrel.gov/electricity/2022/offshore_wind
            
                if MinEstimateCostWithLoss>CostWithLoses :
                    MinEstimateCostWithLoss=CostWithLoses
                    
                    MinCost=AnnualizedCost_DC
                    BestCable=idxCable
                    Best_Efficiency=EfficiencyDC
                    Mode="HVDC"
                    Tmp_NumConductors=NumConductors
                    Tmp_MaxCableCapacity=MaxCableCapacity
        
        #Run with minimum number of conductors +1
        for idxCable in range(len(DC_CableData)):
            AnnualizedCost_DC, NumConductors, EfficiencyDC, MaxCableCapacity=TL_AnnualizedCost_DC(DC_CableData, idxCable, D_SL, RatedPower_Generation, AddExtraCable=1)
            

            if AnnualizedCost_DC!=-1:
                CostWithLoses=AnnualizedCost_DC*10**6 + 83*RatedPower_Generation*0.5*8760*(1-EfficiencyDC) #83$ from https://atb.nrel.gov/electricity/2022/offshore_wind
            
                if MinEstimateCostWithLoss>CostWithLoses :
                    MinEstimateCostWithLoss=CostWithLoses
                    
                    MinCost=AnnualizedCost_DC
                    BestCable=idxCable
                    Best_Efficiency=EfficiencyDC
                    Mode="HVDC"
                    Tmp_NumConductors=NumConductors
                    Tmp_MaxCableCapacity=MaxCableCapacity
        
        
        S_BestCable.append(BestCable)
        S_BestAnnualizedCost.append(MinCost)
        S_Efficiency.append(Best_Efficiency)
        S_Mode.append(Mode)
        S_NumConductors.append(Tmp_NumConductors)
        S_MaxCableCapacity.append(Tmp_MaxCableCapacity)
    

    S_BestACost=np.array(S_BestAnnualizedCost)
    S_BestCable=np.array(S_BestCable)
    S_Efficiency=np.array(S_Efficiency)
    S_Mode=np.array(S_Mode)
    S_NumConductors=np.array(S_NumConductors)
    S_MaxCableCapacity=np.array(S_MaxCableCapacity)
    
    # LCOE_SimpleApproximation=S_BestACost*10**6/(RatedPower_Generation*0.5*8760*S_Efficiency) # LCOE assuming a CF of 0.5 #[$/MWh*year]    
    LCOE_SimpleApproximation=(S_BestACost*10**6 + 83*RatedPower_Generation*0.5*8760*(1-S_Efficiency))/(RatedPower_Generation*0.5*8760) # LCOE assuming a CF of 0.5 #[$/MWh*year]    
     
    #Remove LatLong points that do not have a feasible transmission line
    #Maximum deph of 2500 m is also considered
    FilterIdx=(S_BestCable!=-1) * (TL_Depth<2500)
    
    S_Mode=S_Mode[FilterIdx]
    S_BestCable=S_BestCable[FilterIdx]
    S_BestACost=S_BestACost[FilterIdx]
    S_Efficiency=S_Efficiency[FilterIdx]
    S_NumConductors=S_NumConductors[FilterIdx]
    S_MaxCableCapacity=S_MaxCableCapacity[FilterIdx]
    LCOE_SimpleApproximation=LCOE_SimpleApproximation[FilterIdx]
    TL_LatLong=TL_LatLong[FilterIdx,:]
    TL_ShoreDistance=TL_ShoreDistance[FilterIdx]
    TL_Depth=TL_Depth[FilterIdx]
    
    
    TransmissionLineParameters={"S_Mode":S_Mode,
                                "S_BestCable":S_BestCable,
                                "S_BestACost":S_BestACost,
                                "S_Efficiency":S_Efficiency,
                                "S_NumConductors":S_NumConductors,
                                "LCOE_SimpleApproximation":LCOE_SimpleApproximation,
                                "TL_LatLong":TL_LatLong,
                                "TL_ShoreDistance":TL_ShoreDistance,
                                "TL_Depth":TL_Depth,
                                "DC_CableData": DC_CableData,
                                "AC_CableData": AC_CableData,
                                "RatedPowerMW": RatedPower_Generation , 
                                "MaxCableCapacity": S_MaxCableCapacity}#Considering only the cables, not the transformers (Transformers rated for RatedPower_Generation)
    
    
    if SavePath!=None:
        np.savez(SavePath, TransmissionLineParameters=TransmissionLineParameters)
        
    return TransmissionLineParameters

