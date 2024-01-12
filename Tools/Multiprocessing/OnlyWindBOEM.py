#env Gurobi
import numpy as np
from datetime import datetime, timedelta
import sys
from tqdm import tqdm
from multiprocessing import Pool

# # sys.path.append('../Tools')
sys.path.append('./Tools')
from  Port_Opt_MaxGeneration import SolvePortOpt_MaxGen_LCOE_Iterator


GeneralPathResources="./OutputData/"
PathWindDesigns=[]

PathWindDesigns.append(GeneralPathResources+"Wind/BOEM_Upscale24h_0.02Degree_GenCost_ATB_8MW_2020_Vestas.npz")
PathWindDesigns.append(GeneralPathResources+"Wind/BOEM_Upscale24h_0.02Degree_GenCost_ATB_12MW_2030.npz")
PathWindDesigns.append(GeneralPathResources+"Wind/BOEM_Upscale24h_0.02Degree_GenCost_ATB_15MW_2030.npz")
PathWindDesigns.append(GeneralPathResources+"Wind/BOEM_Upscale24h_0.02Degree_GenCost_ATB_18MW_2030.npz")

PathKiteDesigns=[]
PathWaveDesigns=[]


PathTransmissionDesign=[]
PathTransmissionDesign.append(GeneralPathResources+"Transmission/Transmission_1200MW.npz")
PathTransmissionDesign.append(GeneralPathResources+"Transmission/Transmission_1000MW.npz")
PathTransmissionDesign.append(GeneralPathResources+"Transmission/Transmission_600MW.npz")
PathTransmissionDesign.append(GeneralPathResources+"Transmission/Transmission_300MW.npz")




LCOE_RANGE=list(range(200,119,-5))+list(range(120,30,-2))
Max_CollectionRadious=30
MaxDesingsKite=1
MaxDesignsWind=1
MaxDesingsWave=0
MinNumWindTurb=0
MinNumWaveTurb=0
MinNumKiteTrub=0

def Iterator_WindBOEM(target):
    i=-1
    for PathTransmissionDesign_i in tqdm(PathTransmissionDesign):
        for PathWindDesigns_i in PathWindDesigns:
            TurbineCaseName=PathWindDesigns_i.rsplit(r"/")[-1][:-4]
            TransmissionCaseName=PathTransmissionDesign_i.rsplit(r"/")[-1][:-4]
            
            SavePath="./OutputData/Portfolios/Wind_"+TurbineCaseName+"_"+TransmissionCaseName+".npz"
            ReadMe="Case with wind on BOEM regions, considering a 1.2GW, 1.0, 0.6, 0.3 or 0.1GW transmission system, 30km radious and 1 design for each tech\
                \n Wind designs: 8MW Vestas 2020, 12MW 2030, 15MW 2030, 18MW 2030"
            
            i=i+1
            if i==target:
                #Create and solve the optimization problem
                SolvePortOpt_MaxGen_LCOE_Iterator([PathWindDesigns_i], PathWaveDesigns, PathKiteDesigns, PathTransmissionDesign_i, LCOE_RANGE\
                    ,Max_CollectionRadious,MaxDesignsWind, MaxDesingsWave, MaxDesingsKite,MinNumWindTurb,MinNumWaveTurb,MinNumKiteTrub\
                    ,ReadMe,SavePath=SavePath)


