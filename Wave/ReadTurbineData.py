#Function to read the wave turbine data

import numpy as np
import csv

def ReadTurbineData(TurbineFile):
    
    File_Turbine=open(TurbineFile, "r")
    File_Turbine_csv=csv.reader(File_Turbine,delimiter=',')
    MP_Matrix=[]    
    
    for EachLine in File_Turbine_csv:
    
        
        if File_Turbine_csv.line_num==2:
            RatedPower=float(EachLine[1])
            
        if File_Turbine_csv.line_num==3:
            E_Mec2El=float(EachLine[1])#mechanical to electrical energy conversion efficiency
    
        if File_Turbine_csv.line_num==4:
            E_Av=float(EachLine[1])# losses due to device availability
            
        if File_Turbine_csv.line_num==5:
            E_Tr=float(EachLine[1])# transmission efficiency
 
        if File_Turbine_csv.line_num==6:
            MinDepth=float(EachLine[1])

        if File_Turbine_csv.line_num==7:
            MaxDepth=float(EachLine[1])
            
        
        if File_Turbine_csv.line_num==10:
            Te_Bins=EachLine[1:]
            Te_Bins=np.array(Te_Bins,dtype=float)
            Te_Bins=Te_Bins.reshape(len(Te_Bins),1)
            
        #Mechanical power matrix
        if File_Turbine_csv.line_num>=11:
            MP_Matrix.append(np.array(EachLine,dtype=float))
    
    File_Turbine.close()
            
    MP_Matrix=np.array(MP_Matrix,dtype=float)
    Hs_Bins=MP_Matrix[:,0]
    Hs_Bins=Hs_Bins.reshape(len(Hs_Bins),1)
    
    MP_Matrix=MP_Matrix[:,1:]

    Turbine={"RatedPower":RatedPower,
         "E_Mec2El":E_Mec2El,
         "E_Av":E_Av,
         "E_Tr":E_Tr,
         "MinDepth":MinDepth,
         "MaxDepth":MaxDepth,         
         "Te_Bins":Te_Bins,
         "Hs_Bins":Hs_Bins,
         "MP_Matrix":MP_Matrix,
          }    


    return Turbine