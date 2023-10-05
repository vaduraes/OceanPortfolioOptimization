#Read turbine data for ocean current
import csv

#Get Turbine Data
def ReadTurbineData(TurbineFile):
    
    File_Turbine=open(TurbineFile, "r")
    File_Turbine_csv=csv.reader(File_Turbine,delimiter=',')
    
    for EachLine in File_Turbine_csv:
    
        if File_Turbine_csv.line_num==4:
            RotorDepth=float(EachLine[1])
            
        if File_Turbine_csv.line_num==5:
            NumTurbinesPerSet=float(EachLine[1])
            
        if File_Turbine_csv.line_num==6:
            RatePower=float(EachLine[1])
    
        if File_Turbine_csv.line_num==7:
            RotorDiamater=float(EachLine[1])
            
        if File_Turbine_csv.line_num==8:
            C_Performance=float(EachLine[1])#Coeficient of performance
            
        if File_Turbine_csv.line_num==9:
            E_Gearbox=float(EachLine[1])#Gearbox Efficiency
            
        if File_Turbine_csv.line_num==10:
            E_Generator=float(EachLine[1])#Generator Efficiency
            
        if File_Turbine_csv.line_num==11:
            E_Transformer=float(EachLine[1])
                    
        if File_Turbine_csv.line_num==12:
            E_Inverter=float(EachLine[1])
            
        if File_Turbine_csv.line_num==13:
            Combined_E_Turbine=float(EachLine[1])#Combined power chain conversion efficiency
            
        if File_Turbine_csv.line_num==14:
            Availability=float(EachLine[1])
            
        if File_Turbine_csv.line_num==15:#Efficiency of the collection system
            E_Collection=float(EachLine[1])
            
        if File_Turbine_csv.line_num==16:
            MinDepth=float(EachLine[1])

        if File_Turbine_csv.line_num==17:
            MaxDepth=float(EachLine[1])
            
        if File_Turbine_csv.line_num==18:
            Water_rho=float(EachLine[1])
            
            
    File_Turbine.close()
            
    Turbine={"RotorDepth":RotorDepth,
             "NumTurbinesPerSet":NumTurbinesPerSet,
             "RatePower":RatePower,
             "RotorDiamater":RotorDiamater,
             "C_Performance":C_Performance,
             "E_Gearbox":E_Gearbox,
             "E_Generator":E_Generator,
             "E_Transformer":E_Transformer,
             "E_Inverter":E_Inverter,
             "Combined_E_Turbine":Combined_E_Turbine,
             "Availability": Availability,
             "E_Collection": E_Collection,
             "MinDepth": MinDepth,
             "MaxDepth": MaxDepth,
             "Water_rho": Water_rho
              }    
    return Turbine