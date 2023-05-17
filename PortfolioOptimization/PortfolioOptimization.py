#Portfolio optimization all three energy sources
#When running case 50-100-100: Change mip gap to 0.02 and num iterations to 10000: Non convergence error other
from __future__ import division
from pyomo.environ import *
import numpy as np
import time

from GetIdxInRadious import GetIdxInRadious, GetIdxOutRadious#Get the index of all site locations inside a specific radious
from Upscale import Upscale, NewDataForStep2
from PrepareData import PrepareData #Assure consistency of the data in terms of time interval

start_time = time.time()
S_CPU_TIME=time.perf_counter()
######################

#LCOE range we want to investigate
#LCOE_Max=range(500,50,-10)
ResultsFileName='PortfolioOptimizationWindWaveOcean(50_0_50)'
LCOE_Max=range(500,50,-2)

Nt=[50,0,50]# Total number of installed turbine units Wind/Wave/Ocean

#Upper bound for the number of installed turbines in a single cell Wind/Wave/Ocean
NuWiWaOcean=[4,50,32]

Data=PrepareData()

WindEnergy=Data["WindEnergy"]*10**-6 #[MW]
WindLatLong=Data["WindLatLong"]
AnnualizedCostWind=Data["AnnualizedCostWind"]*10**6 #[$/Year]

NumWindSites=WindEnergy.shape[0]
RatedPowerWind=Data["RatedPowerWind"]

WaveEnergy=Data["WaveEnergy"]*10**-6 #[MW]
WaveLatLong=Data["WaveLatLong"]
NumWaveSites=WaveEnergy.shape[0]
RatedPowerWave=Data["RatedPowerWave"]
AnnualizedCostWave=Data["AnnualizedCostWave"]*10**6 #[$/Year-unit]


OceanEnergy=Data["OceanEnergy"]*10**-6 #[MW]
OceanLatLong=Data["OceanLatLong"]
AnnualizedCostOcean=Data["AnnualizedCostOcean"] #[$/Year]
#OceanEnergy=OceanEnergy[0:500,:]
#OceanLatLong=OceanLatLong[0:500,:]
#AnnualizedCostOcean=AnnualizedCostOcean[0:500]

NumOceanSites=OceanEnergy.shape[0]      
RatedPowerOcean=Data["RatedPowerOcean"]


RatedPowerSystem=(Nt[0]*RatedPowerWind)*10**(-6)+200 #300MW wind/ 200MW Kite

#Eliminate the wind data if the number of wind turbines is zero
if Nt[0]==0:
    WindEnergy=OceanEnergy[[],:]
    WindLatLong=OceanLatLong[[],:]
    NumWindSites=WindEnergy.shape[0]
    AnnualizedCostWind=np.array([])
    

#Eliminate the wave data if the number of wave turbines is zero
if Nt[1]==0:
    WaveEnergy=OceanEnergy[[],:]
    WaveLatLong=OceanLatLong[[],:]
    NumWaveSites=WaveEnergy.shape[0]
    AnnualizedCostWave=np.array([])

#Eliminate the ocean current data if the number of ocean current turbines is zero
if Nt[2]==0:
    OceanEnergy=OceanEnergy[[],:]
    OceanLatLong=OceanLatLong[[],:]
    NumOceanSites=OceanEnergy.shape[0]
    AnnualizedCostOcean=np.array([])

#Vectorize maximum number of turbines per site location per technology
Nu=np.concatenate((np.full((NumWindSites), NuWiWaOcean[0]),np.full((NumWaveSites), NuWiWaOcean[1]),np.full((NumOceanSites), NuWiWaOcean[2])))

#Vectorize total number of installed turbines per technology
RHS_S=np.concatenate((np.full((NumWindSites), Nt[0]*RatedPowerWind),np.full((NumWaveSites), Nt[1]),np.full((NumOceanSites), 200*10**6)))

#Vectorize total number of installed turbines per technology
if NumOceanSites==0:
    RC=np.concatenate((np.full((NumWindSites), RatedPowerWind),np.full((NumWaveSites), Nt[1]),np.full((NumOceanSites), 200*10**6)))
else:
    RC=np.concatenate((np.full((NumWindSites), RatedPowerWind),np.full((NumWaveSites), Nt[1]),RatedPowerOcean[:,0]))

#Vectorize annualized cost for each site location and per technology
AnnCost=np.concatenate((AnnualizedCostWind,AnnualizedCostWave,AnnualizedCostOcean)) #Annualized cost [$/Year]

#Vectorize energy generation in each site location and energy resource
EnergyGeneration=np.concatenate((WindEnergy,WaveEnergy,OceanEnergy),axis=0)
EnergyGeneration=np.average(EnergyGeneration,axis=1)*8766#[MWh/Year] (8766 is the number of hours in a year)

Sigma=np.cov(np.concatenate((WindEnergy,WaveEnergy,OceanEnergy),axis=0))#Variance covariance matrix [MWh**2]

#Step to guarantee that the matrix get SDP
for i in range(Sigma.shape[0]):
    Sigma[i,i]=Sigma[i,i]+10**-7•

NumSites=Sigma.shape[0]

#Upscale the data
_,GrpupsWind=Upscale(WindLatLong, Step=0.1)

_,GrpupsWave=Upscale(WaveLatLong,Translate_IDX=WindLatLong.shape[0], Step=0.06)

_,GrpupsOcean=Upscale(OceanLatLong, Translate_IDX = (WindLatLong.shape[0] + WaveLatLong.shape[0]), Step=0.1)

RHS_K=np.concatenate((np.full((len(GrpupsWind)), Nt[0]*RatedPowerWind),np.full((len(GrpupsWave)), Nt[1]),np.full((len(GrpupsOcean)), 200*10**6)))

#Concatenate all groups
K_Groups=GrpupsWind+GrpupsWave+GrpupsOcean
NumK_var=len(GrpupsWind)+len(GrpupsWave)+len(GrpupsOcean)

#Store set of relaxed solutions 
RelaxedSolutionsLCOE=[]
RelaxedSolutionsVar=[]
RelaxedLCOETarget=[]
OptimalXSolutionsWind=[]
OptimalXSolutionsWave=[]
OptimalXSolutionsOcean=[]

OptimalSSolution=[]
OptimalKSolution=[]

RelaxedSolutionsLatLongWind=[]
RelaxedSolutionsLatLongWave=[]
RelaxedSolutionsLatLongOcean=[]

#Store set of MINLP solutions
MINLPSolutionsLCOE=[]
MINLPSolutionsVar=[]

OptimalYSolutionsWind=[]
OptimalYSolutionsWave=[]
OptimalYSolutionsOcean=[]
MINLPSolutionsLatLongWind=[]
MINLPSolutionsLatLongWave=[]
MINLPSolutionsLatLongOcean=[]
FeasibilityStep1=[]
FeasibilityStep2=[]

#In order to do the optimization process we divided the model in two parts 

#----------------------------Relaxed optimization model  ----------------------------Start
RelaxedMINLP = ConcreteModel()

RelaxedMINLP.SiteWind = RangeSet(0,NumWindSites-1)# Set of all wind sites
RelaxedMINLP.SiteWave = RangeSet(NumWindSites,NumWindSites+NumWaveSites-1)# Set of all wave sites
RelaxedMINLP.SiteOcean = RangeSet(NumWindSites+NumWaveSites,NumSites-1)# Set of all ocean current sites

RelaxedMINLP.Site = RangeSet(0,NumSites-1)# Set of all sites

RelaxedMINLP.x = Var(RelaxedMINLP.Site, domain=NonNegativeReals)# x is the number of turbines in each site location
RelaxedMINLP.k =  Var(range(NumK_var), domain=Binary)#k is an auxiliaty variable to track the region where the deplyments happened

#Objective is to minimize the total variance in the energy generation
def objective_rule(RelaxedMINLP):   
    xtS=[sum(RelaxedMINLP.x[i]*Sigma[i,j] for i in RelaxedMINLP.Site) for j in RelaxedMINLP.Site]
    
    return summation(xtS,RelaxedMINLP.x)

RelaxedMINLP.OBJ = Objective(rule = objective_rule, sense=minimize)

def NumTurbinesCell_rule(RelaxedMINLP,i):
    return RelaxedMINLP.x[i]<=Nu[i]

RelaxedMINLP.Turbines_Cell = Constraint(RelaxedMINLP.Site, rule=NumTurbinesCell_rule)

#Maximum radious of deployment of a certain technology due to collection system constraints
Radious=20
IdxIn=GetIdxInRadious(Radious, WindLatLong, WaveLatLong, OceanLatLong)


def CountS_for_K(RelaxedMINLP,i):
    InTheBatch=np.sum(IdxIn[K_Groups[i],:], axis=0)>=1
    
    return sum(RelaxedMINLP.x[j]*RC[j] for j in np.where(InTheBatch)[0])>=RelaxedMINLP.k[i]*RHS_K[i]

RelaxedMINLP.Sum_S_KGroup= Constraint(range(len(K_Groups)), rule=CountS_for_K)


if NumWindSites!=0:
    RelaxedMINLP.TMaxTurbinesWind = Constraint(expr=sum(RelaxedMINLP.x[i] for i in RelaxedMINLP.SiteWind)==50)
    
    RelaxedMINLP.Sum_k_Wind= Constraint(expr=sum(RelaxedMINLP.k[i] for i in range(0,len(GrpupsWind),1))==1)

if NumWaveSites!=0:    
    RelaxedMINLP.TMaxTurbinesWave = Constraint(expr=sum(RelaxedMINLP.x[i] for i in RelaxedMINLP.SiteWave)==50)
    
    RelaxedMINLP.Sum_k_Wave= Constraint(expr=sum(RelaxedMINLP.k[i] for i in range(len(GrpupsWind),len(GrpupsWind)+len(GrpupsWave),1))==1)

if NumOceanSites!=0:       
     RelaxedMINLP.TMaxTurbinesOcean = Constraint(expr=sum(RelaxedMINLP.x[i]*RC[i] for i in RelaxedMINLP.SiteOcean)>=200*10**6)
     RelaxedMINLP.TMaxTurbinesOcean = Constraint(expr=sum(RelaxedMINLP.x[i]*RC[i] for i in RelaxedMINLP.SiteOcean)<=205*10**6)
     
     RelaxedMINLP.Sum_k_Ocean= Constraint(expr=sum(RelaxedMINLP.k[i] for i in range(len(GrpupsWind)+len(GrpupsWave), NumK_var, 1))==1)

opt = SolverFactory('gurobi')
# opt.options['max_iter'] = 5000
opt.options['mipgap'] = 0.01

LowestLCOE=999

for LCOE_Idx in range(len(LCOE_Max)):
    
    LCOETarget=LCOE_Max[LCOE_Idx]
    
    if LCOETarget<LowestLCOE:    
        Bypass=0
        #upperbound for the LCOE
        RelaxedMINLP.LCOE_Target = Constraint(expr=(sum(RelaxedMINLP.x[i]*AnnCost[i] for i in RelaxedMINLP.Site) - LCOETarget*sum(RelaxedMINLP.x[i]*EnergyGeneration[i] for i in RelaxedMINLP.Site)) <= 0)
       
        print("Running Relaxed Model With LCOE= %.2f" % LCOETarget)
        
        try:
            results=opt.solve(RelaxedMINLP, tee=True)
        except:
            Bypass=1
            RelaxedMINLP.del_component(RelaxedMINLP.LCOE_Target)  
    
        if Bypass==0:
            if (results.solver.status == SolverStatus.ok) and (results.solver.termination_condition == TerminationCondition.optimal):
                FeasibilityStep1.append(1)
                RelaxedLCOETarget.append(LCOETarget)
                
                Optimal_X=RelaxedMINLP.x.get_values()
                Optimal_k=RelaxedMINLP.k.get_values()
                
                Optimal_X=np.reshape(np.array([Optimal_X[i] for i in RelaxedMINLP.Site]),(1,NumSites))
                Optimal_k=np.reshape(np.array([Optimal_k[i] for i in range(NumK_var)]),(1,NumK_var))
                
                OptimalKSolution.append(Optimal_k)
                
                Variance=float(np.dot(np.dot(Optimal_X,Sigma),Optimal_X.T))/RatedPowerSystem**2#Optimal solution. Minimum variance
        
                RelaxedSolutionsVar.append(Variance)
                CurrentLCOE=(np.sum(Optimal_X*AnnCost)/np.sum(Optimal_X*EnergyGeneration))
                RelaxedSolutionsLCOE.append(np.sum(Optimal_X*AnnCost)/np.sum(Optimal_X*EnergyGeneration))
                
                if LowestLCOE>CurrentLCOE:
                    LowestLCOE=CurrentLCOE
                    
                if NumWindSites!=0:
                    TempXWind=Optimal_X[0, RelaxedMINLP.SiteWind]
                    IdxWind=np.reshape(np.array([TempXWind> 10**-2]),-1)
                    OptimalXSolutionsWind.append(TempXWind[IdxWind])
                    RelaxedSolutionsLatLongWind.append(WindLatLong[IdxWind,:])  
                    
                if NumWaveSites!=0:  
                    TempXWave=Optimal_X[0, RelaxedMINLP.SiteWave]
                    IdxWave=np.reshape(np.array([TempXWave> 10**-2]),-1)
                    OptimalXSolutionsWave.append(TempXWave[IdxWave])
                    RelaxedSolutionsLatLongWave.append(WaveLatLong[IdxWave,:])
        
                if NumOceanSites!=0:        
                    TempXOcean=Optimal_X[0, RelaxedMINLP.SiteOcean]
                    IdxOcean=np.reshape(np.array([TempXOcean> 10**-2]),-1)
                    OptimalXSolutionsOcean.append(TempXOcean[IdxOcean])
                    RelaxedSolutionsLatLongOcean.append(OceanLatLong[IdxOcean,:])
                
                #Delete constraint for its modification in the next step of the for loop
                RelaxedMINLP.del_component(RelaxedMINLP.LCOE_Target)
            

            else:# Something else is wrong
                FeasibilityStep1.append(0)
                RelaxedMINLP.del_component(RelaxedMINLP.LCOE_Target)    
                OptimalSSolution.append(None)
                OptimalKSolution.append(None)
                OptimalXSolutionsWind.append(None)
                RelaxedSolutionsLatLongWave.append(None)                   
                OptimalXSolutionsWave.append(None)
                RelaxedSolutionsLatLongWave.append(None)               
                OptimalXSolutionsOcean.append(None)
                RelaxedSolutionsLatLongOcean.append(None)                
                RelaxedSolutionsVar.append(None)
                RelaxedSolutionsLCOE.append(None)
                RelaxedLCOETarget.append(None)
                break
                

##----------------------------Step2----------------------------Start

print("\n\n\n***********************************************" )  
for LCOEscenario in range(len(RelaxedLCOETarget)):
    
    LCOETarget=RelaxedLCOETarget[LCOEscenario]
    
    if LCOETarget!=None:
        
        All_K=OptimalKSolution[LCOEscenario]      
        All_IdxRadious=IdxIn
        
        WindEnergy_s2, WindLatLong_s2, AnnualizedCostWind_s2, WaveEnergy_s2, WaveLatLong_s2, AnnualizedCostWave_s2,\
        OceanEnergy_s2, OceanLatLong_s2, AnnualizedCostOcean_s2, RatedPowerOcean_s2=NewDataForStep2(All_K, All_IdxRadious, GrpupsWind, GrpupsWave, GrpupsOcean,
                            WindEnergy, WindLatLong, AnnualizedCostWind,
                            WaveEnergy, WaveLatLong, AnnualizedCostWave,
                            OceanEnergy, OceanLatLong, AnnualizedCostOcean, RatedPowerOcean)
        
        AnnCost_s2=np.concatenate((AnnualizedCostWind_s2,AnnualizedCostWave_s2,AnnualizedCostOcean_s2)) #Annualized cost [$/Year]
        
        EnergyGeneration_s2=np.concatenate((WindEnergy_s2,WaveEnergy_s2,OceanEnergy_s2),axis=0)
        EnergyGeneration_s2=np.average(EnergyGeneration_s2,axis=1)*8766#[MWh/Year]
        
        Sigma=np.cov(np.concatenate((WindEnergy_s2,WaveEnergy_s2,OceanEnergy_s2),axis=0)) #Variance covariance matrix [MWh**2]
        
        #Step to guarantee that the matrix get SDP
        for i in range(Sigma.shape[0]):
            Sigma[i,i]=Sigma[i,i]+10**-8

        
        print("Running MINLP With LCOE= %.2f" % LCOETarget)  
        
        MINLP = ConcreteModel()

        NumWindSites=WindEnergy_s2.shape[0]
        NumWaveSites=WaveEnergy_s2.shape[0]
        NumOceanSites=OceanEnergy_s2.shape[0] 
        NumSites=NumWindSites+NumWaveSites+NumOceanSites
        
        Nu=np.concatenate((np.full((NumWindSites), NuWiWaOcean[0]),np.full((NumWaveSites), NuWiWaOcean[1]),np.full((NumOceanSites), NuWiWaOcean[2])))

        #Vectorize total number of installed turbines per technology
        RHS_S=np.concatenate((np.full((NumWindSites), Nt[0]*RatedPowerWind),np.full((NumWaveSites), Nt[1]),np.full((NumOceanSites), 200*10**6)))


        #Vectorize total number of installed turbines per technology
        if NumOceanSites==0:
            RC=np.concatenate((np.full((NumWindSites), RatedPowerWind),np.full((NumWaveSites), Nt[1]),np.full((NumOceanSites), 200*10**6)))
        else:
            RC=np.concatenate((np.full((NumWindSites), RatedPowerWind),np.full((NumWaveSites), Nt[1]),RatedPowerOcean_s2[:,0]))


        MINLP.SiteWind = RangeSet(0,NumWindSites-1)
        MINLP.SiteWave = RangeSet(NumWindSites,NumWindSites+NumWaveSites-1)
        MINLP.SiteOcean = RangeSet(NumWindSites+NumWaveSites,NumSites-1)
        MINLP.Site = RangeSet(0,NumSites-1)

        MINLP.y = Var(MINLP.Site, domain=NonNegativeIntegers)# Integer variable to track the number of turbines used per site location

        MINLP.s = Var(MINLP.Site, domain=Binary)# Integer variable to track the site location used for each technology

        def objective_rule(MINLP):
            ytS=[sum(MINLP.y[i]*Sigma[i,j] for i in MINLP.Site) for j in MINLP.Site]
            return summation(ytS,MINLP.y)

        MINLP.OBJ = Objective(rule = objective_rule, sense=minimize)

        MINLP.LCOE_Target = Constraint(expr=(sum(MINLP.y[i]*AnnCost_s2[i] for i in MINLP.Site)-LCOETarget*sum(MINLP.y[i]*EnergyGeneration_s2[i] for i in MINLP.Site)) <= 0)

        def NumTurbinesCell_rule(MINLP,i):
            return MINLP.y[i]<=Nu[i]
    
        MINLP.Turbines_Cell = Constraint(MINLP.Site, rule=NumTurbinesCell_rule)


        if NumWindSites!=0:
            MINLP.TMaxTurbinesWind = Constraint(expr=sum(MINLP.y[i] for i in MINLP.SiteWind)==50)
            MINLP.ChooseOneCircleWind1= Constraint(expr=sum(MINLP.s[i] for i in MINLP.SiteWind)==1)
            
        if NumWaveSites!=0:
            MINLP.TMaxTurbinesWave = Constraint(expr=sum(MINLP.y[i] for i in MINLP.SiteWave)==Nt[1])
            MINLP.ChooseOneCircleWave1= Constraint(expr=sum(MINLP.s[i] for i in MINLP.SiteWave)==1)
        
        if NumOceanSites!=0:
            MINLP.TMaxTurbinesOcean = Constraint(expr=sum(MINLP.y[i]*RC[i] for i in MINLP.SiteOcean)>=200*10**6)
            MINLP.TMaxTurbinesOcean = Constraint(expr=sum(MINLP.y[i]*RC[i] for i in MINLP.SiteOcean)<=205*10**6)
            
            MINLP.ChooseOneCircleOcean1= Constraint(expr=sum(MINLP.s[i] for i in MINLP.SiteOcean)==1)
            
        Radious=20
        IdxIn_s2=GetIdxInRadious(Radious, WindLatLong_s2, WaveLatLong_s2, OceanLatLong_s2)
        
        def MaximumRadious(MINLP,i):  
            
            return sum(MINLP.y[j]*RC[j]  for j in np.where(IdxIn_s2[i,:])[0])>=MINLP.s[i]*RHS_S[i]
        
        MINLP.Maximum_Radious = Constraint(MINLP.Site, rule=MaximumRadious)
        
       
 
        opt = SolverFactory('gurobi')
        opt.options['mipgap'] = 0.005
        
        
        Bypass=0
        try:
            resultsMINLP=opt.solve(MINLP, tee=True)            
        
        except:
            Bypass=1
            FeasibilityStep2.append(0)
            MINLPSolutionsVar.append(None)
            MINLPSolutionsLCOE.append(None)
            OptimalYSolutionsWind.append(None)
            MINLPSolutionsLatLongWind.append(None)                   
            OptimalYSolutionsWave.append(None)
            MINLPSolutionsLatLongWave.append(None)               
            OptimalYSolutionsOcean.append(None)
            MINLPSolutionsLatLongOcean.append(None)           
            
        if Bypass==0:
            #solution in optimal and feasible
            if (resultsMINLP.solver.status == SolverStatus.ok) and (resultsMINLP.solver.termination_condition == TerminationCondition.optimal):
                FeasibilityStep2.append(1)
                Optimal_Y=MINLP.y.get_values()
                
                Optimal_Y=np.reshape(np.array([Optimal_Y[i] for i in MINLP.Site]),-1)
                
                VarianceMINLP=float(np.dot(np.dot(Optimal_Y,Sigma),Optimal_Y.T))#Optimal solution. Minimum variance
                
                VarianceMINLP=VarianceMINLP/RatedPowerSystem**2
                MINLPSolutionsVar.append(VarianceMINLP)
                
                MINLPSolutionsLCOE.append(np.sum(Optimal_Y*AnnCost_s2)/np.sum(Optimal_Y*EnergyGeneration_s2))#CF in the optimal solution.
                  
                
                if NumWindSites!=0:
                    TempYWind=Optimal_Y[ MINLP.SiteWind]
                    IdxWind=np.reshape(np.array([TempYWind> 10**-2]),-1)
                    OptimalYSolutionsWind.append(TempYWind[IdxWind])
                    MINLPSolutionsLatLongWind.append(WindLatLong_s2[IdxWind,:])
                    
                    
                if NumWaveSites!=0:  
                    TempYWave=Optimal_Y[MINLP.SiteWave]
                    IdxWave=np.reshape(np.array([TempYWave> 10**-2]),-1)
                    OptimalYSolutionsWave.append(TempYWave[IdxWave])
                    MINLPSolutionsLatLongWave.append(WaveLatLong_s2[IdxWave,:])
        
                if NumOceanSites!=0:        
                    TempYOcean=Optimal_Y[MINLP.SiteOcean]
                    IdxOcean=np.reshape(np.array([TempYOcean> 10**-2]),-1)
                    OptimalYSolutionsOcean.append(TempYOcean[IdxOcean])
                    MINLPSolutionsLatLongOcean.append(OceanLatLong_s2[IdxOcean,:])
               
            else:# Something else is wrong
                FeasibilityStep2.append(0)
                MINLPSolutionsVar.append(None)
                MINLPSolutionsLCOE.append(None)
                OptimalYSolutionsWind.append(None)
                MINLPSolutionsLatLongWind.append(None)                   
                OptimalYSolutionsWave.append(None)
                MINLPSolutionsLatLongWave.append(None)               
                OptimalYSolutionsOcean.append(None)
                MINLPSolutionsLatLongOcean.append(None)   
                
    else:
        FeasibilityStep2.append(0)
        MINLPSolutionsVar.append(None)
        MINLPSolutionsLCOE.append(None)
        OptimalYSolutionsWind.append(None)
        MINLPSolutionsLatLongWind.append(None)                   
        OptimalYSolutionsWave.append(None)
        MINLPSolutionsLatLongWave.append(None)               
        OptimalYSolutionsOcean.append(None)
        MINLPSolutionsLatLongOcean.append(None)   
  
elapsed_time = time.time() - start_time
E_CPU_TIME = time.perf_counter()-S_CPU_TIME
      
        
np.savez(ResultsFileName + ".npz", TotalNumTurbines=Nt,MaxNumTurbines=NuWiWaOcean,FeasibilityStep1=FeasibilityStep1,
         RelaxedLCOETarget=RelaxedLCOETarget, RelaxedSolutionsLCOE=RelaxedSolutionsLCOE, RelaxedSolutionsVar=RelaxedSolutionsVar,
         OptimalXSolutionsWind=OptimalXSolutionsWind, OptimalXSolutionsWave=OptimalXSolutionsWave,
         OptimalXSolutionsOcean=OptimalXSolutionsOcean, OptimalSSolution=OptimalSSolution,OptimalKSolution=OptimalKSolution,
         RelaxedSolutionsLatLongWind=RelaxedSolutionsLatLongWind,RelaxedSolutionsLatLongWave=RelaxedSolutionsLatLongWave,
         RelaxedSolutionsLatLongOcean=RelaxedSolutionsLatLongOcean,FeasibilityStep2=FeasibilityStep2,
         MINLPSolutionsLCOE=MINLPSolutionsLCOE, MINLPSolutionsVar=MINLPSolutionsVar,OptimalYSolutionsWind=OptimalYSolutionsWind,
         MINLPSolutionsLatLongWind=MINLPSolutionsLatLongWind,
         OptimalYSolutionsWave=OptimalYSolutionsWave, MINLPSolutionsLatLongWave=MINLPSolutionsLatLongWave, 
         OptimalYSolutionsOcean=OptimalYSolutionsOcean, MINLPSolutionsLatLongOcean=MINLPSolutionsLatLongOcean,
         elapsed_time=elapsed_time, E_CPU_TIME=E_CPU_TIME)
