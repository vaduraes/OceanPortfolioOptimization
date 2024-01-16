#Run on geo_env
import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import sys
from scipy.stats import weibull_min
from scipy.stats import lognorm
from scipy.stats import norm
from tqdm import tqdm
from scipy.spatial import ConvexHull
import h5py

from scipy.optimize import curve_fit
from multiprocessing import Pool

# Here I am basing the theory on https://asmedigitalcollection.asme.org/offshoremechanics/article/137/3/031901/377037/Joint-Distribution-of-Environmental-Condition-at
#Rosenblatt transformation

def weibull_CDF(x, shape, loc, scale):
    return weibull_min.cdf(x, shape, loc=loc, scale=scale)

def ComputeExtremeSurface_Multiprocessing(s, Recurrency=50):
    Data=np.load("./InputData/Wave/EVA_Params.npz", allow_pickle=True)["Data"].item()
    ParamWS_Weibull=Data["ParamWS_Weibull"]
    ParamHs_Ws=Data["ParamHs_Ws"]
    HsUw_RefRange=Data["HsUw_RefRange"]
    ParamTpHs=Data["ParamTpHs"]
    Tp_Hs_RefRange=Data["Tp_Hs_RefRange"]
    
    #Definition on P(Ws)
    Param_Ws=ParamWS_Weibull

    #Definition on P(Hs|Ws)
    Param_Hs_Ws=ParamHs_Ws
    Bins_Hs_Ws=HsUw_RefRange

    #Definition on P(Tp|Hs)
    Param_Tp_Hs=ParamTpHs
    Bins_Tp_Hs=Tp_Hs_RefRange

    NSites=len(Param_Hs_Ws)


    #range of values we will numerically investigate
    Ws_EVA_Range=np.arange(0,40,0.2)
    Hp_EVA_Range=np.arange(0,30,0.2)
    Tp_EVA_Range=np.arange(0,60,0.5)

    #Definition on P(Ws)
    Param_Ws_s=Param_Ws[s]

    #Definition on P(Hs|Ws)
    Param_Hs_Ws_s=ParamHs_Ws[s]
    Bins_Hs_Ws_s=np.array(Bins_Hs_Ws[s])

    #Definition on P(Tp|Hs)
    Param_Tp_Hs_s=Param_Tp_Hs[s]
    Bins_Tp_Hs_s=np.array(Bins_Tp_Hs[s])

    ProbSurfaceVariables=[]
    ProbSurfaceValue=[]
    MaxRadius=norm.ppf(1-1/(365*24/3*Recurrency))
    MaxRadius2=MaxRadius**2
    for Ws_r in Ws_EVA_Range:
        
        #Wind
        y1=norm.ppf(weibull_CDF(Ws_r, Param_Ws_s[0], Param_Ws_s[1], Param_Ws_s[2]))
        if y1**2<=MaxRadius2:
            for Hp_r in Hp_EVA_Range:
                
                #Wave height
                #get index of the bin
                Idx_Hs_Ws_closest=np.max([np.argmax(Bins_Hs_Ws_s>=Ws_r)-1,0])
                Param_Hs_Ws_s_i=Param_Hs_Ws_s[Idx_Hs_Ws_closest]

                y2=norm.ppf(weibull_CDF(Hp_r, Param_Hs_Ws_s_i[0], Param_Hs_Ws_s_i[1], Param_Hs_Ws_s_i[2]))

                if y1**2+y2**2<=MaxRadius2:
                    for Tp_r in Tp_EVA_Range:
                        #Wave period
                        #get index of the bin
                        Idx_Tp_Hs_closest=np.max([np.argmax(Bins_Tp_Hs_s>=Hp_r)-1,0])
                        Param_Tp_Hs_s_i=Param_Tp_Hs_s[Idx_Tp_Hs_closest]

                        y3=norm.ppf(weibull_CDF(Tp_r, Param_Tp_Hs_s_i[0], Param_Tp_Hs_s_i[1], Param_Tp_Hs_s_i[2]))

                        r2=y1**2+y2**2+y3**2
                        if r2<=MaxRadius2:
                            ProbSurfaceVariables.append([Ws_r, Hp_r, Tp_r])
    
    return ProbSurfaceVariables

