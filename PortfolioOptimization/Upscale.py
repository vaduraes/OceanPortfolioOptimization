import numpy as np

#Upscale the grid but track which old site locations are in the new upscaled resolution
def Upscale(LatLong, Translate_IDX=0, Step=0.2):
    SiteGroups=[]
    SizeEachGroup=0
    if np.shape(LatLong)[0]!=0:
        
        #NewGridResolution, the smaller more time it takes to simulate
        #Step=0.1# In degrees, about 11km
        #Step=0.2# In degrees, about 20km
         
        LatMax=np.max(LatLong[:,0])+1
        LatMin=np.min(LatLong[:,0])-1
        
        LongMax=np.max(LatLong[:,1])+1
        LongMin=np.min(LatLong[:,1])-1
        
        
        for IdxLat in range(int((LatMax-LatMin)/Step)):
            for IdxLong in range(int((LongMax-LongMin)/Step)):
                L_Lat=LatMin+Step*IdxLat
                U_Lat=LatMin+Step*(IdxLat+1)
                
                L_Long=LongMin+Step*IdxLong
                U_Long=LongMin+Step*(IdxLong+1)
                
                LatTrueFalse=(LatLong[:,0]>=L_Lat)*(LatLong[:,0]<U_Lat)      
                LongTrueFalse=(LatLong[:,1]>=L_Long)*(LatLong[:,1]<U_Long)
                InOut=LatTrueFalse*LongTrueFalse
                if np.sum(InOut)!=0:
                    SiteGroups.append(np.where(InOut)[0]+Translate_IDX)
        
        SizeEachGroup=[]
        for i in range(len(SiteGroups)):   
            SizeEachGroup.append(len(SiteGroups[i]))
            
        SizeEachGroup=np.array(SizeEachGroup)
    
    return  SizeEachGroup, SiteGroups

#Define the new data we will use in step 2. More compact data considering only the circles
def NewDataForStep2(All_K, All_IdxRadious, GrpupsWind, GrpupsWave, GrpupsOcean,
                    WindEnergy, WindLatLong, AnnualizedCostWind,
                    WaveEnergy, WaveLatLong, AnnualizedCostWave,
                    OceanEnergy, OceanLatLong, AnnualizedCostOcean, RatedPowerOcean):
    
    NumWindSites=np.shape(WindEnergy)[0]
    NumWaveSites=np.shape(WaveEnergy)[0]
    NumOceanSites=np.shape(OceanEnergy)[0]
    
    IdxRadiousWind=All_IdxRadious[0:NumWindSites,0:NumWindSites]
    IdxRadiousWave=All_IdxRadious[NumWindSites:(NumWindSites+NumWaveSites),NumWindSites:(NumWindSites+NumWaveSites)]
    IdxRadiousOcean=All_IdxRadious[(NumWindSites+NumWaveSites):,(NumWindSites+NumWaveSites):]
    
    if NumWindSites!=0:
        KWind=All_K[0, 0:len(GrpupsWind)]
        SIdxInWind=GrpupsWind[np.where(KWind)[0][0]]
        IdxWindIn_Y=np.sum(IdxRadiousWind[SIdxInWind,:], axis=0)>=1
        
        WindEnergy=WindEnergy[IdxWindIn_Y,:]
        WindLatLong=WindLatLong[IdxWindIn_Y,:]
        AnnualizedCostWind=AnnualizedCostWind[IdxWindIn_Y]
        
    if NumWaveSites!=0:
        KWave=All_K[0, len(GrpupsWind):(len(GrpupsWind)+len(GrpupsWave))]
        SIdxInWave=GrpupsWave[np.where(KWave)[0][0]]-NumWindSites
        IdxWaveIn_Y=np.sum(IdxRadiousWave[SIdxInWave,:], axis=0)>=1
        
        WaveEnergy=WaveEnergy[IdxWaveIn_Y,:]
        WaveLatLong=WaveLatLong[IdxWaveIn_Y,:]
        AnnualizedCostWave=AnnualizedCostWave[IdxWaveIn_Y]
    
    if NumOceanSites!=0:
        KOcean=All_K[0, (len(GrpupsWind)+len(GrpupsWave)):(len(GrpupsWind)+len(GrpupsWave)+len(GrpupsOcean))]
        SIdxInOcean=GrpupsOcean[np.where(KOcean)[0][0]]-(NumWindSites+NumWaveSites)
        IdxOceanIn_Y=np.sum(IdxRadiousOcean[SIdxInOcean,:], axis=0)>=1
        
        OceanEnergy=OceanEnergy[IdxOceanIn_Y,:]
        OceanLatLong=OceanLatLong[IdxOceanIn_Y,:]
        AnnualizedCostOcean=AnnualizedCostOcean[IdxOceanIn_Y]
        RatedPowerOcean=RatedPowerOcean[IdxOceanIn_Y]
        
    return WindEnergy, WindLatLong, AnnualizedCostWind, WaveEnergy, WaveLatLong,\
             AnnualizedCostWave, OceanEnergy, OceanLatLong, AnnualizedCostOcean, RatedPowerOcean
             
    
