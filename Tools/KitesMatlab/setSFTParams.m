
SFOTParams.gammaw = 0.95; 
SFOTParams.eLw = 0.71;
SFOTParams.Clw0 = 0.16; 
SFOTParams.Lscale = 1; 
SFOTParams.x_g = 0*SFOTParams.Lscale;
SFOTParams.h_sm = -0.2*SFOTParams.Lscale; 
SFOTParams.Cd0_h = 1.7*(10^-4);
SFOTParams.Cdw_visc = 0.03; 
SFOTParams.Cdw_ind = 1/(pi*0.98);
SFOTParams.Cdh_ovrall =  0.16;  
SFOTParams.Cfuse = 0.006; 
SFOTParams.stallAoA = 14; 
SFOTParams.ClHStall = 1.3;

% Tether
SFOTParams.CDthr = 1; 
SFOTParams.Diathr = 0.05;


SFOTParams.eta = 0.3; 
SFOTParams.netaV = 0.8; 
SFOTParams.rho = 1000.0;

% Optimal depth calculation
SFOTParams.NtimeSteps = 20;

SFOTParams.fminconDisplay = 'iter'; 
