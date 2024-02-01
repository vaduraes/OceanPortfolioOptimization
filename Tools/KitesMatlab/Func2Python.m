%% Loading: Run ONCE at start ONLY
xInd = 50;
yInd = 50;
xSite = [xInd yInd];

load('DataSetPlatform.mat'); %platform for ocean surface doployment, seabed for seabed deployment
setSFDTParams; 
setSWDTParams; 
setSFTParams; 

clc

%% Specify site index



global Fmax
Fmax = 1;

%% Plant decision variable bounds and initial guess

% uPlant = [s AR L D tsh tsp Nsp tfs dThr]
lb = [7   3  7 0.4 1e-4 0.01 1 1e-4 5e-3];
ub = [11 10  10 0.6 0.01 0.15  3 0.05 4e-2]; 

u0 = 0.5*(lb+ub);
u0(1:2) = [10 4]; 

%% Run multiObj optimization 

rng default
opts = optimoptions("gamultiobj","PlotFcn","gaplotpareto","ParetoFraction",0.3); 

J = @(u)Plant_cost(SFOTParams,SWDTParams,SFDTParams,xSite,DataSet,u);
C = @(u)Plant_constraint(SFDTParams,SWDTParams,u);

[uopt,CostOpt,conv_flag] = gamultiobj(J, 9, [], [], [], [], lb, ub, C, opts);

%% Get best P/M value 

power = (-1*CostOpt(:,1));
yearlyEnergy = power*24*365.25;
mass = (CostOpt(:,2)); 
cost = getCostFromMass(mass);
[LCOEBest,LCOEInd] = min(cost./yearlyEnergy); 

uGeo = uopt(LCOEInd,1:4);
 
uDep = [DataSet.meanTh(xSite(1),xSite(2)) DataSet.meanL(xSite(1),xSite(2))]; 
vel = DataSet.vOpD(xSite(1),xSite(2)); 

[~,~,~,CLCD] = PerfCalc2(SFOTParams,uGeo,uDep,vel);

%% Display results

mFuse = SWDTParams.rhow*2*pi*uopt(LCOEInd,4)*uopt(LCOEInd,3)*uopt(LCOEInd,8);
mWing = mass(LCOEInd) - mFuse;
Power = power(LCOEInd);
uopt = uopt(LCOEInd,:);
uopt([5 6 8 9]) = uopt([5 6 8 9])*1000;


%% Actual Display

disp(' u     = [ s(m)   AR    L(m)   D(m)   t_sh(mm)   t_sp(mm)   Nsp    T_fs(mm)   d_thr(mm)]')
fprintf(' u_opt = [%.2f  %.2f   %.2f   %.2f    %.2f        %.2f    %.2f    %.2f        %.2f  ]\n\n',uopt);


disp('    mFuse     mWing     Power(kW)   LCOE($/kWh)    C_L        C_D      eta')
fprintf('   %.2f    %.2f      %.2f       %.3f       %.4f    %.4f    %.3f\n',[mFuse mWing Power LCOEBest CLCD])

%% Functions

function [Cost] = Plant_cost(SFOTParams,SWDTParams,SFDTParams,xSite,DataSet,u)
global Fmax mKite

if isnan(Fmax)
    Fmax = 1; 
end 

% Estimate Power 
uGeo = u(1:4);
 
uDep = [DataSet.meanTh(xSite(1),xSite(2)) DataSet.meanL(xSite(1),xSite(2))]; 
vel = DataSet.vOpD(xSite(1),xSite(2)); 

SFOTParams.Diathr = 0.012; 

[PowerLF,aOpt,Fwing,~] = PerfCalc(SFOTParams,uGeo,uDep,vel);

Pavg = PowerLF; 
Fmax = Fwing;

% Calculate Ixx for a given chord and  wing structural variables
nsp = round(u(7));
Chord = u(1)/u(2);
tsh = u(5);
tsp = u(6); 

[~, area_tot]= Wing_MoIArea(SWDTParams,Chord,tsh,nsp,tsp);
 

% Calculate wing mass 
mWing = SWDTParams.rhow*(area_tot)*u(1);

% Calculate fuselage mass 
mFuse = SWDTParams.rhow*2*pi*u(4)*u(3)*u(8); 

% Calculate kite mass 
mKite = mWing + mFuse; 

% Cost 
Cost = [-Pavg,mKite]; 
end 

function [c_ineq,c_eq] = Plant_constraint(SFDTParams,SWDTParams,u)
global Fmax mKite

if isnan(Fmax)
    Fmax = 1; 
end 

%% Wing structural constraints 

% Calculate required Ixx for the maximum foces on wing
halfSpan = u(1)/2;
delx = SWDTParams.defper*halfSpan/100.0;
Ixx_req = (39.37^4)*5*(0.9)*Fmax*(halfSpan^3)/(48*SWDTParams.E*delx);

% Calculate Ixx for a given chord and  wing structural variables
nsp = round(u(7));
Chord = u(1)/u(2);
tsh = u(5);
tsp = u(6); 
[Ixx_tot, ~]= Wing_MoIArea(SWDTParams,Chord,tsh,nsp,tsp);


ineq1 =  Ixx_req - (Ixx_tot*(39.37^4));


%% Fuselage calculations (messy!)

Forces.FzW = Fmax*0.9; 
Forces.FzH = Fmax*0.1; 

Inp.HStab_c = SFDTParams.C_Hstab;          %HStab chord 

% Positions of tether attachment points 
Inp.x1 = SFDTParams.x_1; 
Inp.x2 = SFDTParams.x_2;               %from the tail 

Inp.xW = SFDTParams.x_W;              %position (%) of the wing 

Inp.fos = SFDTParams.FOS;              %factor of safety 
Inp.Syield = SFDTParams.S_yield;         %yield stress 

% Internal, external and dynamic pressures
Inp.IntP = SFDTParams.Int_P; 
Inp.ExtP = SFDTParams.Ext_P; 
Inp.DynP = SFDTParams.Dyn_P; 

% Densities  
Inp.rhow = 1e3; 
Inp.rhoAl = 2710; 

Inp.tarBuoy = 0.98; 

% Fuselage is not a cylinder 
effV = SFDTParams.effV; 


% Reaction forces 
SolRHS = [1 1; (Inp.x1 - Inp.xW)*u(3) (1-Inp.x2-Inp.xW)*u(3)]; 
SolLHS = [Forces.FzW + Forces.FzH ; Forces.FzH*(1-Inp.xW-0.75*Inp.HStab_c)*u(3)]; 
Forces.Rx_vec = inv(SolRHS)*SolLHS; 


r = 0.5*u(4); 
l = u(3); 

% Length from shear stress 
ineq2 = (Forces.FzW+Forces.FzH)/(u(8)*l) - (1/Inp.fos)*Inp.Syield; 

% Radius for hoop stress due to shear  
Pdiff = Inp.ExtP + Inp.DynP - Inp.IntP; 
ineq3 = (Pdiff*r)/u(8) - (1/Inp.fos)*Inp.Syield; 

% Buckling of fuselage  
BMmax = (Inp.x1*l*Forces.Rx_vec(1))-(Inp.xW*l*Forces.FzW)-((1-Inp.x2)*l*Forces.Rx_vec(2)); %max bending moment
SMod = pi*r^2*u(8); 
ineq4 = abs(BMmax/SMod) - (1/Inp.fos)*Inp.Syield;

% Neutral buoyancy (Makes solver significantly slower)
ineq5 = mKite - 1e3*((0.15*Chord*u(1)) + pi*u(4)^2*u(3));

% Tether strength constraint
% fos = 1.5;
% ineq3 = abs(Fmax*fos - 1e9*pi*0.25*(u(9)^2)) - 10;

c_ineq = [ineq1; ineq2; ineq3; ineq4; ineq5]; 
c_eq = []; 
end 