%function[l_vec,theta_vec,Jopt_vec] = findOptDepth(uGeo,xSite,DataSet,SFOTParams)
%Uncomment above line if you want to use this as a function. Be sure to
%also uncomment specified end (it is labeled)
%Jopt_vec is the power output indexed by time step

%% Load Data Set
addpath("./functions/")
addpath("./scripts/")
addpath("./functions/SWDT/")
File='DataSetPlatform.mat';
load(File) %%%%%Input your data set here
% Set xSite and uGeo %%%%%Can comment out if already have a vector for xsite, can also use uopt(PMInd,1:4) to populate uGeo.
xSite = [1,1];
uGeo = [11,3,7,0.4,18.879664767349897]; %See 3 lines below for what these are


% Geometric decision vairable set
Span = uGeo(1);         %wingspan (in m) 
AR = uGeo(2);           %wing aspect ratio 
Len = uGeo(3);          %fuselage length (in m)
Dia = uGeo(4);          %fuselage diameter (in m)

%% Set Rated Power%%%%% From Ben's Code
ratedPWR = 150000.1097;
uGeo(5)=ratedPWR;

%% Parameter setup


%%%%%%SFOT SETUP (COMMENT OUT IF YOU HAVE SFOTParams/want to use this as a function)
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
SFOTParams.Diathr = 0.012; %Updated For Mary


SFOTParams.eta = 0.3; 
SFOTParams.netaV = 0.8; 
SFOTParams.rho = 1000.0;

% Optimal depth calculation
SFOTParams.NtimeSteps = 87; %Number of hours in each step of simulation

SFOTParams.fminconDisplay = 'iter'; 

%%%%% END SFOT Setup Comment Out


% Define dmax for site and send it in through SFOTParams
SFOTParams.dmax = DataSet.dmax(xSite(1),xSite(2)); %%%%%Choose from Data set loaded from above



%IMPORTANT% Define theta (AoA) upper bound READ INSTRUCTIONS BELOW
SFOTParams.thetaUB = 50; %%%%%Make sure proper angle of attack bound. User/geometry decision but could be in Preexisting SFOT parameters, if using a predefined SFOT array, check to make sure this is included


NtimeSteps = SFOTParams.NtimeSteps; 
deltaT = round(8761/NtimeSteps);
t_hr = 1; 

% pre-allocation
l_vec = zeros(1,NtimeSteps); 
theta_vec = zeros(1,NtimeSteps); 
Jopt_vec = zeros(1,NtimeSteps); 

%% Optimization setup

options = optimoptions('fmincon','Display',SFOTParams.fminconDisplay,'Algorithm','sqp',...
                    'MaxFunctionEvaluations',1e4,'MaxIterations',1e6,'OptimalityTolerance',1e-12);
%uDep = [theta lThr]
theta_ub = SFOTParams.thetaUB; 
% lb = [15 40]; %Noticed problems with these feasibility bounds lb and ub
lb = [15 20/sind(50)];

%ub = [theta_ub SFOTParams.dmax/sind(theta_ub)]; 
% ub = [theta_ub 2e3];
ub = [50 300/sind(15)];

Jmulti = zeros(1,3);
uDepOptmulti = zeros(2,3);

%uDep0 = [20 0.1*ub(2)/2];

%% Calculate uDep = [theta* lthr*] for each time period 
    for t = 1:NtimeSteps
        J = @(uDep)costCalc(uDep, uGeo, xSite, DataSet, t_hr, SFOTParams);
       
        %Implement multi-start
        %1. Try lb
        uDep0 = lb;
        [uDepOpt,Jopt] = fmincon(J,uDep0,[],[],[],[],lb,ub,[],options);
        Jmulti(1) = Jopt; 
        uDepOptmulti(:,1) = uDepOpt'; 
        
        %2. Try ub
        uDep0 = ub;
        [uDepOpt,Jopt] = fmincon(J,uDep0,[],[],[],[],lb,ub,[],options);
        Jmulti(2) = Jopt; 
        uDepOptmulti(:,2) = uDepOpt'; 
        
        %3. Try somehwere in between
        uDep0 = [20 ub(2)/4];
        [uDepOpt,Jopt] = fmincon(J,uDep0,[],[],[],[],lb,ub,[],options);
        Jmulti(3) = Jopt; 
        uDepOptmulti(:,3) = uDepOpt'; 
        
        % Choose best
        [val ind] = min(Jmulti);
        
        Jopt_vec(t) = -val; 
        theta_vec(t) = uDepOptmulti(1,ind);
        l_vec(t) = uDepOptmulti(2,ind); %tether length (in m)
        
        t_hr = t_hr + deltaT; 
    end
    
%end %%%%%To be uncommented if using as a function

%% Function
function [Cost] = costCalc(uDep, uGeo, xSite, DataSet, t_hr, SFOTParams)%%%%%Needs data from data set from above
xInd = xSite(1); 
yInd = xSite(2);

% define operating depth

%if File=="DataSetSeabed.mat"
%d_m = SFOTParams.dmax - uDep(2)*sind(uDep(1));


%"DataSetPlatform.mat"
d_m = uDep(2)*sind(uDep(1));


% obtain velocity at depth
[vel,flag] = getVelocityForXYDT(xInd,yInd,d_m,t_hr,DataSet); 

% Power generated
[Power,aOpt,Fwing,~] = PerfCalc2(SFOTParams,uGeo,uDep,vel);

%Limit generation to rated power
% Power(Power>uGeo(5))=uGeo(5);
Cost = -Power;

end
