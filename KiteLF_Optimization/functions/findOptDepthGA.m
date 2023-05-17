function [l_vec,theta_vec,Jopt_vec] = findOptDepthGA(uGeo,xSite,DataSet,SFOTParams)

% Geometric decision vairable set
% Span = uGeo(1);         %wingspan (in m)
% AR = uGeo(2);           %wing aspect ratio
% Len = uGeo(3);          %fuselage length (in m)
% Dia = uGeo(4);          %fuselage diameter (in m)
%% Parameter setup
% Define dmax for site and send it in through SFOTParams
SFOTParams.dmax = DataSet.dmax(xSite(1),xSite(2)); 

NtimeSteps = SFOTParams.NtimeSteps;
deltaT = round(8761/NtimeSteps);
t_hr = 1; 

% pre-allocation
l_vec = zeros(1,NtimeSteps); 
theta_vec = zeros(1,NtimeSteps); 
Jopt_vec = zeros(1,NtimeSteps); 
%% Optimization setup
% options = optimoptions('ga','PlotFcn',@gaplotbestf);
options = optimoptions('ga','Display', 'none');
options.FunctionTolerance = 1e-3;
options.MaxGenerations = 150;
options.PopulationSize = 50;
% uDep = [theta lThr]
theta_ub = 30; 
lb = [15 50];
ub = [theta_ub SFOTParams.dmax/sind(theta_ub)]; 
uDep0 = [20 ub(2)/2];
%% Calculate uDep = [theta* lthr*] for each time period 
    for t = 1:NtimeSteps
        J = @(uDep)costCalc(uDep, uGeo, xSite, DataSet, t_hr, SFOTParams);
        [uDepOpt,Jopt] = ga(J,length(uDep0),[],[],[],[],lb,ub,[],[],options);
        
        theta_vec(t) = uDepOpt(1);
        l_vec(t) = uDepOpt(2);
        Jopt_vec(t) = -Jopt; 
        
        t_hr = t_hr + deltaT; 
    end 

end 

function [Cost] = costCalc(uDep, uGeo, xSite, DataSet, t_hr, SFOTParams)

xInd = xSite(1); 
yInd = xSite(2);

% define operating depth
d_m = SFOTParams.dmax - uDep(2)*sind(uDep(1));

% obtain velocity at depth
[vel,flag] = getVelocityForXYDT(xInd,yInd,d_m,t_hr,DataSet); 

% Power generated
[Power,aOpt,Fwing] = PerfCalc(SFOTParams,uGeo,uDep,vel); 

Cost = -Power; 
end
