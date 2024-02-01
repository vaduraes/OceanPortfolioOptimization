function [l_vec,theta_vec,Jopt_vec] = findOptDepthMulti(uGeo,xSite,DataSet,SFOTParams)

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
options = optimoptions('fmincon','Display',SFOTParams.fminconDisplay,'Algorithm','sqp',...
                    'MaxFunctionEvaluations',1e4,'MaxIterations',1e6,'OptimalityTolerance',1e-12);
% uDep = [theta lThr]
theta_ub = SFOTParams.thetaUB; 
lb = [15 40];
% ub = [theta_ub min(2e3,SFOTParams.dmax/sind(50))]; 
ub = [theta_ub 3.5e3]; 
Jmulti = zeros(1,3);
uDepOptmulti = zeros(2,3);
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
        l_vec(t) = uDepOptmulti(2,ind);
        
        t_hr = t_hr + deltaT; 
    end 

end 

function [Cost] = costCalc(uDep, uGeo, xSite, DataSet, t_hr, SFOTParams)

xInd = xSite(1); 
yInd = xSite(2);

t_hr = min(8761,t_hr);

% define operating depth
d_m = SFOTParams.dmax - uDep(2)*sind(uDep(1));
d_m = max(0,d_m); 
% obtain velocity at depth
[vel,flag] = getVelocityForXYDT(xInd,yInd,d_m,t_hr,DataSet); 

% Power generated
[Power,aOpt,Fwing] = PerfCalc(SFOTParams,uGeo,uDep,vel); 

Cost = -Power; 
end
