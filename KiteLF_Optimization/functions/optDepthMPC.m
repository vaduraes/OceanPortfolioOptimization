function [l_vec,theta_vec,uspl_vec,Jopt_vec] = optDepthMPC(uGeo,xSite,DataSet,SFOTParams)

% Geometric decision vairable set
% Span = uGeo(1);         %wingspan (in m)
% AR = uGeo(2);           %wing aspect ratio
% Len = uGeo(3);          %fuselage length (in m)
% Dia = uGeo(4);          %fuselage diameter (in m)
%% Parameter setup
% Define dmax for site and send it in through SFOTParams
SFOTParams.dmax = DataSet.dmax(xSite(1),xSite(2)); 

NtimeSteps = SFOTParams.NtimeSteps;
deltaT = 1;
t_hr = 1; 

% pre-allocation
l_vec = zeros(1,NtimeSteps); 
theta_vec = zeros(1,NtimeSteps); 
uspl_vec = zeros(1,NtimeSteps); 
Jopt_vec = zeros(1,NtimeSteps); 
%% Optimization setup
options = optimoptions('fmincon','Display',SFOTParams.fminconDisplay,'Algorithm','sqp',...
                    'MaxFunctionEvaluations',1e4,'MaxIterations',1e6,'OptimalityTolerance',1e-12);
% uDep = [thetaf lThrf uspl]
theta_ub = SFOTParams.thetaUB; 
lb = [15 50 0.1];
ub = [theta_ub min(2e3,SFOTParams.dmax/sind(50)) 1]; 

Jmulti = zeros(1,3);
uDepOptmulti = zeros(3,3);

%% Choose first time step at optimal (l, theta)
% u0 = [theta0 lThr0]
SFOTParams.NtimeSteps = 1;
[l_0,theta_0,Jopt_0] = findOptDepthMulti(uGeo,xSite,DataSet,SFOTParams);
u0 = [theta_0 SFOTParams.thetaUB]; 

SFOTParams.Pin = 5;
SFOTParams.NtimeSteps = NtimeSteps; 
%% Calculate uDep = [theta* lthr*] for each time period 
    for t = 1:NtimeSteps
        
        J = @(uDep)costCalc(uDep, u0, uGeo, xSite, DataSet, t_hr, SFOTParams);
        C = @(uDep)constraint(uDep, u0);
        if SFOTParams.MPCSolver == 1
                
        %Implement multi-start
        %1. Try lb
        uDep0 = lb;
        [uDepOpt,Jopt] = fmincon(J,uDep0,[],[],[],[],lb,ub,C,options);
        Jmulti(1) = Jopt; 
        uDepOptmulti(:,1) = uDepOpt'; 
        
        %2. Try ub
        uDep0 = ub;
        [uDepOpt,Jopt] = fmincon(J,uDep0,[],[],[],[],lb,ub,C,options);
        Jmulti(2) = Jopt; 
        uDepOptmulti(:,2) = uDepOpt'; 
        
        %3. Try somehwere in between
        uDep0 = [20 ub(2)/4 0.5];
        [uDepOpt,Jopt] = fmincon(J,uDep0,[],[],[],[],lb,ub,C,options);
        Jmulti(3) = Jopt; 
        uDepOptmulti(:,3) = uDepOpt'; 
        
        % Choose best
        [val,ind] = min(Jmulti);
        
        Jopt_vec(t) = -val; 
        theta_vec(t) = uDepOptmulti(1,ind);
        l_vec(t) = uDepOptmulti(2,ind);
        uspl_vec(t) = uDepOptmulti(3,ind);
        
        
        else 
        options = optimoptions('ga','Display', 'none');
        options.FunctionTolerance = 1e-3;
        options.MaxGenerations = 150;
        options.PopulationSize = 50; 
        uDep0 = [20 ub(2)/4 0.5];
        
        [uDepOpt,Jopt] = ga(J,length(uDep0),[],[],[],[],lb,ub,C,[],options);

        Jopt_vec(t) = -Jopt; 
        theta_vec(t) = uDepOpt(1);
        l_vec(t) = uDepOpt(2);
        uspl_vec(t) = uDepOpt(3);
        
        end 

        
        t_hr = t_hr + deltaT; 
        u0 = [uDepOpt(1) uDepOpt(2)]; 

        
    end 

end 

function [Cost] = costCalc(uDep, u0, uGeo, xSite, DataSet, t_hr, SFOTParams)

xInd = xSite(1); 
yInd = xSite(2);

% define operating depth
d_m = SFOTParams.dmax - uDep(2)*sind(uDep(1));
d_m = max(0,d_m); 
% obtain velocity at depth
[vel,flag] = getVelocityForXYDT(xInd,yInd,d_m,t_hr,DataSet); 

% Power generated
[Pgen,aOpt,Fwing] = PerfCalc(SFOTParams,uGeo,uDep,vel); 

% uDep = [thetaf lThrf uspl]
% u0 = [theta0 lThr0]

w = uDep(3); 

if w == 0
    uspl = 0; 
else
    uspl = (1/3600)*((uDep(2)-u0(2))/w);
end 

% w = (1/3600)*(abs(uDep(2)-u0(2))/uDep(3)); 

eta = 1-(0.85*(uspl^2)); 

if uDep(2)>=u0(2)
    Pspl = 0.85*Pgen; 
else 
    Pspl = (eta-1)*Pgen; 
end 

Power = (1-w)*Pgen + w*Pspl; 

Cost = -Power; 
end

function [c_ineq, c_eq] = constraint(uDep, u0)

% uDep = [thetaf lThrf uspl]
% u0 = [theta0 lThr0]
c_ineq = abs(u0(2)-uDep(2))-(uDep(3)*3600); 
c_eq = []; 
end 

