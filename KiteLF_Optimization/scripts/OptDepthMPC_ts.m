
xSiteVec = [linspace(17,49,33)' linspace(1,33,33)';linspace(16,49,34)' linspace(1,34,34)'; ...
            linspace(17,49,33)' linspace(3,35,33)'];

%% Testig functions

setSFTParams; 

uGeo = [9.2 4.6 7 0.48];
SFOTParams.Diathr = 0.012;

for i = 1


% xInd = xSiteVec(i,1);
% yInd = xSiteVec(i,2);


xInd = 46;
yInd = 33;

xSite = [xInd yInd];
uGeo = [sVec(i) ARVec(i) 7 0.5];
d_m = 40;
t_hr = 100;

% Remove while running sweep
% uGeo = [10 3.5 7 0.5];
SFOTParams.fminconDisplay = 'none';

%% HF: Static opt
SFOTParams.Diathr = 0.012;
SFOTParams.thetaUB = 50;
tic
SFOTParams.NtimeSteps = 8761;
SFOTParams.MPCSolver = 1;   %1 = fmin, 2 = ga    
[l_vec,theta_vec,uspl_vec,Jopt_vec] = optDepthMPC(uGeo,xSite,DataSet,SFOTParams);
mean(Jopt_vec)
toc

%% Save data

MPCResults.l_vec = l_vec; 
MPCResults.theta_vec = theta_vec; 
MPCResults.uspl_vec = uspl_vec;
MPCResults.Jopt_vec = Jopt_vec; 

figure(1)
hold on 
plot(Jopt_vec)
% 
pMPC = mean(Jopt_vec); 

%% MF: Optimum chasing 1

tic
xSite = [xInd yInd];
SFOTParams.NtimeSteps = 150;
[l_vec,theta_vec,Jopt_vec] = findOptDepthMulti(uGeo,xSite,DataSet,SFOTParams);
mean(Jopt_vec)
toc


%% Make MF the right length

NSteps = length(Jopt_vec); 
Total = 8761; 
deltaT = round(Total/NSteps);
clear t; 

Jvec = zeros(1,Total);
lvec = zeros(1,Total); 
thvec = zeros(1,Total); 

t = 1; 
for tInd = 1:NSteps
    
    if tInd == NSteps
        Jvec(t:end) = Jopt_vec(tInd);
        lvec(t:end) = l_vec(tInd);
        thvec(t:end) = theta_vec(tInd);
    else 
        Jvec(t:t+deltaT) = Jopt_vec(tInd);
        lvec(t:t+deltaT) = l_vec(tInd);
        thvec(t:t+deltaT) = theta_vec(tInd);
    end 
    t = t + deltaT; 
end


%% Save data

DOptResults.l_vec = lvec; 
DOptResults.theta_vec = thvec; 
DOptResults.Jopt_vec = Jvec; 

%% Plot 

figure(1)
hold on 
plot(Jvec,'LineWidth',1.5)

pDOpt = mean(Jvec); 


%% LF: Constant vel

tic
uDep = [DataSet.meanTh(xInd,yInd) DataSet.meanL(xInd,yInd)]; 
vel = DataSet.vOpD(xInd,yInd); 
[PowerLF,aOpt,Fwing] = PerfCalc(SFOTParams,uGeo,uDep,vel);
pLFOpt = PowerLF*ones(1,8761);
toc
%% Plot 

figure(1)
hold on 
plot(pLFOpt)

%% Save results
LFResults.Jopt_vec = pLFOpt; 
LFResults.velocity = vel; 
LFResults.meanPowerVec = [pMPC pDOpt PowerLF]; 
LFResults.LFHFerr = rmseErr(MPCResults.Jopt_vec,pLFOpt); 
LFResults.LFMFerr = rmseErr(MPCResults.Jopt_vec,DOptResults.Jopt_vec); 

disp([pMPC pDOpt PowerLF])
disp([LFResults.LFHFerr/pMPC LFResults.LFMFerr/pMPC])

%% Save results

fileName = sprintf('x%dy%d_s%0.1fAR%0.1f_mltiFid.mat',xInd,yInd,uGeo(1),uGeo(2));
save(fileName,'MPCResults','DOptResults','LFResults');

fprintf('Evaluation %d of 100 \n',i)

end 