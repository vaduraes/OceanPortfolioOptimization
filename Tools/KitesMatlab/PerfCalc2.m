function [Power,aOpt,Fwing,CLCD] = PerfCalc2(SFOTParams,uGeo,uDep,vel)

% Geometric decision vairable set
Span = uGeo(1);         %wingspan (in m)
AR = uGeo(2);           %wing aspect ratio
Len = uGeo(3);          %fuselage length (in m)
Dia = uGeo(4);          %fuselage diameter (in m)

% Depth control decision variable set
theta = uDep(1);        %elevation angle (in deg)
lthr = uDep(2);         %tether length (in m)

% Velocity(vel) is an input

%% Hydrodynamic calculations 

lamda = Len/Dia;
L = Len*0.7; 

% Angle of atack 
AoA = -12.5:0.5:SFOTParams.stallAoA;                                               % upto stall angles (User defined) 
AoA = AoA.*(pi/180); 

% Areas  
% Sh = 2;                 %(for AR 8)
Sw = (Span^2)/AR; 
Swet = pi*Len*Dia*((1 - (2/lamda))^(2/3))*(1 + 1/(lamda^2)); 

% Choose tether diameter based on wing area 
%Diathr = min(SFOTParams.Diathr,0.8*thrDiaMap(Sw));
Diathr = SFOTParams.Diathr;
lthr_eff = lthr;
% if lthr < 600
%     lthr_eff = lthr;
% else
%     lthr_eff = 600 + (lthr - 600)/4;
% end
Sthr = Diathr*(lthr_eff/4);    %0.25 tether length straight assummption

% Wing 
slopeW = ((2*pi*SFOTParams.gammaw)./(1+((2.*pi.*SFOTParams.gammaw)./(pi.*SFOTParams.eLw.*AR))));
Clw = SFOTParams.Clw0 + slopeW.*AoA ; 

ClwD0 = 0.0026.*AR; 
Cd0w = (0.00015.*AR + 0.0053);
Cd0w = ones(1,numel(AoA)).*Cd0w; 
Cdw = (SFOTParams.Cdw_ind./AR+0.0334).*(Clw-ClwD0).^2 + Cd0w;  

% HStab 
Sh = ((SFOTParams.x_g-SFOTParams.h_sm)/(L+SFOTParams.h_sm-SFOTParams.x_g))*Sw;    % Area gotten by satisfying stability margins
Sh = max(Sh,2.5); % was min(Sh, 2.5); Mary 2/2/24
xeta = (SFOTParams.ClHStall*Sh)/((Sw*max(Clw)) + (SFOTParams.ClHStall*Sh));       % maximizing lift from H.Stabilizer 
Clh = (Sw*Clw*xeta)./(Sh*(1-xeta));                                               % Lift of H.Stabilizer 
Cdh = SFOTParams.Cdh_ovrall.*((Clh*(Sh/Sw)).^2) + SFOTParams.Cd0_h; 

% Fuselage 
CD_fuse = SFOTParams.Cfuse*(Swet/Sw);

% Net Lift 
CL = Clw + (Clh*(Sh/Sw));

% Tether Drag
CD_thr = SFOTParams.CDthr*(Sthr/Sw);

% Net Drag 
CD = Cdw + (Cdh*(Sh/Sw)) + CD_fuse + CD_thr ; 


%% EFFICIENCY MAP
p00 =       9.909;%  (-20.01, 39.83)
p10 =      -3.517;% (-13.43, 6.397)
p01 =      0.2191;%  (-0.453, 0.8911)
p20 =      0.4333;%  (-0.6641, 1.531)
p11 =    -0.05303;%  (-0.1889, 0.08285)
p02 =    0.003564;%  (-0.01284, 0.01997)
p30 =    -0.01758;%  (-0.05813, 0.02297)
p21 =    0.003423;%  (-0.003847, 0.01069)
p12 =  -0.0005983;%  (-0.002416, 0.001219)

u = [AR Span];

Eff = p00 + p10*u(2) + p01*u(1) + p20*u(2)^2 + p11*u(2)*u(1) + p02*u(1)^2 + p30*u(2)^3 + p21*u(2)^2*u(1) ...
            + p12*u(2)*u(1)^2;% + p03*u(1)^3 + p40*u(2)^4 + p31*u(2)^3*u(1) + p22*u(2)^2*u(1)^2 ... 
            %+ p13*u(2)*u(1)^3 + p04*u(1)^4;

% SFOTParams.eta = 0.85*Eff; %on the seabed
SFOTParams.eta = 1; %deploying on the surface

%% Calculate power and performance

[Perf, Ind] = max((CL.^3./CD.^2));                                             %find optimal AoA
aOpt = AoA(Ind)*(180/pi);                                                      %optimal AoA

CLCD = [CL(Ind) CD(Ind) SFOTParams.eta, ];

Power = (2/27)*SFOTParams.rho*SFOTParams.eta*Perf*Sw*(vel^3)*1e-3*((cosd(theta))^3); 
Fwing = 0.5*(Power/(vel/3))*(1e3);
end

