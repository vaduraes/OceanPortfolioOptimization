%% Set airfoil shape and geometry

[xaf,usaf,lsaf] = airfoil_data();
SWDTParams.xaf = xaf;
SWDTParams.usaf = usaf;
SWDTParams.lsaf = lsaf;

SWDTParams.ChrdT = 0.12;

% Material properties
SWDTParams.E = 69*(10^9); %Aluminium
SWDTParams.rhow = 2710.0; %Al

% Percentage Deflection at centroid
SWDTParams.defper = 5.0;
SWDTParams.NSparmax = 3.0;