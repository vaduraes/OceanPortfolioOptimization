function [Ixx_tot, area_tot]= Wing_MoIArea(SWDTParams,Chord,tsh,nsp,tsp)
% OLD: tsh, tle,tte, nsp, tsp1, tsp2, tsp3, locsp1,locsp2,locsp3
% OLD: Wing_MoIArea(MTParams,u(2)/u(1),u(3),0,0,NSp,u(4),u(4),u(4),15,30,60);

global xaf usaf lsaf
global Naf
Naf = 800;

xaf = SWDTParams.xaf;
usaf = SWDTParams.usaf;
lsaf = SWDTParams.lsaf;

locsp1 = 15; 
locsp2 = 30;
locsp3 = 60;

tsp1 = tsp;
tsp2 = tsp;
tsp3 = tsp;

C = Chord;

%% First step: Shell Area and MoI calculation
tsh = tsh*C;
area_sh = 0;
Ixx_shout = 0;
Ixx_shin = 0;

for j = 1:1:Naf-1
    b_sh = C*(xaf(j+1) - xaf(j));
   
    y_us = C*usaf(j+1); y_ls = C*lsaf(j+1);
    if usaf(j) > usaf(j+1)
        y_us = C*usaf(j);
    end
    if abs(lsaf(j)) > abs(lsaf(j+1))
        y_ls = C*lsaf(j);
    end
    d_shout = y_us - y_ls;
    Ixx_shout = Ixx_shout + (b_sh.*(d_shout^3.0)/12.0);
    
    d_shin = d_shout - 2*(tsh);
    
    if d_shout > tsh
        Ixx_shin = Ixx_shin + (b_sh.*(d_shin^3.0)/12.0);
        area_sh = area_sh + 2*(b_sh*tsh);
    else
        area_sh = area_sh + 2*(b_sh*d_shout);
    end
end

Ixx_sh = (Ixx_shout - Ixx_shin);


%% Second Step: Spar Area and MoI calculations
Sp1T = int16(tsp1*Naf); % front spar
Sp2T = int16(tsp2*Naf); % middle spar
Sp3T = int16(tsp3*Naf); % aft spar

% Spar MoI calculation function
if nsp == 0
    Sp1T = 0; Sp2T = 0; Sp3T = 0;
elseif nsp == 1
    Sp1T = 0;Sp3T = 0; 
end


x_fore_ori = locsp1*Naf/100;
x_mid_ori  = locsp2*Naf/100;
x_aft_ori  = locsp3*Naf/100;

[Ixx_spar1, area_spar1] = SparCalc(C,Sp1T,x_mid_ori,tsh);
[Ixx_spar2, area_spar2] = SparCalc(C,Sp2T,x_fore_ori,tsh);
[Ixx_spar3, area_spar3] = SparCalc(C,Sp3T,x_aft_ori,tsh);

%% Third Step: Total MoI and Area Calculation
Ixx_tot = (Ixx_sh + Ixx_spar1 + Ixx_spar2 + Ixx_spar3);

area_tot = (area_spar1 + area_spar2 + area_spar3 + area_sh);
end

function [Ixx_spar,area_spar] = SparCalc(C,SpT,x_sp_ori,tsh)
global xaf usaf lsaf

% Assumption - spars are made of 2 trapezoids
x_sp_ll  = x_sp_ori  - (SpT/2); 
x_sp_ul  = x_sp_ori  + (SpT/2);

a_ll =  usaf(x_sp_ll)  -  lsaf(x_sp_ll)  - 2*tsh;
a_ul =  usaf(x_sp_ul)  -  lsaf(x_sp_ul)  - 2*tsh;
b    =  usaf(x_sp_ori) -  lsaf(x_sp_ori) - 2*tsh;
h1   =  xaf(x_sp_ori)  -  xaf(x_sp_ll);
h2   =  xaf(x_sp_ul)   -  xaf(x_sp_ori);

Ixx_spar1 = (C^4)*(h1^3)*(3*a_ll + b)/12.0;
Ixx_spar2 = (C^4)*(h2^3)*(3*a_ul + b)/12.0;

area_spar1 = 0.5*(a_ll + b)*h1*(C^2);
area_spar2 = 0.5*(a_ul + b)*h2*(C^2);

Ixx_spar = (Ixx_spar1 + Ixx_spar2);
area_spar = area_spar1 + area_spar2;
end

