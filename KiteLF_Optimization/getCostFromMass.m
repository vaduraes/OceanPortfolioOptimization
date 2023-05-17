function cost = getCostFromMass(mass)
%using numbers from SHARKS 200kW Utility scale, computes yearly cost as a
%function of system mass

Csref = 2;
f1 = [8/Csref; 2/Csref*ones(3,1); 1; 0.6/Csref];
scaleRat = [0.20263518487989; 0.508112208638595; 1; ... 
    0.20263518487989; 0.999999703428232; 0.318932811324866];
f2 = scaleRat.*[2.80066; 2.80066; 8.03236; ...
    7.83725; 1.72207; 1.72207];
f3 = scaleRat.*[ones(4,1)*0.06365; ones(2,1)*0.11747];
opex = 293.673*200;
FCR = 0.082;

cost = nan(size(mass));
for i = 1:length(mass)
    m = [mass(i); 13445.51*2; 0; 6587.09; 2500; 106100.8];
    meq = sum(m.*f1.*(1+f2+f3));
    capex = Csref*meq;
    cost(i) = opex + capex*FCR;
end

end