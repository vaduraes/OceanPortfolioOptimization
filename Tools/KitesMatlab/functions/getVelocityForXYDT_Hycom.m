function [vMag,flag] = getVelocityForXYDT(iSite,d_m,t_hr,DataSet)

flag = 0; 
% Inputs ranges: 
%     x,y = 1,50
%     d = int 0,dmax for x,y
%     t = int 0,8760 for x,y,d

% Flags: 
% 1 = depth > max depth at location 
% 2 = location unusable, max depth ~ 0


% Check if d is within range of d at (x,y) and location is usable
    if d_m > DataSet.Dmax(iSite)
        flag = 1;
        fprintf('Depth > max depth at location')
        
        
    elseif DataSet.Dmax(iSite) < 5
        flag = 2;
        fprintf('Unusable location, max depth ~0')
        
    end 

% Get d,v vectors for x,y,t and interpolate based on specified depth
dVec = squeeze(DataSet.depth(:));% From hycom same depth profile for all sites
vVec = squeeze(DataSet.OCSpeed(iSite,:,t_hr));

vMag = interpn(dVec,vVec,d_m);
% vMag = min(vMag,2.5); %what this is doing mary? Limiting current speed to 2.5?
    
end 