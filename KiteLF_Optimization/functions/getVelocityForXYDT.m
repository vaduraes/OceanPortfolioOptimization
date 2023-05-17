function [vMag,flag] = getVelocityForXYDT(xInd,yInd,d_m,t_hr,DataSet)

flag = 0; 
% Inputs ranges: 
%     x,y = 1,50
%     d = int 0,dmax for x,y
%     t = int 0,8760 for x,y,d

% Flags: 
% 1 = depth > max depth at location 
% 2 = location unusable, max depth ~ 0


% Check if d is within range of d at (x,y) and location is usable
    if d_m > DataSet.dmax(xInd,yInd)
        flag = 1;
        fprintf('Depth > max depth at location')
        
        
    elseif DataSet.dmax(xInd,yInd) < 5
        flag = 2;
        fprintf('Unusable location, max depth ~0')
        
    end 

% Get d,v vectors for x,y,t and interpolate based on specified depth
dVec = squeeze(DataSet.depth(xInd,yInd,:));
vVec = squeeze(DataSet.vel(xInd,yInd,:,t_hr));

vMag = interpn(dVec,vVec,d_m);
vMag = min(vMag,2.5);
    
end 