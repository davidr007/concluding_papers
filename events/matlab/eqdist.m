function distance = eqdist(lat1,lon1,depth1,lat2,lon2,depth2)
% computes the euclidean distance between two earthquakes
%
% INPUTS:
% lat1      latitude of earthquake 1 (degrees)
% lon1      longitude of earthquake 1 (degrees)
% depth1    depth of earthquake 1 from ellipsoid (meters)
%               0 =>    surface of ellipsoid
%               +ve =>  below ellipsoid (most common)
%               -ve =>  above elipsoid   
% lat2      latitude of earthquake 2 (degrees)
% lon2      longitude of earthquake 2 (degrees)
% depth2    depth of earthquake 2 from ellipsoid 
%
% OUTPUTS
% distance  euclidean distance (meters)
[x1,y1,z1] = llh2xyz(lat1,lon1,-depth1);
[x2,y2,z2] = llh2xyz(lat2,lon2,-depth2);
distance = sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2);
