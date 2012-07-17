function [x,y,z]=llh2xyz(lat,lon,height)
% function to convert latitude and longitude and height from an ellipsoid
% into cartesian coordiantes (x,y,z). 
%
% Uses WGS84 for:
%       AE = 6378137.0
%       FLAT = 1./298.257222101
%
% INPUTS:
% lat       latitude on surface of ellipsoid (degrees)
% lon       longitude on surface of ellipsoid (degrees)
% height    height above centroid (meteres)
%           height<0 => depth
%
% Outputs:
% x,y,z     cartesian coordiantes (meters)

%setup
AE = 6378137.0;
FLAT = 1./298.257222101;
rlat = lat*pi/180;
rlon = lon*pi/180;

% do it
e2 = 2.0 * FLAT - FLAT^2;
s2 = sin(rlat)^2;
N = AE / (sqrt(1.0 - (e2 * s2)));
x = cos(rlat) * (height + N) * cos(rlon);
y = cos(rlat) * (height + N) * sin(rlon);
z = sin(rlat) * (height + N * (1.0 - e2));


