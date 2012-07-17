function [profile] = make_profiles_from2DPDF(grid_xy, PDF_2D)
% This function extracts some profiles form a 2D PDF. 
% Note that it is does not conduct and normalization
% Note that it assumes the grid in the x and y directions are identical
%
% INPUTS:
% grid_xy   [vector 1 x n] row vector containing grid values in both the 
%           x and y direction 
% PDF_2D    [matrix nxn] 2D PDF at x and y grid values given by grid_xy
%
% OUTPUTS: 
% profile   [matrix n x4] matrix containing two profiles. Profile1 is
%           horizontal through the center and Profile2 is diagonal 
%           through the center
%               profile(:,1) contains x values for profile1
%               profile(:,2) contains PDF values for profile1
%               profile(:,3) contains position along diagonal for profile2
%               profile(:,4) contains PDF values for profile2
%
% David Robinson
% 1 April 2008

% find where the horizontal through the center is located
ind = find(grid_xy==0) 
% Have a first pass at making the profiles
profile = [grid_xy',PDF_2D(ind,:)', sqrt(2*grid_xy'.^2),diag(flipud(PDF_2D))];
% Correct the third colum for negative 'separations' 
ind2 = find(profile(:,1)<0);
profile(ind2,3) = -profile(ind2,3);


