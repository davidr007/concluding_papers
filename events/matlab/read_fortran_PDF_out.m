function [PDF_2D,grid_x2,grid_y2] = read_fortran_PDF_out(fname,dname)
% Reads the fortran output of a 2D PDF (columnwise) from 
% fname and converts it to a 2D PDF. 
%
% INPUTS: 
% fname     [string] filename
%               default = 'example_3eq_1.txt'
% dname     [string] directory containing fname
%               default =
%               '/export/storage/davidr/sandpit/davidr/events/fortran/'
% 
% OUTPUTS:
% PDF_2D    [matrix ny x nx] The desired 2D PDF
% grid_x2   [vector 1 x nx] grid values for x-dimension
% grid_y2   [vector ny x 1] grid values for y-dimension
%
% David Robinson 
% 18 March 2008

%% Sort out the default inputs
if nargin ==0
    fname = 'example_3eq_1.txt';
    dname = pwd;%'/export/storage/davidr/sandpit/davidr/events/fortran/';
end

if isempty(dname)
    dname = pwd;%'/export/storage/davidr/sandpit/davidr/events/fortran/';
end

%% Now do the work

tmp = load([dname,fname]);
nx = tmp(1); 
deltax = tmp(2);
ny = tmp(3);
deltay = tmp(4);

PDF_2D = reshape(tmp(5:end), ny,nx);

nx_half = floor(nx/2);
upper_xlim = nx_half*deltax;
grid_x2 = linspace(-upper_xlim,upper_xlim,nx);

ny_half = floor(ny/2);
upper_ylim = ny_half*deltay;
grid_y2 = linspace(-upper_ylim,upper_ylim,ny);






