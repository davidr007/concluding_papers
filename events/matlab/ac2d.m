function [frames,t,receiver1] = ac2d(grid)
% Simple finite difference solver for the  Acoustic wave equation  
%               p_tt   = c^2 p_xx + src
% on a 2-D regular grid.
%
% Note that the grid consists of a sub block centered within a larger
% computational block. The sub-block defines the region of interest. The
% larger computational block is used to ensure that undesired reflections
% do not disturb the solution. 
%
%INPUTS
% grid     [structure] defining the grid parameters as follows. 
%   grid.dx             [scalar] delta x spacing in meters
%   grid.sublength_x    [scalar] total x length of sub-block in meters
%   grid.total_length_x [scalar] total x length (must be larger than grid.sublength)
%   grid.x              [vector] x locations [0:grid.dx:grid.total_length_x]
%   grid.nx             [scalar] number of x grid points [length(grid.x)] 
%   grid.dz             [scalar] delta z spacing in meters
%   grid.sublength_z    [scalar] total z length of sub-block in meters
%   grid.total_length_x [scalar] total z length (must be larger than grid.sublength)
%   grid.x              [vector] z locations [0:grid.dx:grid.total_length_x]
%   grid.nz             [scalar] number of z grid points [length(grid.x)] 
%   grid.time = 15;     [scalar] desired interval of time to be sampled (in sec)
%   grid.c0             [matrix(grid.nz,grid.nx)] velocity model in m per sec
%                          often produced by GEN_RANDOM_MEDIUM
%   grid.source_xfrac   [scalar] see notes on source location below    
%   grid.source_dx_step [scalar] see notes on source location below
%   grid.source_zfrac   [scalar] see notes on source location below
%   grid.source_dz_step [scalar] see notes on source location below
%
% OUTPUTS:
%
%
% NOTES:
% 1) Source location 
% The source position is computed by the local function 
% LOC_LOCATE_SRC as follow:
%   x index position = 
%       left of sub block 
%               + index corresponding to grid.source_xfrac*grid.sublength_x 
%                   + index corresponding to grid.source_dx_step
%   z index position = 
%       0 (i.e. top of sub block at surface) 
%               + index corresponding to grid.source_xfrac*grid.sublength_x 
%                   + index corresponding to grid.source_dx_step
%
% David Robinson
% 27 June 2005
% Modified from Matlab script of Heiner Igel (see e-mail 16 June 2005). 




if nargin ==0
    grid.dx = 20;
    grid.sublength_x = 20000;
    grid.total_length_x = 90000;
    grid.x = 0:grid.dx:grid.total_length_x;
    grid.nx = length(grid.x);
    grid.dz = 20;
    grid.sublength_z = 20000;
    grid.total_length_z = 55000;
    grid.z = 0:grid.dz:grid.total_length_z;
    grid.nz = length(grid.z);
    grid.time = 15;
    grid.c0 = 6000-750+1500*rand(grid.nz,grid.nx); % gaussain 
    grid.source_xfrac = 1/4;
    grid.source_dx_step = 0;
    grid.source_zfrac = 1/4;
    grid.source_dz_step = 0;
end




% Basic parameters

% dx=20;  % define space increment to be 20m
% glength = 90000;
% x = 0:dx:glength;
% z = x;
% nx = length(x);
% nz = length(z);

% DR


eps=1.;     % wave to travel maximum of one (or half-see below) box at a time. 
isnap=10;   % images updated every isnap time step


p=zeros([grid.nz grid.nx]); pnew=p; pold=p; d2px=p; d2pz=p; 
c=p; c=c+grid.c0;% initialization of pressure fields
dt=.5*eps*min([grid.dx,grid.dz])/max(max(c)); % calculate time step from stability criterion (wave travel max of 0.5eps box per time step)
%dt = 0.005;
%dt =5*dt;  % we can increase the time step due to the strong correlation across grid cells
nt=grid.time/dt   %compute total number of time steps


% source time function
f0=8;     					% dominant frequency (DR - same as peak frequency??)
t=(1:nt)*dt; 				% initialize time axis
t0=2/f0;                    % time delay as a multiple of the dominant freq. 

%second attempt at a Ricker Wavelet (this time we do not differentiate). (not finished results too high) 
% src = 2*f0^2*pi^2*exp(-f0^2*pi^2*(t-t0).^2).*(2*f0^2*pi^2*(t-t0).^2-1);
% src = diff(src);

% First attempt at a Ricker Wavelet (event1_3 and event2_3)
src = (1-2*pi^2*f0^2*(t-t0).^2).*exp(-pi^2*f0^2*(t-t0).^2);
src=diff(src);							% time derivative to obtain Gaussian

% Original from Heiner Igel (event1_4 and event2_4)
% src=exp(-f0^2*(t-t0).*(t-t0));	% source time function
% src=diff(src);

[ind4srcloc_x, ind4srcloc_z,sub_x, sub_z]= LOC_locate_src(grid)

count =1;
% Time stepping
for it=1:nt-1,
      
   % FD
   
   disp(sprintf(' Time step : %i',it));
   
   for j=3:grid.nz-2,  %DR loop over x
      for k=3:grid.nx-2,  %DR loop over z
         d2pz(j,k)=(-1/12*p(j+2,k)+4/3*p(j+1,k)-5/2*p(j,k)+4/3*p(j-1,k)-1/12*p(j-2,k))/grid.dx^2; % space derivative DR in x
         d2px(j,k)=(-1/12*p(j,k+2)+4/3*p(j,k+1)-5/2*p(j,k)+4/3*p(j,k-1)-1/12*p(j,k-2))/grid.dx^2; % space derivative DR in z
      end
   end
        
   pnew=2*p-pold+c.*c.*(d2px+d2pz)*dt^2;                % time extrapolation
   %pnew(floor(nx/4)+4,floor(nz/4)+4)=pnew(floor(nx/4)+4,floor(nz/4)+4)+src(it)*dt^2;     % add perturbed source term
   pnew(ind4srcloc_z, ind4srcloc_x)=pnew(ind4srcloc_z, ind4srcloc_x)+src(it)*dt^2;     % add un-perturbed source term 

   pold=p;											% time levels
   p=pnew;
   center_x = floor(grid.nx/2)
   receiver1(it) = p(3,center_x);
   
   if rem(it,isnap)== 0,
   % Display
   fh = figure(1)
   subplot(2,2,1)
   imagesc(grid.x,grid.z,p), axis equal
   title(' FD ')
   
   subplot(2,2,3)
   imagesc(grid.x(sub_x(1):sub_x(2)),grid.z(sub_z(1):sub_z(2)),p(sub_z(1):sub_z(2),sub_x(1):sub_x(2))), axis equal
   title(' FD ')
   
   subplot(2,2,2), plot(p(3,:))  % added by DR - note top 2 rows and left 2 & right cols all zero 
   drawnow
   subplot(2,2,4)
   plot(receiver1)
   %input('press enter to continue')
   
   frames(count) =getframe(fh);
   count = count+1
   
   end
end



function [ind4srcloc_x, ind4srcloc_z, sub_x, sub_z]=LOC_locate_src(grid)
% Local function to locate the source position in term of an index
% The algorithm can be summarised as follows
% x index position = 
%       left of sub block 
%               + grid.source_xfrac*grid.sublength_x 
%                   + grid.source_dx_step
% x index position = 
%       0 (i.e. top of sub block at surface 
%               + grid.source_xfrac*grid.sublength_x 
%                   + grid.source_dx_step
%
% Also return the index coordinates for the sub block

subblock_lower_x = grid.total_length_x/((grid.total_length_x-grid.sublength_x)/2);
ind_subblock_lower_x = floor(length(grid.x)/subblock_lower_x);
frac_step_x = floor(grid.source_xfrac*grid.sublength_x/grid.dx);
ind4srcloc_x = ind_subblock_lower_x + frac_step_x + grid.source_dx_step;

sub_x = [ind_subblock_lower_x, ind_subblock_lower_x+floor(grid.sublength_x/grid.dx)];

ind_subblock_lower_z = 0;
frac_step_z = floor(grid.source_zfrac*grid.sublength_z/grid.dz);
ind4srcloc_z = ind_subblock_lower_z + frac_step_z + grid.source_dz_step;

sub_z = [1, floor(grid.sublength_z/grid.dz)]
