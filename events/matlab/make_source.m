function [dt,stf] = make_source(fname,tlength,dx, max_vel, pwidth,src_posy,strike,dip,rake,Mo,ttrans,extras)

% INPUTS:
% fname     [string] filename for source time function file
% tlength   [scalar] length of desired waveform (seconds)
% dx        [scalar] grid size (metres)
% max_vel   [scalar] approximate maximum shear wave velocity for
%           calculating stability criterion
% src_posy  [vector 1x3]
%               1) src_posy(1) = x position of source (index)
%               2) src_posy(2) = y position of source (index)
%               3) src_posy(3) = z position of source (index)
% pwidth    [scalar] width of stf pulse (stanadrd deviation for gaussain)
% strike    [scalar] strike in degrees (theta)
% dip       [scalar] dip in degree  (delta)
% rake      [scalar] rake in degrees (lambda)
% Mo        [scalar] scalar moment
% ttrans    [scalar] horizontal translation of pulse. Note this previously
%           had default value of 0.42 sec. 
% extras    [structure] for creating different src type
%               e.g. for ricker
%               extras.type = 'ricker'
%               extras.domfreq = 8;  % i.e. 8Hz
%
% OUTPUTS:
% dt        [scalar] computed time step to obey stability criterion
% stf       [matrix nX6]    
%               Row 1: [src_posy(1), src_posy(2), src_posy(3), 0,0,0]
%               Others: [M(1,1), M(2,2), M(3,3), M(1,3), M(2,3), M(1,2)]
%                           with M computed by make_mtensor

if nargin ==10
    ttrans = 0.42
end


if isempty(src_posy)
    src_posy = [14 12 22]
end

if nargin <7
    strike = 0;
    dip = 90;
    rake = 0;
    Mo = 1;
end  

if nargin ==11
    type = 'gaussian'
else
    type = extras.type;
    domfreq = extras.domfreq;
end

dt=.25*dx/max_vel; % calculate time step from stability criterion (wave travel max of 0.5eps box per time step)
disp(['the somputed time step is: ', num2str(dt)])

nt=tlength/dt   %compute total number of time steps		
t=(1:nt)*dt; 			% initialize time axis

% First attempt at a Ricker Wavelet (event1_3 and event2_3)
%src = (1-2*pi^2*f0^2*(t-t0).^2).*exp(-pi^2*f0^2*(t-t0).^2);
%src=diff(src);		
%src = (1-2*pi^2*f0^2*(t-t0).^2).*exp(-pi^2*f0.^2*(t-t0).^2);

switch type
    case 'gaussian'
        src = exp(-(t-ttrans).^2/(2*(pwidth)^2));
    case 'ricker'
        src = (1-2*pi^2*domfreq^2*(t-ttrans).^2).*exp(-pi^2*domfreq^2*(t-ttrans).^2);
end

M = make_mtensor(strike,dip,rake,Mo);
display('Moment Tensor: ')
display(M)
a = [src_posy(1) src_posy(2) src_posy(3) 0 0 0];  % initialise first row of stf.txt
for i = 1: length(t)
    % Note top is wrong because it is in the Kennett coordinates (for Olsen
    % must swap x and y
    %a = [a; M(1,1)*src(i), M(2,2)*src(i), M(3,3)*src(i), M(1,3)*src(i), M(2,3)*src(i), M(1,2)*src(i)];  %taxx(i,j),tayy(i,j),tazz(i,j),taxz(i,j),tayz(i,j),taxy(i,j)
    a = [a; M(2,2)*src(i), M(1,1)*src(i), M(3,3)*src(i), M(2,3)*src(i), M(1,3)*src(i), M(2,1)*src(i)];  %taxx(i,j),tayy(i,j),tazz(i,j),taxz(i,j),tayz(i,j),taxy(i,j)
end
save(fname,'a','-ascii')

stf = a;
figure
plot(t, a(2:end,6))
