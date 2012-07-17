function taper = ztaper_tsignal(veclength, tlength, vecmax)

% INPUTS:
% veclength  [scalar] length of vector to be tapered. Note that the taper
%            will have this length.
% tlength    [scalar] length for tapering on each side of the taper i.e 
%            tlength << veclength
% vecmax     [scalar] Optional - maximum value of taper - should be the same as the
%               maximum value of the vector to be tapered.
%
% OUTPUTS:
% taper      [vector] contains 1 in the centre and the ramps down to zero 
%            at each end 
%
%
% Note that vecmax is optional. If the amplitude of the waveform is
% important than vecmax should be set to 1. Otherwise setting vecmax to the
% maximum of the waveform means that the taper will plot nicely over the
% waveform. 



if nargin ==2
    vecmax = 1
end

lof_constz = floor(0.2*tlength);
lof_ramp = tlength - lof_constz;
lof_center = veclength - 2*(lof_constz+lof_ramp);

x = [0:lof_ramp-1];
b = pi/(2*x(end));
ramp = sin(b*x);

ans = iseven(veclength);
taper = [zeros(1, lof_constz), vecmax*ramp, vecmax*ones(1,lof_center), vecmax*fliplr(ramp), zeros(1,lof_constz)]';
% if ans ==1 % length of vector = length of taper is even
%     taper = [zeros(1, lof_constz), vecmax*ramp, vecmax*ones(1,lof_center), vecmax*fliplr(ramp), zeros(1,lof_constz)]';
% elseif ans ==0 % length of vector = length of taper is odd
%         taper = [zeros(1, lof_constz), vecmax*ramp, vecmax*ones(1,lof_center), vecmax*fliplr(ramp), zeros(1,lof_constz)]';
% else
%     error(['ERROR: incorrect value for veclength - is should odd or even (i.e. integer)'])
% end
