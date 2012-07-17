function M = make_mtensor(strike,dip,rake,Mo)
% makes the monent tensor (not the time dependant part).
%
% INPUTS:
% strike    [scalar] strike in degrees (theta)
% dip       [scalar] dip in degree  (delta)
% rake      [scalar] rake in degrees (lambda)
% Mo        [scalar] scalar moment
%
% OUTPUT:
% M         [matrix 3x3] moment tensor (without the source time function)
%
% CAUTION: Moment tensor elements are limited in precision to 10^10. 
%
% NOTE:
%   x is positive North
%   y is positive East
%   z is positive Down


% convert to randians and change variable to name for brevity;
t = strike*pi/180; %theta
d = dip*pi/180; % delta
l = rake*pi/180; % lambda

%initialise M
M = zeros(3,3);

M(1,1) = -Mo*(sin(d)*cos(l)*sin(2*t) + sin(2*d)*sin(l)*(sin(t))^2);

M(1,2) = Mo*(sin(d)*cos(l)*cos(2*t)+sin(2*d)*sin(l)*sin(2*t)/2);
M(2,1) = M(1,2);

M(1,3) = -Mo*(cos(d)*cos(l)*cos(t) + cos(2*d)*sin(l)*sin(t));
M(3,1) = M(1,3);

M(2,2) = Mo*(sin(d)*cos(l)*sin(2*t) - sin(2*d)*sin(l)*(cos(t))^2);

M(2,3) = -Mo*(cos(d)*cos(l)*sin(t) - cos(2*d)*sin(l)*cos(t));
M(3,2) = M(2,3);

M(3,3) = Mo*sin(2*d)*sin(l);

M = round(M*10^10)/10^10;


