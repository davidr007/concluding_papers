function [loc_out] = convert_eventcoordinates(E1,E2,E3,convert_type, loc_in)
% The purpose of this function is to convert the event coordinates between
% cartesian (as in original locations) with unit vector i,j and k 
% to and from the local coordinate system with l,m,n. 
% 
% The transformation is defined by three points in the original cartesian
% coordinates 
%   E1 = x1i + y1j + z1k
%   E2 = x2i + y2j + z2k
%   E3 = x3i + y3j + z3k
%
% Which are converted to a new coordinate system such that 
% n is normal to the plane going through E1, E2 and E3
% l is in direction defined by E1 to E2
% m is perpendicular to l in the plane going through E1, E2 and E3
%
% Note that in the new coordinate system we have 
% E1 = ol + 0m + 0n
% E2 = dE1E2l + 0m + 0n
% E3 = ????


% Example 1: - cartesian input E1 already at origin - no translation
% [loc_out] = convert_eventcoordinates([0 0 0],[ 1 2 3],[ -2 3 4], 'cart2loc', [0 0 0;1 2 3; -2 3 4; 5 6 7])
% [loc_out] = convert_eventcoordinates([0 0 0],[ 1 2 3],[ -2 3 4], 'loc2cart', loc_out)
%
% Example 2: - cartesian input E1 already at origin - no translation
% [loc_out] = convert_eventcoordinates([1 -2 3],[ 1 2 3],[ -2 3 4], 'cart2loc', [1 -2 3;1 2 3; -2 3 4; 5 6 7])
% [loc_out] = convert_eventcoordinates([1 -2 3],[ 1 2 3],[ -2 3 4], 'loc2cart', loc_out)
%

if nargin ==0
    E1 = [0 0 0];
    E2 = [ 1 1 1 ];
    E3 = [ 5 1 0]; 
    convert_type= 'cart2loc';
    loc_in = [E1;E2;E3];
end

[n1 m1] = size(loc_in);
translation = E1(:);
x1 = E1(1)-translation(1);
x2 = E2(1)-translation(1);
x3 = E3(1)-translation(1);
y1 = E1(2)-translation(2);
y2 = E2(2)-translation(2);
y3 = E3(2)-translation(2);
z1 = E1(3)-translation(3);
z2 = E2(3)-translation(3);
z3 = E3(3)-translation(3);
switch convert_type
    case 'cart2loc' % If in cartesian coordinates may have to shift the events so that the first event is at the origin
        loc_in = loc_in - repmat(translation(:)',n1,1);
end
        
% d12 = sqrt( (x2-x1)^2 + (y2-y1)^2 +(z2-z1)^2);
% d13 = sqrt( (x3-x1)^2 + (y3-y1)^2 +(z3-z1)^2);
% d23 = sqrt( (x3-x2)^2 + (y3-y2)^2 +(z3-z2)^2);
% s1 = (d12^2+d13^2-d23^2)/(2*d12);
% s2 = sqrt(d13^2-s1^2);
% 
% B = [x1 y1 z1; x2 y2 z2; x3 y3 z3];
% Row1_t = inv(B)*[0;d12;s1];
% Row2_t = inv(B)*[0;0;s2];
% Row3_t = inv(B)*[0;0;0];
%
%R = [Row1_t';Row2_t';Row3_t']



%% First we define the new basis vectors in terms of the original coordinate system i,j,k...

% l is defined by the vector pointing along E1E2
d12 = sqrt( (x2-x1)^2 + (y2-y1)^2 +(z2-z1)^2);
l = 1/d12 * [x2-x1,y2-y1 , z2-z1];

% n is defined by the normal to the plane passing through E1,E2,E3
A = y1*(z2 - z3) + y2*(z3 - z1) + y3*(z1 - z2);
B = z1*(x2 - x3) + z2*(x3 - x1) + z3*(x1 - x2);
C = x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2);
D = -(x1*(y2*z3 - y3*z2) + x2*(y3*z1 - y1*z3) + x3*(y1*z2 - y2*z1));
n = [A,B,C];
n = 1/sqrt(sum(n.^2)) * n;

% m is the mutual cross product between l and n
m = [ l(2)*n(3) - l(3)*n(2), l(3)*n(1) - l(1)*n(3), l(1)*n(2) - l(2)*n(1)];
m = m;


A = [l',m', n'];



loc_out = NaN*ones(size(loc_in));
switch convert_type
    case 'cart2loc'

        for i = 1: n1
            loc_in_tmp = loc_in(i,:)';
            loc_out_tmp = inv(A)*loc_in_tmp;
            loc_out(i,:) = loc_out_tmp';
        end

    case 'loc2cart'
        for i = 1: n1
            loc_in_tmp = loc_in(i,:)';
            loc_out_tmp = A*loc_in_tmp;
            loc_out(i,:) = loc_out_tmp';
        end
        loc_out = loc_out + repmat(translation(:)',n1,1);
end


