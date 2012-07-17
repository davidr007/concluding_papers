function [circ_plus,circ_minus] = make_circle(xval, cent,r)
% This function returns the top (circ_plus) and bottom (circ_minus) halves
% of a circle with radius r and center cent.
%
% David Robinson
% 13 August 2008

circ_plus = cent(2) + sqrt(r^2-(xval-cent(1)).^2);
circ_minus = cent(2) - sqrt(r^2-(xval-cent(1)).^2);





