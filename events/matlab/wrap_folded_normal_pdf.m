function [Y] = wrap_folded_normal_pdf(x,xdata)
% This is a wrapper so we can use folded_normal_pdf in 
% lsqcurvefit. 
%
% Note that it returns the the PDF of a folded normal distribution
%
% INPUTS: 
% x         x(1) = mean of normal pdf which is folded
%           x(2) = std of normal pdf which is folded
%           x(3) vertical scale parameter
% xdata     [vector] values at which to determine the folded normal 
%           pdf

[Y,PLO] = folded_normal_pdf(xdata,x(1),x(2));
Y = x(3)*Y;