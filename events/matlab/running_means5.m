function y = running_means5(x)
%
% Robinson, D., Dhu, T. and Schneider, J. 
% "SUA: A Computer Program to Compute Regolith Site-response and Estimate Uncertainty for
% Probabilistic Seismic Hazard Analyses".
% Computers and Geosciences
%
% y = RunningMeans5(x)
%
% RunningMeans5 computes the the five point running mean of x 
% i.e. it smooths x. 
%
% Inputs:
% x         Vector to be smoothed
%
% Outputs
% y         Smoothed Vector
%
%                                           David Robinson
%                                           Last Updated: 15/07/02

x = x(:);

len = length(x);
x_n2=x(1:len-4);  x_n1=x(2:len-3);  x_p2=x(5:len);  x_p1=x(4:len-1);
y = [x(1); (x(1)+x(2)+x(3))/3; (x_n2+x_n1+x(3:len-2)+x_p1+x_p2)/5; (x(len-2)+x(len-1)+x(len))/3;  x(len)];
