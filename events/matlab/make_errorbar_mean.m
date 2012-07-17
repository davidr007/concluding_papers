function ydata = make_errorbar_mean(x, xdata)
% Defines the shape of the errorbar mean function (also the sigma function)
% Used by fit_errorbars
%
% Inputs
% x     [vector 3x1] defines the parameters of the function
% xdata [vector nx1] defines the values at which the fitting function is
%       evaluated
% 
% David Robinson
% 5 Feb 2008

% orig
ydata = x(1).*  (x(2)*xdata.^x(4)+x(3).*xdata.^x(5))./ ...
                (x(2)*xdata.^x(4)+x(3).*xdata.^x(5) + 1);


