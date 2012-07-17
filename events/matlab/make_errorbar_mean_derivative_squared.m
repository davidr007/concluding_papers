function ydata = make_errorbar_mean_derivative_squared(x, xdata)
% Defines the shape of the errorbar mean function (also the sigma function)
% Used by fit_errorbars
%
% Inputs
% x     [vector 5x1] defines the parameters of the function
% xdata [vector nx1] defines the values at which the fitting function is
%       evaluated
% 
% David Robinson
% 5 Feb 2008

%%xdata(xdata==0) = 0.001;

% orig
ydata = x(1).* ...
                (   x(2)*x(4)/2*xdata.^(x(4)/2-1)+ ...
                    x(3)*x(5)/2.*xdata.^(x(5)/2-1))./ ...
                (   x(2)*xdata.^(x(4)/2)+ x(3).*xdata.^(x(5)/2) +1).^2;

%ydata(isinf(ydata))=0;
