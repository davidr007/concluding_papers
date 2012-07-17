function ydata = make_errorbar_mean_derivative(x, xdata)
% Computes the derivative of MU1 with respect to delta_t
% Note that the functional form is identical for SIGMA1
%
% Inputs
% x     [vector 5x1] defines the parameters of the function
% xdata [vector nx1] defines the values at which the fitting function is
%       evaluated
% 
% David Robinson
% 5 Feb 2008

% orig
ydata = x(1).* ...
                (   x(2)*x(4)*xdata.^(x(4)-1)+ ...
                    x(3)*x(5).*xdata.^(x(5)-1))./ ...
                (   x(2)*xdata.^x(4)+ x(3).*xdata.^x(5) +1).^2;


