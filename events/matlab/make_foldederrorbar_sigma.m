function ydata = make_foldederrorbar_sigma(x, xdata)
% Defines the shape of the errorbar sigma function 
% Used by fit_errorbars
%
% Inputs
% x     [vector 4x1] defines the parameters of the function
% xdata [vector nx1] defines the values at which the fitting function is
%       evaluated
% 
% David Robinson
% 5 Feb 2008

% orig
% ydata = x(1).*  (x(2)*xdata.^x(4)+x(3).*xdata)./ ...
%                 (x(2)*xdata.^x(4)+x(3).*xdata + 1);


%ydata = 0.02 + x(1).*  (x(2)*xdata.^x(4)+x(3).*xdata.^x(5))./ ...
%                 (x(2)*xdata.^x(4)+x(3).*xdata+x(6));

% This was the latest try at making them different
           ydata =  0.017 + ...%x(5).*xdata + ... %0.0167+...
               x(1).*  (x(2)*xdata.^x(4)+x(3).*xdata.^x(5))./ ...
               (x(2)*xdata.^x(4)+x(3).*xdata.^x(5) + 1);
            

           
%        ydata = abs(x(1).*  (x(2)*xdata.^x(4)+x(3).*xdata)./ ...
%                 (x(2)*xdata.^x(4)+x(3).*xdata + 1));     
%             
            
         %   ydata = x(1)+x(2)*exp(-x(3)*xdata);
