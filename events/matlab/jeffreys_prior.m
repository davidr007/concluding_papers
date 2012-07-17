function [y] = jeffreys_prior(x)
% This function is designed to evaluate the jeffreys prior P(H|I)
% for the hypotheses H=x. For example, if we are using with the 
% CWI problem H=x=delta_t. 
%
% INPUTS:
% x     [vector 1xm]  containing the different hypotheses
%
% OUTPUTS:
% y     [vector 1xm] contains the jeffreys prior evaluated at each of the
%       hypothese
%
% David Robinson
% 7 February 2008

% check that size of x is correct
[n m] = size(x)
if n~=1|m<1
    error('x must be a row vector')
end

y = 1./x;


