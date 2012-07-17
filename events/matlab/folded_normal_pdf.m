function [Y,PL0] = folded_normal_pdf(X,mu,sigma)
% Computes the value of the folded normal PDF given 
% by mean (mu) and standard deviation (sigma). 
%
% Note this is defined by creating a normal distribution and 
% folding those with negative X back into the positive. 
%
% To compute the folded normal PDF we simply consider the Normal PDF 
% P(X=x) = f(mu,sigma,x) for X>0 and multiple it by the the recipricol of 
% the cummulative Normal distribution function at X=0 i.e. 1/P(X<0).
%
% INPUTS: 
%
% OUTPUTS: 
%
% David Robinson 
% .... 

% First check that all values in X >= 0
ind = find(X<0);
if ~isempty(ind)
    error('ERROR in folded_normal_pdf.m: All values in input variable X must be >=0 ')
end

% Compute P(X<0)

%zerovec = zeros(size(mu));
%PL0= normcdf(zerovec,mu,sigma);
PL0= normcdf(0,mu,sigma);

Y = normpdf(X,mu,sigma)./(1-PL0);


