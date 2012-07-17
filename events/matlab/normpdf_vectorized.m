function [y] = normpdf_vectorized(mu, sigma,x)
% Evaluates the m Gaussian probability density functions
% at n values each. 
%
% INPUTS:
% mu        [vector 1xm] the mean for m different Gaussian density
%           functions. 
% sigma     [vector 1xm] the sigma for m different Gaussian density
%           functions
% x         [matrix: nxm] Evaluates the m different Gaussian density 
%           functions at the values defined by each x. 
% ***   That is, for m =2  we evaluate the Gaussian PDF defined by mean = mu(2)
%       and sigma = sigma(2) at all the values defined in x(2,:). 
% 
% OUTPUTS:
% y         [matrix: mxn] matrix of densities following format described above. 
%
% NOTES: 
%   *   When using this for determining the likelihood functions with CWI estimates 
%       then m corresponds to m different hypothesis=delta_t and n corresponds
%       to n different data=delta_CWI
%
% David Robinson
% 7 Febuary 2008

%% First let's check the inputs
[n_mu, m_mu] = size(mu);
[n_sigma, m_sigma] = size(sigma);
[n_x, m_x] = size(x);
% check1 - lengths are correct
if m_mu~=m_sigma 
    error('Problem with input sizes: the number of rows in mu must equal the number in sigma')
elseif m_mu~=m_x
    error('Problem with input sizes: the number of rows in mu must equal the number in x')
end
% check2 - check shape
if n_mu~=1 | m_mu<1  % i.e. must be a row vector
    error('mu must be a column vector')
elseif n_sigma~=1 | m_sigma<1 % i.e. must be a rwo vector
    error('sigma must be a column vector')
end

%% copy mu and sigma vector so thjat they are the same size as x
MU = repmat(mu,n_x,1);
SIGMA = repmat(sigma,n_x,1); 


%% Evaluate the Gaussian PDF
y = exp(-0.5 * ((x - MU)./SIGMA).^2) ./ (sqrt(2*pi) .* SIGMA);


