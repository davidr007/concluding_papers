function [P_CWIN_deltat] = vectorized_build_likelihood4_CWIN(delta_cwi,mu_1,sigma_1,mu_N,sigma_N)
% Builds the likelihood function: P(delta_CWIN|delta_t,I)
%
% INPUTS:
% delta_CWI [vector 1xm] m values at which to consider non-noisy CWI
%           values. Note that we integrate ocver these.
% mu_1      [vector 1xn] mean estimate of PDF (bounded Gaussian with delta>0)
%           for the CWI estimates defined for the n event pairs of interest
% sigma_1   [vector 1xn] standard deviation estimate of PDF (bounded Gaussian
%           with delta>0) for the CWI estimates defined for each of the n 
%           event pairs of interest.
% mu_N      [vector 1xn] noisy data value (i.e. delta_CWIN). Note that
%           typically this is given by the sample mean of all the noisy CWI
%           estimates that we have (Hence the terminology mu_N). There is
%           one of these for each of the n event pairs of interest.
% sigma_N   [vector 1xn] standard deviation of noisy data (i.e. delta_CWIN). Note that
%           typically this is given by the sample standard deviation of all the
%           noisy CWI estimates that we have (Hence the terminology mu_N). There is
%           one of these for each of the n event pairs of interest.
%
% David Robinson
% 16 June 2008

% Make sure delta_cwi is a column vector
delta_cwi = delta_cwi(:);

% Get everything into a large matrix size
%       number of rows = length(delta_cwi)
%       number of columns = length(delta_t) = length(mu_1) = length(mu_N)
MU_1 = repmat(mu_1, length(delta_cwi), 1);
SIGMA_1 = repmat(sigma_1, length(delta_cwi), 1);
DELTA_cwi = repmat(delta_cwi,1,length(mu_1));
SIGMA_N = repmat(sigma_N,length(delta_cwi),1);
MU_N = repmat(mu_N,length(delta_cwi),1);

% compute A
[theoretical_likelihood.Y,Phi_mu1_sigma1_0] = ...
    folded_normal_pdf(zeros(size(mu_1)),mu_1,sigma_1);
A = 1./((1-Phi_mu1_sigma1_0).*sigma_1*sqrt(2*pi));

% compute B
B = exp( -(DELTA_cwi- MU_1).^2 ./ (2*SIGMA_1.^2));

% Compute C
[rubish,Phi_mu2_sigma2_0] = ...
    folded_normal_pdf(0,DELTA_cwi,SIGMA_N);
C = 1./((1-Phi_mu2_sigma2_0).*SIGMA_N*sqrt(2*pi));

% Compute D
D = exp( -(MU_N-DELTA_cwi).^2 ./ (2.*SIGMA_N.^2));

integrand = B.*C.*D;
P_CWIN_deltat = A.* trapz(delta_cwi,integrand);


Lstar = -sum(log(P_CWIN_deltat));

%==================================================================


