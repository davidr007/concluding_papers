function [P_CWIN_deltat] = build_likelihood4_CWIN(delta_t,delta_cwi,mu_1,sigma_1,mu_N,sigma_N)
% Builds the likelihood function: P(delta_CWIN|delta_t,I)
%
% INPUTS: 
% delta_t   [vector 1xn] n values at which to evaluate P(delta_CWIN|delta_t,I)
% delta_CWI [vector 1xm] m values at which to consider non-noisy CWI
%           values. Note that we integrate ocver these. 
% mu_1      [vector 1xn] mean estimate of PDF (bounded Gaussian with delta>0) 
%           for the CWI estimates defined for each value in delta_t. 
% sigma_1   [vector 1xn] standard deviation estimate of PDF (bounded Gaussian
%           with delta>0) for the CWI estimates defined for each value in delta_t. 
% mu_N      [scalar] noisy data value (i.e. delta_CWIN). Note that
%           typically this is given by the sample mean of all the noisy CWI
%           estimates that we have (Hence the terminology mu_N) 
% sigma_N   [scalar] standard deviation of noisy data (i.e. delta_CWIN). Note that
%           typically this is given by the sample standard deviation of all the 
%           noisy CWI estimates that we have (Hence the terminology mu_N) 
%
% David Robinson
% 7 March 2008

% Make sure delta_cwi is a column vector
delta_cwi = delta_cwi(:);




%for i = 1:length(delta_t)   % loop over all of the delta_t
    % compute A
    [theoretical_likelihood.Y,Phi_mu1_sigma1_0] = ...
        folded_normal_pdf(delta_t,mu_1,sigma_1);
    A = 1./((1-Phi_mu1_sigma1_0).*sigma_1*sqrt(2*pi)); 
    
    % Get everything into a large matrix size 
    %       number of rows = length(delta_cwi)
    %       number of columns = length(mu_1)
    MU_1 = repmat(mu_1, length(delta_cwi), 1);
    SIGMA_1 = repmat(sigma_1, length(delta_cwi), 1);
    DELTA_cwi = repmat(delta_cwi,1,length(mu_1));
    
    % compute B
    B = exp( -(DELTA_cwi- MU_1).^2 ./ (2*SIGMA_1.^2));
    
    % Compute C
%        [rubish,Phi_mu2_sigma2_0] = ...
%         folded_normal_pdf(0,delta_cwi,sigma_N);
        [rubish,Phi_mu2_sigma2_0] = ...
        folded_normal_pdf(0,mu_N,sigma_N);
       C = 1./((1-Phi_mu2_sigma2_0).*sigma_N*sqrt(2*pi)); % This is a scalar
       %C = repmat(c,1,length(mu_1));
    
        % Compute D
       D = exp( -(mu_N-DELTA_cwi).^2 ./ (2.*sigma_N ^2));
            
    integrand = B.*D;
    P_CWIN_deltat = A.*C.* trapz(delta_cwi,integrand);
%end