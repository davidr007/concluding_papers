function [PDF, EV, CI68, CI95] = PDF4sep_between2GaussianEvents(mu1,sigma1,mu2,sigma2,sep_oi,diag_switch)
% This function is designed to compute the PDF for separation between two 
% events with locations that follow a Multi Variate Gaussian PDF
%
% Note that the function handles both the 2D and 3D cases depending on the
% length of input variables: mu1,sigma1, mu2, sigma2
%
% INPUTS: 
% mu1       [vectror] 
%               2D:  [x,y] location of mean for event1
%               3D:  [x,y,z] location of mean for event1
% sigma1    [vector]
%               2D:  [sigmax,sigmay] sigmas for event1
%               3D:  [sigmax,sigmay,sigmaz] sigmas for event1
% mu2       [vectror] 
%               2D:  [x,y] location of mean for event2
%               3D:  [x,y,z] location of mean for event2
% sigma2    [vector]
%               2D:  [sigmax,sigmay] sigmas for event2
%               3D:  [sigmax,sigmay,sigmaz] sigmas for event2
% sep_oi    [vector] contains the separation of interest for which the PDF
%           is returned
% diag_switch   [scalar]
%                   1 => darw diagnostics
%                   0 => do not draw diagnostics
%
%
% OUTPUTS: 
% PDF       [vector] PDF defined at the spearations given in sep_oi
% EV        [scalar] expected value for separation
% CI68      [vector] lower and upper bound of 68% confidence interval
% CI95      [vector] lower and upper bound of 95% confidence interval
%
% DEMO:
% 1) [PDF, EV, CI68, CI95] = PDF4sep_between2GaussianEvents([1 1],[0.5 0.5],[1 1],[.5 0.5],[0:0.01:6],1)
% 2) [PDF, EV, CI68, CI95] = PDF4sep_between2GaussianEvents([20 40],[10 10],[40 40],[10 10],[0:0.5:100],1)  
% 3) [PDF, EV, CI68, CI95] = PDF4sep_between2GaussianEvents([1 1 1],[0.5 0.5 0.5],[1 1 1],[0.5 0.5 0.5],[0:0.01:6],1)
% 4) [PDF, EV, CI68, CI95] = PDF4sep_between2GaussianEvents([50 50 50],[12 12 25],[100 50 50],[15 15 20],[0:10:800],1)
%
% David Robinson 
% 12 May 2008

% First we must work out if it is 2D or 3D
if length(mu1) == length(sigma1) & length(sigma1) == length(mu2) & length(mu2) == length(sigma2)
    nD = length(mu1); 
else
    error('lengths of mu1, sigma1, mu2 and sigma2 must be identical (length 2 for 2D or 3 for 3D)')
end

if nD ==2
    [PDF] = LOC_2DPDF(mu1,sigma1,mu2,sigma2,sep_oi);
elseif nD ==3
    [PDF] = LOC_3DPDF(mu1,sigma1,mu2,sigma2,sep_oi);
end

[EV,CI68, CI95,EV2] = LOC_getstats(mu1,mu2,sep_oi,PDF);

if diag_switch ==1
   figure
   plot(sep_oi,PDF)
   hold on
   % plot the expected value
   %h1 = plot(EV*[1 1], get(gca,'ylim'),'g');
   %plot(EV2*[1 1], get(gca,'ylim'),'g--');
   h2 = plot(CI68(1)*[1 1], get(gca,'ylim'),'r--');
   plot(CI68(2)*[1 1], get(gca,'ylim'),'r--');
   h3 = plot(CI95(1)*[1 1], get(gca,'ylim'),'c--');
   plot(CI95(2)*[1 1], get(gca,'ylim'),'c--');
   legend([h2,h3],{'68% CI','95% CI'})
end

%% ==================================================================
function [PDF] = LOC_2DPDF(mu1,sigma1,mu2,sigma2,sep_oi)
% local function for the 2D case
muX = mu1(1)-mu2(1);
muY = mu1(2)-mu2(2);
sigmaX = sqrt(sigma1(1)^2+sigma2(1)^2);
sigmaY = sqrt(sigma1(2)^2+sigma2(2)^2);

R = sep_oi;
theta = linspace(0, 2*pi, 200);
for i = 1:length(sep_oi)
   integrand = R(i)/(2*pi*sigmaX*sigmaY).*exp(-1/(2*sigmaX^2*sigmaY^2).*...
                  ( sigmaY^2.*(R(i).*cos(theta)-muX).^2  + sigmaX^2.*(R(i).*sin(theta)-muY).^2));
   PDF(i) = trapz(theta, integrand);     
end

% Do a quick re-normalisation in case the sampling is not good enough
% That if there are not enought theta 
% you could simply make more theta but this also increases computation
% A quick normalization here fixes the problem without causing any major issue 
PDF = PDF/trapz(sep_oi,PDF);

%% ==================================================================
function [PDF] = LOC_3DPDF(mu1,sigma1,mu2,sigma2,sep_oi)
% local function for the 2D case
muX = mu1(1)-mu2(1);
muY = mu1(2)-mu2(2);
muZ = mu1(3)-mu2(3);
sigmaX = sqrt(sigma1(1)^2+sigma2(1)^2);
sigmaY = sqrt(sigma1(2)^2+sigma2(2)^2);
sigmaZ = sqrt(sigma1(3)^2+sigma2(3)^2);

% initialise
R = sep_oi;
n_phi = 300;
n_theta = 150; 
phi = linspace(0, 2*pi, n_phi);
theta = linspace(0,pi,n_theta)';
PHI = repmat(phi,n_theta,1); 
THETA = repmat(theta, 1, n_phi);

for i = 1:length(sep_oi)
    A = -1/(2*sigmaX^2*sigmaY^2*sigmaZ^2);
    B = sigmaY^2*sigmaZ^2*(R(i).*sin(THETA).*cos(PHI) - muX).^2;
    C = sigmaX^2*sigmaZ^2*(R(i).*sin(THETA).*sin(PHI) - muY).^2;
    D = sigmaX^2*sigmaY^2*(R(i).*cos(THETA) - muZ).^2;
    integrand = 1/((2*pi)^(3/2)*sigmaX*sigmaY*sigmaZ) * exp(A.*(B+C+D)).*(R(i)).^2.*sin(THETA);
    tmp = trapz(theta, integrand); 
    PDF(i) = trapz(phi, tmp);     
end

% Do a quick re-normalisation in case the sampling is not good enough
% That iis if n_theta and n_phi is not large enough
% you could simply make n_theta and n_phi but this also increases computation
% A quick normalization here fixes the problem without causing any major issue 
PDF = PDF/trapz(sep_oi,PDF);


%% ==================================================================

function [EV,CI68, CI95,EV2] = LOC_getstats(mu1,mu2,sep_oi,PDF);
% get the important statistics
 

%EV = sqrt(sum( (mu1-mu2).^2 ) );
EV = trapz(sep_oi,sep_oi.*PDF)

CI68 = NaN*ones(1,2); 
CI95 = NaN*ones(1,2); 
CDF = cumtrapz(sep_oi,PDF);
ind68 = find(CDF>0.16 & CDF<1-0.16); 
ind95 = find(CDF>0.025 & CDF<1-0.025);

p1 = polyfit([sep_oi(ind68(1)-1), sep_oi(ind68(1))], [CDF(ind68(1)-1), CDF(ind68(1))], 1);
CI68(1) = (0.16-p1(2))/p1(1);
p2 = polyfit([sep_oi(ind68(end)-1), sep_oi(ind68(end))], [CDF(ind68(end)-1), CDF(ind68(end))], 1);
CI68(2) = (1-0.16-p2(2))/p2(1);
p3 = polyfit([sep_oi(ind95(1)-1), sep_oi(ind95(1))], [CDF(ind95(1)-1), CDF(ind95(1))], 1);
CI95(1) = (0.025-p3(2))/p3(1);
p4 = polyfit([sep_oi(ind95(end)-1), sep_oi(ind95(end))], [CDF(ind95(end)-1), CDF(ind95(end))], 1);
CI95(2) = (1-0.025-p4(2))/p4(1);

% CI68 = [sep_oi(ind68(1)), sep_oi(ind68(end))];
% CI95 = [sep_oi(ind95(1)), sep_oi(ind95(end))];

ind_EV = find(CDF>0.5); 
EV2 = sep_oi(ind_EV(1));




