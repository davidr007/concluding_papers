function [ynew_mu, ynew_sigma] = make_synthetic_muN(delta_t)

% The purpose of this function is to make the muN for the simple synthetic
% earthquake location problems....
%
% Note that it work in normalized space and that it assumes that mu_N is coming
% from the 1 to 5HZ all errorbar plot (i.e. the mean delta_cwi estimate)
%
% INPUTS:
% delta_t   [vector] true separations of interest. 

if nargin ==0
   delta_t = [0.05, 0.1, 0.15,0.2];  
end

if strcmp(filesep,'\')
    pathfordata = 'c:/datafiles/workstuff/sandpit/davidr/thesis_version2/diags/revised_2dacoustics/';
elseif strcmp(filesep,'/')
    pathfordata = '/export/storage/davidr/sandpit/davidr/thesis_version2/diags/revised_2dacoustics/';
end
load([pathfordata,'errorbar_data_pdftype_BoundedGaussian_col3_normtype1.mat']);

xorig = errorbar_data.known_sep_perfreq.for_1to5Hz;
yorig = errorbar_data.fitted_par_data.for_1to5Hz(:,1);  %the mean fits
yorig2 = errorbar_data.fitted_par_data.for_1to5Hz(:,2);  % the sigma fits

ynew_mu = spline(xorig,yorig,delta_t);
ynew_sigma = spline(xorig,yorig2,delta_t);


disp ('====================================')
for i = 1:length(delta_t)
    disp(['delta_t = ', num2str(delta_t(i)), ' : mu_N = ', num2str(ynew_mu(i)), ' : sigma_N = ', num2str(ynew_sigma(i))])
end
disp ('====================================')