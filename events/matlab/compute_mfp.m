function [mfp,maxp_hom,maxp_grm] = compute_mfp(fit_lines, maxp_hom, maxp_grm,plot_struct)
%
% This function is used to generate the mean free path.
% The user must supply intensity-distance data for propagation
% in homogeneous (maxp_hom) and heterogeneous (e.g. maxp_grm)
% media. Note that intensity = max(abs(amplitude))^2.
% A scatter plot of distance versus log(maxp_grm./maxp_hom) is 
% created. Straight lines are fitted to the to the scatter using
% (1) least squares and (2) weighted least squares with weights
% defined by 1/distance. The mean free path is defined as -1/slope. 
%
% INPUTS:
% fit_lines     [scalar]
%                   0 => do not fit lines
%                   1 => fit weighted least squares line
%                   2 => fit least squares
%                   3 = > fit 1 and 2
% maxp_hom      [char or matrix] homogeneous distance-intensity data
%                   char   =>   pathname to file containing the distnace 
%                               (col 1) and intensity (col 2) data. 
%                   matrix =>   2 column matrix stored in workspace
% maxp_hom      [char or matrix] heterogeneous distance-intesnity data
%               with the same format as maxp_hom.
% plot_struct   [structure] containing obtional plotting options e.g.
%                   plot_struct.title = title for plot
%                   plot_struct.xlabel = xlabel
%                   plot_struct.ylabel = ylabel
%                   plot_struct.ylim = ylim
%                   plot_struct.xlim = xlim
%                   plot_struct.fsize = fontsize for all labels.
%
%
% OUTPUTS: 
% mfp           [vector 2x1]
%               mfp(1) =    mean free path determined by weighted least
%                           squares
%               mfp(2) =    mean free path determined by least squares
%
% David Robinson
% 12 April 2006

mfp = [NaN NaN];  % intialise output

if ischar(maxp_hom)
    maxp_hom = load(maxp_hom);
end

if ischar(maxp_grm)
    maxp_grm = load(maxp_grm);
end

%check that maxp_hom and maxp_grm are for the same points
tmp = maxp_hom(:,1) - maxp_grm(:,1); 
if min(tmp)~=0 & max(tmp)~=0
    error('maxp_hom and maxp_grm appear to be for different points (distances)')
end

% Remove any NaN values
x =  maxp_hom(:,1); 
y = log(maxp_grm(:,2)./maxp_hom(:,2)); 
x2 = x(~isnan(x));
y2 = y(~isnan(x));
x3 = x2(~isnan(y2));
y3 = y2(~isnan(y2));
% Remove zero values for distance
x4 = x3(x3~=0);
y4 = y3(x3~=0);

% Plot sctatter
figure
plot(x3, y3 ,'.')
hold on
% Try  and cath to reflect optional nature of plotting options
try, title(plot_struct.title,'fontsize',plot_struct.fsize),catch,  end
try, xlabel(plot_struct.xlabel,'fontsize',plot_struct.fsize),catch,end
try, ylabel(plot_struct.ylabel,'fontsize',plot_struct.fsize), catch, end
try, set(gca, 'ylim',plot_struct.ylim), catch, end
try, set(gca, 'xlim',plot_struct.xlim), catch, end

xlimits = get(gca,'xlim');

% fit weighted least squares
if fit_lines == 1 | fit_lines == 3
    %w = ones(1,length(x4)); % for testing against normal least squares
    w = 1./x4;  % define weights
    det = sum(w.^2)*sum(w.^2.*x4.^2) - (sum(w.^2.*x4))^2; %compute determinant
    yint = 1/det * (sum(w.^2.*x4.^2)*sum(w.^2.*y4)-sum(w.^2.*x4)*sum(w.^2.*x4.*y4)); % compute y-intercept
    slope = 1/det * (-sum(w.^2.*x4)*sum(w.^2.*y4) + sum(w.^2)*sum(w.^2.*x4.*y4));  
    y_fl1 = yint + slope*[xlimits(1),xlimits(2)];
    h1 = plot([xlimits(1),xlimits(2)],y_fl1,'r','linewidth',2); 
    mfp(1) = -1/slope;
end

% fit least squares (not weighted)
if fit_lines == 2 | fit_lines == 3
    p = polyfit(x4,y4,1);
    y_fl2 = polyval(p,[xlimits(1),xlimits(2)]);
    h2 = plot([xlimits(1),xlimits(2)],y_fl2,'y','linewidth',2);  
    mfp(2) = -1/p(1);
end