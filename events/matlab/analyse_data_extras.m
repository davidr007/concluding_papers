
function [data, hf] = analyse_data_extras(medinp, atype,plot_struct,startval,diag,known_sep)

%INPUTS:
% medinp    [structure] multi-tiered structure containing all of the data for
%           the medium of interest i.e. separation vectors for all event
%           pairs at each station. Structure tiers are best explained by 
%           the following example for the medium exp001: 
%               exp001 = 
%                       xcorrStruc: [1x1 struct]
%                       sepStruc: [1x1 struct]
%               exp001.sepStruc =
%                       receiver000468_receiver000469: [1x1 struct]
%                       receiver000468_receiver000470: [1x1 struct]
%                       receiver000468_receiver000471: [1x1 struct]
%                       receiver000468_receiver000472: [1x1 struct]
%                       receiver000468_receiver000473: [1x1 struct]
%                       receiver000468_receiver000474: [1x1 struct]
%                       receiver000468_receiver000475: [1x1 struct]
%                       receiver000468_receiver000476: [1x1 struct]
%                       receiver000468_receiver000477: [1x1 struct]
%                       receiver000468_receiver000478: [1x1 struct]
%                       receiver000468_receiver000479: [1x1 struct]
%                       receiver000468_receiver000480: [1x1 struct]
%                       receiver000468_receiver000481: [1x1 struct]
%                       receiver000468_receiver000482: [1x1 struct]
%                       receiver000468_receiver000483: [1x1 struct]
%               exp001.sepStruc.receiver000468_receiver000469
%                       station1: [22x4 double]
%                       station2: [22x4 double]
%                       station3: [22x4 double]
%                       station4: [22x4 double]
%                       station5: [22x4 double]
%                       station6: [22x4 double]
%                       station7: [22x4 double]
%                       station8: [22x4 double]
%                       station9: [22x4 double]
%                       station10: [22x4 double]
%                       station11: [22x4 double]
% atype         [scalar - switch]
%                   1 => do analysis with separation estimates sepStruc
%                   2 => do analysis with maximum cross correlations xcorrStruc
% plot_struct   [structure with plotting options] - all are optional
%                  plot_struct.ylims => [vector 2x1] limits for y-axis 
%                  plot_struct.title => string to be appended to figure
%                                       specific title
%               If empty the ylimits are chose automatically
% startval      [scalar] startval used to discard intial separation estimates
%               for statistics and plots 
% diag      [   vector 1x5] switches for controlling diagnostics
%                   diag(1) == 1 => display key outputs to screen 
%                   diag(2) == 1 => plot errorbars mean +/- 1std 
%                   diag(3) == 1 => plot median +/1 1.483*MAD 
%                   diag(4) == 1 => plot boxplots
%                   diag(5) == 1 => construct plots showing mean/median and all
%                                   estimates as a function of window
%                                   centroid.
% known_sep [vector nx1] Vector containing known separation values. Note
%           that n must be the number of event pairs in eacn of the
%           shallow, mid and deep. 
%                   
% NOTE the are many processes within loops in this program that could be
% joined to reduce run time. However this would result in code that is more
% difficult to read (particularly given the multiple run options). A
% conscious decision has been made to keep the code more readable.


% Important hard wired parameters
MADscale = 1.483;

%Initialisations
%initialisation of vector for figure handles
hf.diag2 = 0; 
hf.diag3 = 0;
hf.diag4 = 0;
hf.diag5 = [];


if atype == 1  %determine structure for analysis (i.e. separation of maxXcorr)
    astruc = medinp.sepStruc;
elseif atype == 2
    astruc = medinp.xcorrStruc;
else
error('Invalid value for atype in ANALYSE_DATA')
end

% first we need to find out how many event pairs there are. At this point 
% we assume that the number of event pairs is equal for all three of the 
% depth ranges (i.e. shallow, middle, deep).
eventpairs = fieldnames(astruc);
nep = length(eventpairs);

% Now we can repeat the same for the station names and number of stations.
% Once again we assume that the stations are the same for all event-pairs
% and depth scenarios. 
stations = fieldnames(astruc.(eventpairs{1}));
nst = length(stations);


% here we put all of the data into a big matrix so that we can easily
% compute mean and median values
for i =1: nep  % loop over event pairs
    count =1;
    for j = 1: nst % loop over stations
        % Here we look at how many estimates there are for each
        % station for each of the shallow, mid and deep scnearios. This is
        % required for the counters to place data into a single matrix.
        swl = length(astruc.(eventpairs{i}).(stations{j})(startval:end,1));
        % Now we can put the data in the matrix
        data(count:count+swl-1,i) = astruc.(eventpairs{i}).(stations{j})(startval:end,2); 
        % Update the counter
        count = count+swl-1;
    end
end

% Intialisation to permute data1 and data2 values into single vector format for boxplots
known_sep_vec = [];
data_vec = [];
% 
%    
for i = 1:nep % loop over columns to create the mean and std values for each sim
    % this loop is required due to the presence of NaN in some entries of
    % data_shallow, data_med or data_deep   
    ind = ~isnan(data(:,i)); % find all vaules that are not NaN in column i
    mean_data(i) = mean(data(ind,i));  
    std_data(i) = std(data(ind,i));
    median_data(i) = median(data(ind,i)); 
    MAD_data(i) = median(abs(data(ind,i)-median_data(i)));
    % Here we permute the data1 and data2 values into single vector format for boxplots
    known_sep_vec = [known_sep_vec; known_sep(i)*ones(sum(ind),1)];
    data_vec = [data_vec; data(ind,i)];   
    clear ind 
end
% 
% Do plot of mean and +/- 1 std
if diag(2) ==1
    hf.diag2 = figure;
    plot(known_sep, known_sep,'k','linewidth',2)
    hold on
    for i = 1: length(mean_data)
        h = plot(known_sep(i)-5, mean_data(i),'go');
        errorbar(known_sep(i)-5, mean_data(i), ...
            std_data(i),'g')% 
    end
    % cd diags
    xlabel('known separation (m)')
    ylabel('computed separation (m)')
    try, title(['mean and std: ',plot_struct.title]), 
    catch, title('mean and std'), end
    %set(gca,'xticklabel',round(known_sep))
    try, set(gca,'ylim', plot_struct.ylimits), catch, end
    try, set(gca,'xlim',plot_struct.xlimits), catch, end
    % print(h1, '-depsc',['mean_',fname(1:end-4),'.eps'])   
end

% Do plot of median and +/- 1.483MAD
if diag(3) ==1
   error('Functionality not yet available in analyse_data_extras')
end

% Do boxplots
if diag(4) ==1
    error('Functionality not yet available in analyse_data_extras')
end  

% Do plots of mean and median and all traces
% Note that there are some similar lines to those above. We keep it all
% separate here so the code is easier to read.
if diag(5) ==1 
  error('Functionality not yet available in analyse_data_extras')
end
