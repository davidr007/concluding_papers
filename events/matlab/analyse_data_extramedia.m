
function [data_eg1, mean_eg1, std_eg1, domfreq_mean4rpair,hf] = analyse_data(medinp, atype,plot_struct,startval,diag,known_sep, col, normtype,existing_stats_asft,pdftype)

%INPUTS:
% medinp    [structure] multi-tiered structure containing all of the data for
%           the medium of interest i.e. separation vectors for all event
%           pairs at each station for multiple depths (e.g. shallow, middle
%           deep). Structure tiers are best explained by the following
%           example for the medium grm004:
%               grm004 =
%                       xcorrStruc: [1x1 struct]
%                       sepStruc: [1x1 struct]
%               grm004.sepStruc =
%                       mid: [1x1 struct]
%                       deep: [1x1 struct]
%                       shallow: [1x1 struct]
%               grm004.sepStruc.mid =
%                       receiver000036_receiver000037: [1x1 struct]
%                       receiver000036_receiver000038: [1x1 struct]
%                       receiver000036_receiver000039: [1x1 struct]
%                       receiver000036_receiver000040: [1x1 struct]
%                       receiver000036_receiver000076: [1x1 struct]
%                       receiver000036_receiver000077: [1x1 struct]
%                       receiver000036_receiver000078: [1x1 struct]
%                       receiver000036_receiver000041: [1x1 struct]
%                       receiver000036_receiver000082: [1x1 struct]
%                       receiver000036_receiver000083: [1x1 struct]
%                       receiver000036_receiver000084: [1x1 struct]
%                       receiver000036_receiver000042: [1x1 struct]
%                       receiver000036_receiver000043: [1x1 struct]
%                       receiver000036_receiver000044: [1x1 struct]
%                       receiver000036_receiver000045: [1x1 struct]
%           grm004.sepStruc.mid.receiver000036_receiver000037   %
%                       station1: [12x4 double]
%                       station2: [12x4 double]
%                       station3: [12x4 double]
%                       station4: [12x4 double]
%                       station5: [12x4 double]
%                       station6: [12x4 double]
%                       station7: [12x4 double]
%                       station8: [12x4 double]
%                       station9: [12x4 double]
%                       station10: [12x4 double]
%                       station11: [12x4 double]
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
%                   diag(6) ==1 =>  construct plots showing mean Rmax and all
%                                   estimates as a function of sliding window
%                                   centroid.
% known_sep [vector nx1] Vector containing known separation values. Note
%           that n must be the number of event pairs in eacn of the
%           shallow, mid and deep.
% col       [OPTIONAL - scalar]
%               []  => col = 2 - use Taylor series approximation
%               2 - Use Taylor Series approximation
%               3 - Use autucorrelation of unperturbed wave
%               4 - Use autucorrelation of perturbed wave
% normtype  [Optional - Scalar]
%               [] => normtype = 0 - no normalisation
%               0 - No normalisation of separation_true
%               1 - Normalise separation_true (i.e. multiply by domfreq/velocity)
% existing_stats_asft [optional - structure]
%               if existing_stats_asft is given and if diag(5) ==1 the
%               sliding time window mean CWI estimates from an earlier
%               analysis are also added to the plot.
% pdftype       [string] defines the type of pdf to use when estimating the
%               central tendancy and uncertainty
%               'Gaussian' => use a Gaussian PDF
%               'BoundedGaussian' => Use a Bounded Gaussian
%
% NOTE the are many processes within loops in this program that could be
% joined to reduce run time. However this would result in code that is more
% difficult to read (particularly given the multiple run options). A
% conscious decision has been made to keep the code more readable.

vs = 3300;

if ~exist('col')
    col =2
    disp(['Setting col = 2 - i.e. using Taylor series approx.' ])
elseif isempty(col)
    col=2;
    disp(['Setting col = 2 - i.e. using Taylor series approx.' ])
elseif col == 2 | col==3 | col==4
    % no need to do anything
else
    disp(['col = ', num2str(col)])
    error('ERROR: Invalid value of col in analyse_data')
end

if ~exist('normtype')
    normtype =0
    disp(['Setting normtype = 0 - i.e. no normalisation of separation_true' ])
elseif isempty(normtype)
    normtype=0;
    disp(['Setting normtype = 0 - i.e. no normalisation of separation_true' ])
elseif normtype == 0 | normtype==1
    % no need to do anything
else
    disp(['normtype = ', num2str(normtype)])
    error('ERROR: Invalid value of normtype in ANALYSE_DATA')
end

% Bounds ofr use when pdftype = 'BoundedGaussian'
if normtype ==1
    lbounds_start = [0 0 0];
    upperbounds_start = [1 1 10000];
elseif normtype ==0
    lbounds_start = [0 0 0];
    upperbounds_start = [2000 2000 10000];
end


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
eventpairs_eg1 = fieldnames(astruc.eg1);
nep = length(eventpairs_eg1);

% Now we can repeat the same for the station names and number of stations.
% Once again we assume that the stations are the same for all event-pairs
% and depth scenarios.
stations = fieldnames(astruc.eg1.(eventpairs_eg1{1}));
nst = length(stations);


% here we put all of the data into a big matrix so that we can easily
% compute mean and median values
for i =1: nep  % loop over event pairs
    count_eg1 = 1;
    for j = 1: nst % loop over stations
        % Here we look at how many estimates there are for each
        % station for each of the shallow, mid and deep scnearios. This is
        % required for the counters to place data into a single matrix.
        swl_eg1 = length(astruc.eg1.(eventpairs_eg1{i}).(stations{j})(startval.eg1:end,1));
        % Now we can put the data in the matrix
        data_eg1(count_eg1:count_eg1+swl_eg1-1,i) = astruc.eg1.(eventpairs_eg1{i}).(stations{j})(startval.eg1:end,col);
        % Compute the dominant frequencies for each station (note all windows considered separately due to vectors)

        domfreq_eg1(j) = mean( ...
            [medinp.domfreq.eg1.(['ref_',eventpairs_eg1{1}(1:14)]).(stations{j}), ...
            medinp.domfreq.eg1.(eventpairs_eg1{i}(16:end)).(stations{j})]);


        % Normalise the known separations if required
        if normtype ==1
            data_eg1(count_eg1:count_eg1+swl_eg1-1,i) = data_eg1(count_eg1:count_eg1+swl_eg1-1,i)*domfreq_eg1(j)/vs;
        end
        % Update the counters
        count_eg1 = count_eg1+swl_eg1-1;

    end  % end of loop over stations

    % compute average dominant frequency for the receiver pair (i.e.
    % average over all stations) -note we join shallow, mid and deep since
    % we only use on value for the plotting.
    domfreq_mean4rpair(i) = mean([domfreq_eg1]);


end % end of loop over receiver pairs

% Intialisation to permute data1 and data2 values into single vector format for boxplots
known_sep_vec_eg1 = [];
data_vec_eg1 = [];
%
%
for i = 1:nep % loop over columns to create the mean and std values for each sim
    % this loop is required due to the presence of NaN in some entries of
    % data_shallow, data_med or data_deep
    ind = ~isnan(data_eg1(:,i)); % find all vaules that are not NaN in column i
    mean_eg1_tmp(i) = mean(data_eg1(ind,i));
    std_eg1_tmp(i) = std(data_eg1(ind,i));
    if strcmp(pdftype,'Gaussian')
        mean_eg1(i) = mean_eg1_tmp(i);
        std_eg1(i) = std_eg1_tmp(i);
    elseif strcmp(pdftype,'BoundedGaussian')
        [N X] = hist(data_eg1(ind,i),10);
        lbounds = [min(X), 0, 0];
        upperbounds = [max([upperbounds_start, max(X)]), max(X)-min(X), 10000];
        [out1,resnorm,residual,exitflag,output] = lsqcurvefit(@wrap_folded_normal_pdf,[mean_eg1_tmp(i),std_eg1_tmp(i),1],X,N,lbounds, upperbounds);
        mean_eg1(i) = out1(1);
        std_eg1(i) = out1(2);
    end
    median_eg1(i) = median(data_eg1(ind,i));
    MAD_eg1(i) = median(abs(data_eg1(ind,i)-median_eg1(i)));
    % Here we permute the data1 and data2 values into single vector format for boxplots
    known_sep_vec_eg1 = [known_sep_vec_eg1; known_sep(i)*ones(sum(ind),1)];
    data_vec_eg1 = [data_vec_eg1; data_eg1(ind,i)];
    clear ind
end
%end
%
% Do plot of mean and +/- 1 std
if diag(2) ==1
    hf.diag2 = figure;
    try, title(['mean and std: ',plot_struct.title]),
    catch, title('mean and std'), end

    if normtype ==0
        known_sep_tmp = known_sep;
        xlabel('known separation (m)')
        ylabel('computed separation (m)')
    elseif normtype ==1
        known_sep_tmp = known_sep.*domfreq_mean4rpair/vs;
        xlabel('dimensionless known separation')
        ylabel('dimensionless computed separation')
        try, plot_struct.ylimits = plot_struct.ylimits*domfreq_mean4rpair(1)/vs; catch, end
        try, plot_struct.xlimits = plot_struct.xlimits*domfreq_mean4rpair(1)/vs; catch, end
    end
    %plot(known_sep_tmp, known_sep_tmp,'k','linewidth',2)
    %hold on
    for i = 1: length(mean_eg1)
        h_eg1 = plot(known_sep_tmp(i), mean_eg1(i),'ko'); % h-translation of 10% of known_sep(1) ~=5m
        hold on
        errorbar(known_sep_tmp(i), mean_eg1(i), ...
            std_eg1(i),'k')%
    end

    %set(gca,'xticklabel',round(known_sep))
    try, set(gca,'ylim', plot_struct.ylimits), catch, end
    try, set(gca,'xlim',plot_struct.xlimits), catch, end
end


% Do plot of median and +/- 1.483MAD
if diag(3) ==1
    hf.diag3 = figure;
    plot(known_sep, known_sep,'k','linewidth',2)
    hold on
    for i = 1: length(mean_eg1)
        h_eg1 = plot(known_sep(i)-5, median_eg1(i),'go');
        errorbar(known_sep(i)-5, median_eg1(i), ...
            1.483*MAD_eg1(i),'g')%

    end
    % cd diags
    xlabel('known separation (m)')
    ylabel('computed separation (m)')
    try, title(['median and 1.483*MAD: ',plot_struct.title]),
    catch, title('median and std'), end
    %set(gca,'xticklabel',round(known_sep))
    try, set(gca,'ylim', plot_struct.ylimits), catch, end
    try, set(gca,'xlim',plot_struct.xlimits), catch, end
    % print(h1, '-depsc',['mean_',fname(1:end-4),'.eps'])
end

% Do boxplots
if diag(4) ==1
    hf.diag4 = figure;
    plot(known_sep, known_sep,'k','linewidth',2)
    hold on
    h_eg1 = boxplot(data_vec_eg1,known_sep_vec_eg1,'position',known_sep-7,'widths', 5,'colors','ggg','symbol','g+');
    xlabel('known separation (m)')
    ylabel('computed separation (m)')
    try, title(['Boxplots: ',plot_struct.title]),
    catch, title(['Boxplots']), end
    try, set(gca,'ylim', plot_struct.ylimits), catch, end
    try, set(gca,'xlim',plot_struct.xlimits), catch, end
    xlimits = get(gca,'xlim')
    set(gca,'xtick',linspace(xlimits(1),xlimits(2),9),'xticklabels',linspace(xlimits(1),xlimits(2),9))
end

% Do plots of mean and median and all traces
% Note that there are some similar lines to those above. We keep it all
% separate here so the code is easier to read.
% if diag(5) ==1
%     % Make a colour matrix of unique colors (one per station for 11
%     % stations)
%     colorchart =[   0.5 1 0.749; 0.172 0.747 1; 0.75 0.75 0; 1 0.997 0.369; ...
%                     0.5 0.6 1; 0.83 0.6 1; 0.58 0.458 0.541; 0.31 0.4 0.58; ...
%                     0.7 0.78 1; 0.113 0.558 0.98; 0.84 0.91 0.85];
%
%     for i =1: nep  % loop over event pairs
%         hf.diag5(i) = figure;
%         for j = 1: nst % loop over stations
%             subplot(3,1,1)
%             time_shallow = astruc.shallow.(eventpairs_shallow{i}).(stations{j})(startval.shallow:end,1);
%             data2_shallow(j,:) = astruc.shallow.(eventpairs_shallow{i}).(stations{j})(startval.shallow:end,2)';
%             plot(time_shallow, data2_shallow(j,:),'color',[0.5 0.5 0.5])%colorchart(j,:))
%             hold on
%
%             subplot(3,1,2)
%             time_mid = astruc.mid.(eventpairs_mid{i}).(stations{j})(startval.mid:end,1)';
%             data2_mid(j,:) = astruc.mid.(eventpairs_mid{i}).(stations{j})(startval.mid:end,2)';
%             plot(time_mid, data2_mid(j,:),'color',[0.5 0.5 0.5])%colorchart(j,:))
%             hold on
%
%             subplot(3,1,3)
%             if ~isempty(astruc.deep)
%                 time_deep = astruc.deep.(eventpairs_deep{i}).(stations{j})(startval.deep:end,1)';
%                 data2_deep(j,:) = astruc.deep.(eventpairs_deep{i}).(stations{j})(startval.deep:end,2)';
%                 plot(time_deep, data2_deep(j,:),'color',[0.5 0.5 0.5])%colorchart(j,:))
%                 hold on
%             end
%         end
%
%         for lvar1  = 1:size(data2_shallow,2);
%             mean_asft_shallow(lvar1) = mean(data2_shallow(~isnan(data2_shallow(:,lvar1)),lvar1));
%             median_asft_shallow(lvar1) = median(data2_shallow(~isnan(data2_shallow(:,lvar1)),lvar1));
%         end
%
%         subplot(3,1,1)
%         h101 = plot(time_shallow, mean_asft_shallow,'g','linewidth',2);
%         if exist('existing_stats_asft') & ~isempty('existing_stats_asft')
%                 plot(existing_stats_asft.(['pair',num2str(i)]).time_shallow, ...
%                     existing_stats_asft.(['pair',num2str(i)]).mean_asft_shallow,...
%                     'r', 'linewidth',1)
%         end
%         %h102 = plot(time_shallow, median_asft_shallow,'g','linewidth',2);
%         try, title([plot_struct.title,': ',eventpairs_shallow{i},...
%                 ': shallow: knownsep = ', num2str(round(known_sep(i))) ]), catch, end
%         ylims_shallow = get(gca,'ylim');
%         xlims_shallow = get(gca,'xlim');
%
%         for lvar1  = 1:size(data2_mid,2);
%             mean_asft_mid(lvar1) = mean(data2_mid(~isnan(data2_mid(:,lvar1)),lvar1));
%             median_asft_mid(lvar1) = median(data2_mid(~isnan(data2_mid(:,lvar1)),lvar1));
%         end
%         subplot(3,1,2)
%         plot(time_mid, mean_asft_mid,'g','linewidth',2);
%         if exist('existing_stats_asft') & ~isempty('existing_stats_asft')
%                 plot(existing_stats_asft.(['pair',num2str(i)]).time_mid, ...
%                     existing_stats_asft.(['pair',num2str(i)]).mean_asft_mid,...
%                     'r', 'linewidth',1)
%         end
%         %plot(time_mid, median_asft_mid,'g','linewidth',2);
%         try, title([plot_struct.title,': ',eventpairs_mid{i},...
%                 ': mid: knownsep = ', num2str(round(known_sep(i))) ]), catch, end
%         ylims_mid = get(gca,'ylim');
%         xlims_mid = get(gca,'xlim');
%
%         if ~isempty(astruc.deep)
%             for lvar1  = 1:size(data2_deep,2);
%                 mean_asft_deep(lvar1) = mean(data2_deep(~isnan(data2_deep(:,lvar1)),lvar1));
%                 median_asft_deep(lvar1) = median(data2_deep(~isnan(data2_deep(:,lvar1)),lvar1));
%             end
%             subplot(3,1,3)
%             plot(time_deep, mean_asft_deep,'g','linewidth',2);
%             if exist('existing_stats_asft')& ~isempty('existing_stats_asft')
%                 plot(existing_stats_asft.(['pair',num2str(i)]).time_deep, ...
%                     existing_stats_asft.(['pair',num2str(i)]).mean_asft_deep,...
%                     'r', 'linewidth',1)
%         end
%             %plot(time_deep, median_asft_deep,'g','linewidth',2);
%             try, title([plot_struct.title,': ',eventpairs_deep{i}, ...
%                 ': deep: knownsep = ', num2str(round(known_sep(i))) ]), catch, end
%             ylims_deep = get(gca,'ylim');
%             xlims_deep = get(gca,'xlim');
%         else
%             ylims_deep = [NaN NaN];
%             xlims_deep = [NaN NaN];
%         end
%
%         % Now make limits the same and add know separation line
%         actual_ylims = [min([ylims_shallow,ylims_shallow,ylims_deep]),...
%                         max([ylims_shallow,ylims_shallow,ylims_deep, ...
%                         known_sep(i) + 50])];
%         actual_xlims = [min([xlims_shallow,xlims_shallow,xlims_deep]),...
%                         max([xlims_shallow,xlims_shallow,xlims_deep])];
%         subplot(3,1,1)
%         set(gca,'ylim', actual_ylims, 'xlim', actual_xlims, ...
%             'ytick',linspace(actual_ylims(1),actual_ylims(2) ,3), ...
%             'yticklabel',round(linspace(actual_ylims(1),actual_ylims(2) ,3)));
%         h103 = plot(actual_xlims, [known_sep(i), known_sep(i)], 'b--');
%         subplot(3,1,2)
%         set(gca,'ylim', actual_ylims, 'xlim', actual_xlims, ...
%             'ytick',linspace(actual_ylims(1),actual_ylims(2) ,3), ...
%             'yticklabel',round(linspace(actual_ylims(1),actual_ylims(2) ,3)));
%         plot(actual_xlims, [known_sep(i), known_sep(i)], 'b--')
%         subplot(3,1,3)
%         set(gca,'ylim', actual_ylims, 'xlim', actual_xlims, ...
%             'ytick',linspace(actual_ylims(1),actual_ylims(2) ,3), ...
%             'yticklabel',round(linspace(actual_ylims(1),actual_ylims(2) ,3)));
%         plot(actual_xlims, [known_sep(i), known_sep(i)], 'b--')
%         %legend([h101,h102,h103],{'mean','median', 'known separation'},'Location','SouthOutside')
%         clear data2_shallow data2_mid data2_deep
%
%         stats_asft.known_sep = known_sep;
%         stats_asft.(['pair',num2str(i)]).mean_asft_mid = mean_asft_mid;
%         stats_asft.(['pair',num2str(i)]).median_asft_mid = median_asft_mid;
%         stats_asft.(['pair',num2str(i)]).time_mid = time_mid;
%         stats_asft.(['pair',num2str(i)]).mean_asft_shallow = mean_asft_shallow;
%         stats_asft.(['pair',num2str(i)]).median_asft_shallow = median_asft_shallow;
%         stats_asft.(['pair',num2str(i)]).time_shallow = time_shallow;
%         stats_asft.(['pair',num2str(i)]).mean_asft_deep = mean_asft_deep;
%         stats_asft.(['pair',num2str(i)]).median_asft_deep = median_asft_deep;
%         stats_asft.(['pair',num2str(i)]).time_deep = time_deep;
%
%     end
%
%     save(['stats_asft_col',num2str(col),'_normtype', num2str(normtype),'.mat'],'stats_asft')
% end


% if diag(1) ==1
%     disp('------------------------------------')
%     disp('Results for shallow depth')
%     disp('known_sep:')
%     disp(round(known_sep))
%     disp('mean_shallow:')
%     disp(round(mean_shallow))
%     disp('std_shallow:')
%     disp(round(std_shallow))
%     disp('  ')
%     disp('median_shallow:')
%     disp(round(median_shallow))
%     disp('1.483*MAD_shallow:')
%     disp(round(1.483*MAD_shallow))
%     disp('  ')
%     disp('------------------------------------')
%     disp('Results for middle depth')
%     disp('known_sep:')
%     disp(round(known_sep))
%     disp('mean_mid:')
%     disp(round(mean_mid))
%     disp('std_mid:')
%     disp(round(std_mid))
%     disp('  ')
%     disp('median_mid:')
%     disp(round(median_mid))
%     disp('1.483*MAD_mid:')
%     disp(round(1.483*MAD_mid))
%     disp('  ')
%     disp('------------------------------------')
%     if ~isempty(astruc.deep)
%         disp('Results for deep depth')
%         disp('known_sep:')
%         disp(round(known_sep))
%         disp('mean_deep:')
%         disp(round(mean_deep))
%         disp('std_deep:')
%         disp(round(std_deep))
%         disp('  ')
%         disp('median_deep:')
%         disp(round(median_deep))
%         disp('1.483*MAD_deep:')
%         disp(round(1.483*MAD_deep))
%         disp('  ')
%     end
% end

