function [data, domfreq_mean4rpair] = combine_smd(medinp, atype,startval,col, normtype,statnum,vel)
%
% COMBINE_SMD combines the shallow, middle and deep estimates of sep or max
% xcorr into the one matrix. 
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
%               If empty the ylimits are chose automatically
% startval      [scalar] startval used to discard intial separation estimates
%               for statistics and plots
% col       [OPTIONAL - scalar]
%               []  => col = 2 - use Taylor series approximation
%               2 - Use Taylor Series approximation
%               3 - Use autucorrelation of unperturbed wave
%               4 - Use autucorrelation of perturbed wave
% normtype  [Optional - Scalar]
%               [] => normtype = 0 - no normalisation
%               0 - No normalisation of separation_true
%               1 - Normalise separation_true (i.e. multiply by domfreq/velocity)% 
% statnum   [Optional]
%               []      => use all stations
%               'all'   => use all stations (default)
%               's1'    => use station 1
%               's12'   => use stations 1 and 2
%
% OUTPUTS:
% 

%
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
    error('ERROR: Invalid value of col in COMBINE_SMD')
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
    error('ERROR: Invalid value of normtype in COMBINE_SMD')
end


% Important hard wired parameters
MADscale = 1.483;

% Initialisations
% initialisation of vector for figure handles
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
eventpairs_mid = fieldnames(astruc.mid);
eventpairs_shallow = fieldnames(astruc.shallow);
if ~isempty(astruc.deep), eventpairs_deep = fieldnames(astruc.deep); end
nep = length(eventpairs_mid);

% Now we can repeat the same for the station names and number of stations.
% Once again we assume that the stations are the same for all event-pairs
% and depth scenarios. 
stations = fieldnames(astruc.mid.(eventpairs_mid{1}));
nst = length(stations);

% here we put all of the data into a big matrix so that we can easily
% compute mean and median values
for i =1: nep  % loop over event pairs
    count =1;
    switch statnum
        case {[],'all'}
            looper = [1:nst];
            usedeep = 1;
        case {'s1'}
            looper = [1:2];
            usedeep = 1;
        case {'s12'}
            looper = [1:5];
            usedeep = 1;
    end
                
    looper
    for j = looper % loop over stations
        
        % Here we look at how many estimates there are for each
        % station for each of the shallow, mid and deep scnearios. This is
        % required for the counters to place data into a single matrix.
        % Note that unlike ANALYSE_DATA we place the shallow, mid and deep
        % values into the one matrix
        swl_shallow = length(astruc.shallow.(eventpairs_shallow{i}).(stations{j})(startval:end,1));
        swl_mid = length(astruc.mid.(eventpairs_mid{i}).(stations{j})(startval:end,1));
        % Now we can put the shallow data in the matrix
        data(count:count+swl_shallow-1,i) = astruc.shallow.(eventpairs_shallow{i}).(stations{j})(startval:end,col); 
        if nargout >1
            domfreq_shallow(j) = mean( ...
                [medinp.domfreq.shallow.(['ref_',eventpairs_shallow{1}(1:14)]).(stations{j}), ...
                medinp.domfreq.shallow.(eventpairs_shallow{i}(16:end)).(stations{j})]);
        end
        if normtype ==1
            data(count:count+swl_shallow-1,i) = data(count:count+swl_shallow-1,i)*domfreq_shallow(j)/vel;
        end
        %Update the counter based on how many values have been taken from shallow
        count = count +swl_shallow-1;      
        
        % Now we can put the mid data in the matrix
        data(count:count+swl_mid-1,i) = astruc.mid.(eventpairs_mid{i}).(stations{j})(startval:end,col); 
        % Compute the dominant frequencies for each station (note all windows considered separately due to vectors) 
        if nargout >1
        domfreq_mid(j) = mean( ...
            [medinp.domfreq.mid.(['ref_',eventpairs_mid{1}(1:14)]).(stations{j}), ...
             medinp.domfreq.mid.(eventpairs_mid{i}(16:end)).(stations{j})]);
        end
       % Normalise the known separations if required
        if normtype ==1
            data(count:count+swl_mid-1,i) = data(count:count+swl_mid-1,i)*domfreq_mid(j)/vel;
        end
        %Update the counter based on how many values have been taken from mid
        count = count +swl_mid-1;

        % Now we do the same for the deep block but we keep it separate
        % because it does not always exist. 
        if ~isempty(astruc.deep) & usedeep==1
            swl_deep = length(astruc.deep.(eventpairs_deep{i}).(stations{j})(startval:end,1));
            data(count:count+swl_deep-1,i) = astruc.deep.(eventpairs_deep{i}).(stations{j})(startval:end,col);
            if nargout > 1
                domfreq_deep(j) = mean( ...
                    [medinp.domfreq.deep.(['ref_',eventpairs_deep{1}(1:14)]).(stations{j}), ...
                    medinp.domfreq.deep.(eventpairs_deep{i}(16:end)).(stations{j})]);
            end
            if normtype ==1
                data(count:count+swl_deep-1,i) = data(count:count+swl_deep-1,i)*domfreq_deep(j)/vel;
            end
            count = count +swl_deep-1;
        end % end of astruc.deep block
    end % end of loop over stations nst  

     % compute average dominant frequency for the receiver pair (i.e.
    % average over all stations) -note we join shallow, mid and deep since
    % we only use on value for the plotting.
    if nargout > 1
        if ~isempty(astruc.deep) & usedeep ==1
            domfreq_mean4rpair(i) = mean([domfreq_shallow,domfreq_mid,domfreq_deep]); 
        else
            domfreq_mean4rpair(i) = mean([domfreq_shallow,domfreq_mid]);
        end
    end
end % end of loop over event pairs nep


