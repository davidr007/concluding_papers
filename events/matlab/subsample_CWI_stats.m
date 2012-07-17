function [CWI_stat_small,CWI_stat_large] = subsample_CWI_stats(CWI_stat,events_oi)
% This function can be used to subsample the CWI_stats matrix into a
% smaller number of events
%
% Note that CWI_stat_small is also re-ordered based on the following
% principles
%       events_oi(1) = master event with location = (0,0,0) (stored in CWI_stat_small(1,1))
%       events_oi(2) = second event with location (x,0,0) (stored in CWI_stat_small(1,2))
%       events_oi(3:end) = events 3,4,5,..,N with free locations. 
%
% INPUTS
% CWI_stat_orig     [char or double]
%                       char   =>   full path to where the CWI_stats matrix of
%                                   interest is saved (note that it assumes it is in a mat 
%                                   file with parameter called CWI_stat
%                       double =>   matrix CWI_stat with more events then those given
%                                   in events_oi  
% events_oi         [double - vector] list of events of interest to pull out of
%                   CWI_stat 
%
% OUTPUTS
% CWI_stat_small   [double] - this is the same as CWI_stat except that it
%                   only includes event pairs with both events listed in events_oi
% CWI_stat_large   [double] - this is the same as CWI_stat_small except
%                   that there are null values assigned to rows where the
%                   CWI constraints do not exist
%
%
% Note the CWI_stat matrices that go in and out of this function have 8
% columns - the ninth column which typically includes the velocity is added
% later. 
%
% David Robinson
% 13 June 2008


if ischar(CWI_stat)
    load(CWI_stat);  % This will load a double (matrix) with name CWI_stat
end
[nin,min] = size(CWI_stat);


ind1 = [];
ind2 = [];
for i = 1:length(events_oi)
    ind1 = [ind1;find(CWI_stat(:,1) == events_oi(i))];
    ind2 = [ind2;find(CWI_stat(:,2) == events_oi(i))];
end
ind = intersect(ind1,ind2);

CWI_stat_small_tmp = CWI_stat(ind,:);

% Now we must correct for the event order when processing
% i.e. we need event_oi(i) on the left and event_oi(j) in the right 
% columns whenever i<j
% This may not have been the same order during the original processing
CWI_stat_small_tmp2 = [];
for i = 1:length(events_oi)-1
    for j = i+1:length(events_oi)
        ind1 = find(CWI_stat_small_tmp(:,1) == events_oi(i) & CWI_stat_small_tmp(:,2) == events_oi(j)); 
        ind2 = find(CWI_stat_small_tmp(:,2) == events_oi(i) & CWI_stat_small_tmp(:,1) == events_oi(j)); 
        if ~isempty(ind1)
            CWI_stat_small_tmp2 = [CWI_stat_small_tmp2;CWI_stat_small_tmp(ind1,:)]; % no swapping required
        elseif ~isempty(ind2)
            CWI_stat_small_tmp2 = [CWI_stat_small_tmp2;CWI_stat_small_tmp(ind2,2),CWI_stat_small_tmp(ind2,1),CWI_stat_small_tmp(ind2,3:end)];
        end           
        
    end   
end


% Here we do the re-rodering
CWI_stat_small = [];
for i = 1:length(events_oi)-1; 
    for j = i+1:length(events_oi);
        %disp([i j])
        ind = find( CWI_stat_small_tmp2(:,1) ==events_oi(i) & ...
                    CWI_stat_small_tmp2(:,2) ==events_oi(j));
        CWI_stat_small = [CWI_stat_small; CWI_stat_small_tmp2(ind,:)];
    end
end

% To use the fortran version of Lstar_3D we can not have blank pairs
% They must be replaced will null values (= -99999) in column 6
nconstr = (length(events_oi)^2-length(events_oi))/2;
CWI_stat_large = zeros(nconstr,min);
[nsmall,msmall] = size(CWI_stat_small);
count_big = 1;
count_small = 1;
for i = 1:length(events_oi)-1; 
    for j = i+1:length(events_oi);
        if count_small <= nsmall
            if (CWI_stat_small(count_small,1) ==events_oi(i) & CWI_stat_small(count_small,2) ==events_oi(j))
                CWI_stat_large(count_big,:) = CWI_stat_small(count_small,:); 
                count_small = count_small+1;
            else
                CWI_stat_large(count_big,:) = [events_oi(i),events_oi(j),-99999*ones(1,6)];
            end
        else
            CWI_stat_large(count_big,:) = [events_oi(i),events_oi(j),-99999*ones(1,6)];
        end
            count_big = count_big +1;           
    end
end
1+1