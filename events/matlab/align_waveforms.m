function [time1,wave1,time2,wave2,details] = align_waveforms(intime1,inwave1,intime2,inwave2,alignval1,alignval2,qa)
% To align two waveforms
% Note that the waves are truncated to the same length as well.
%
%INPUTS:
% intime1   [double - vector] time values for wave 1
% inwave1   [double - vector] wave 1
% intime2   [double - vector] time values for wave 2
% inwave2   [double - vector] wave 2
% aignval1  [double - scalar] time value in time1 for which to align the vector
%           around i.e. tp, ts or some other critical point.
% aignval2  [double - scalar] time value in time2 for which to align the vector
%           around i.e. tp, ts or some other critical point.
% qa        [OPTIONAL]
%               'v'  => verbose comments
%
% OUTPUTS:
% time1   [double - vector] aligned time values for wave 1 (starts at 0 and truncated)
% wave1   [double - vector] aligned wave 1 (truncated as needed)
% time2   [double - vector] aligned time values for wave 2 (starts at 0 and truncated)
% wave2   [double - vector] aligned wave 2 (truncated as needed) 
% details [structure]  contains information about the alignment

% make sure the intime1 and intime2 are row vectors (find_closest bombs if
% they are column vectors)

%% First we must make sure the waves cover roughly the same time bounds. If
% they are vastly different then the following code may not work. 
mintime = max([min(intime1),min(intime2)]);
maxtime = min([max(intime1), max(intime2)]);
% intime1 = intime1(intime1>=mintime & intime1<=maxtime);
% inwave1 = inwave1(intime1>=mintime & intime1<=maxtime);
% intime2 = intime2(intime2>=mintime & intime2<=maxtime);
% inwave2 = inwave2(intime2>=mintime & intime2<=maxtime);
ind1 = find(intime1>=mintime & intime1<=maxtime);
intime1 = intime1(ind1);
inwave1 = inwave1(ind1);
ind2 = find(intime2>=mintime & intime2<=maxtime);
intime2 = intime2(ind2);
inwave2 = inwave2(ind2);

intime1 = intime1(:);
intime2 = intime2(:);
[details.ind1,details.closest1,details.mind1] = find_closest(intime1,alignval1,'euclidean');
[details.ind2,details.closest2,details.mind2] = find_closest(intime2,alignval2,'euclidean');
%disp('ind values')
%disp([details.ind1 details.ind2])

details.trans = details.ind1 -details.ind2;
if details.ind1 == details.ind2
    % Waveforms are already aligned 
    % There is no need to do anything
    time1 = intime1;
    wave1 = inwave1;
    time2 = intime2;    % we must re-set this in case the starting times were different
    wave2 = inwave2;
    details.trans1 = 0;
    details.trans2 = 0;
    details.explanation = 'waveforms were already aligned - no extra alignment needed'
    details.alignval1_out = alignval1 - intime1(1); 
    details.alignval2_out = alignval2 - intime2(1);
    time1 = intime1-intime1(1);
    time2 = intime2-intime2(1);
    % do a quick check 
    if  details.alignval1_out - details.alignval2_out > 0.005
        
        disp(['details.alignval1_out ', num2str(details.alignval1_out)])
        disp(['details.alignval2_out ', num2str(details.alignval2_out)])
        error('you should not be here')
    end
elseif details.ind1 <  details.ind2
    details.explanation = 'wave1 shifted right (end truncated) - wave2 shifted left (strat truncated)';
    details.trans1 = 0;
    details.trans2 = -abs(details.trans);
    details.starttime = intime1(1);
    [time1,wave1,time2,wave2, details.alignval1_out,details.alignval2_out ] ...
        = LOC_doalignment(intime1,inwave1,intime2,inwave2,alignval1,alignval2,details);   
 
    
elseif details.ind1 > details.ind2
        details.explanation = 'wave1 shifted left (start truncated) - wave2 shifted right (end truncated)';
        details.trans1 = -abs(details.trans);
        details.trans2 = 0;
        details.starttime = intime1(1);
        [time2,wave2,time1,wave1,details.alignval2_out, details.alignval1_out] ...
            = LOC_doalignment(intime2,inwave2,intime1,inwave1,alignval2,alignval1,details)   ;
                    
end


function [time2mvright,wave2mvright,time2mvleft,wave2mvleft, ...
                tp4wave2mvright, tp4wave2mvleft] = ...
                                LOC_doalignment(intime2mvright,inwave2mvright,...
                                                intime2mvleft,inwave2mvleft,...
                                                alignval4_wave2mvright,alignval4_wave2mvleft, ...
                                                details)
                                            
% Do the right move first. 
time2mvright = intime2mvright(1:end - abs(details.trans));
wave2mvright = inwave2mvright(1:end - abs(details.trans));
% Do the left move 
time2mvleft = intime2mvleft(abs(details.trans)+1:end);
wave2mvleft = inwave2mvleft(abs(details.trans)+1:end);


%% Re-set intial time value to 0 and check that wavelengths are the same 
% Note the method for this may be clunky - could look at improving
% algorithm later
if time2mvright(1) ~= details.starttime
    deltat = -(time2mvright(1) - details.starttime);
    time2mvright = time2mvright +deltat;
    tp4wave2mvright = alignval4_wave2mvright+deltat;
else
    tp4wave2mvright = alignval4_wave2mvright;
end
if time2mvleft(1) ~= details.starttime
    deltat = -(time2mvleft(1) - details.starttime);
    time2mvleft = time2mvleft +deltat;
    tp4wave2mvleft = alignval4_wave2mvleft+deltat;
else
    tp4wave2mvleft = alignval4_wave2mvleft;
end
% Check that the wavelengths are the same
len_wave2mvright = length(wave2mvright);
len_wave2mvleft = length(wave2mvleft);
if length(time2mvright)~= length(time2mvleft)
    warning(['WARNING: length of time2mvright (', num2str(length(time2mvright)), ...
         ') ~= length of time2mvleft (', num2str(length(time2mvleft)),')'])
    %disp(' The longest wave has been truncated at the end')
    if len_wave2mvright > len_wave2mvleft 
        wave2mvright = wave2mvright(1:length(wave2mvleft));
        mover = time2mvright(1);
        time2mvright = time2mvright(1:length(wave2mvleft))-mover;
        time2mvleft = time2mvright;
    elseif len_wave2mvright < len_wave2mvleft 
        wave2mvleft = wave2mvleft(1:length(wave2mvright));
        mover = time2mvleft(1);
        time2mvleft = time2mvleft(1:length(wave2mvright))-mover;
        time2mvright = time2mvleft;
    end
else % then the waves are already the same length
    mover = 0;
    
end
tp4wave2mvleft = tp4wave2mvleft-mover;
tp4wave2mvright = tp4wave2mvright-mover;

