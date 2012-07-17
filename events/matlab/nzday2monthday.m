function [month day] = nzday2monthday(nzday, year);
% converts nzday (day from start of the year) to month and day

leapyears = 1800:4:2016;

if (year>max(leapyears) | year <min(leapyears))
    error('year input must be between 1800 and 2016')
end
    
dayspermonth = [    31, ... %Jan  
                    28, ... %Feb
                    31, ... %Mar
                    30, ... %Apr
                    31, ... %May
                    30, ... %Jun
                    31, ... %Jul
                    31, ... %Aug
                    30, ... % Sept
                    31, ... % Oct
                    30, ... % Nov
                    31];    %Dec
ind =find(year==leapyears);
if ~isempty(ind)
    dayspermonth(2) = 29; 
end

incsum = cumsum(dayspermonth);

if nzday<=31
    month = 1;
    day = nzday;
elseif nzday>incsum(11)
    month = 12; 
    day = nzday-incsum(11);
else
    indhigher = find(incsum>=nzday);
    month = indhigher(1);
    day = nzday-incsum(indhigher(1)-1);
end
    


