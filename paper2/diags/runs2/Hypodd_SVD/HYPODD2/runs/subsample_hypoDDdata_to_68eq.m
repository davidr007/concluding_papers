function [dtct_data] = subsample_hypoDDdata_to_68eq()
% This function locates the 69 earthquakes of interest then loads all of the 
% hypoDD input data from */orig_69eqonly and sub-samples it to events of
% interest. 

inputdatadir = './orig_example2/';
outputdir = './orig_68eqonly/';

%===================================================
% Find the events of interest

events_oi = load('../../../events_oi.txt');
nE = length(events_oi);
events_oi'
%===================================================
% Now let us load in the event data 
% allevents = load([inputdatadir,'event.dat']);
% [c, ind_allevents, iindevents_oi] = intersect(allevents(:,end), events_oi);
% events69eq = allevents(ind_allevents,:);
% %save([outputdir,'event.dat'],'events69eq','-ascii');
% [n m] = size(events69eq);
% fid = fopen([outputdir,'event.dat'],'w');
% for i = 1:n
%     fprintf(fid, '%d%d%d%d%d%d%d%d%d%d\n', events69eq(i,:)); % write the evnt pair line to file
% end
% fclose(fid)

fid = fopen([inputdatadir,'event.dat']);
event_data = textscan(fid, '%s', 'delimiter', '\n');
fclose(fid)
fid = fopen([outputdir,'event.dat'],'w');
for i = 1:length(event_data{1})
    tmp = event_data{1}(i);
    tmp2 = str2num(tmp{1}); 
    ind1 = find(tmp2(end) == events_oi)
    if ~isempty(ind1)
        %fprintf(fid, '%8.0f %8.0f %7.4f %9.4f %5.3f %3.1f %4.2f %4.2f %4.2f %f \n', tmp2); % write the evnt pair line to file
        fprintf(fid, '%s \n', tmp{1})
    end
end
fclose(fid)


%===================================================
% Now let us copy over the station data
% i.e. no need to sub-sample it
copyfile([inputdatadir,'station.dat'],[outputdir,'station.dat'])

%===================================================
% Now let us have a look at the Calaveras.pha file used by ph2dt
% First read in all the data
create_phasefile('Calaveras.pha',inputdatadir,outputdir,events_oi)

%===================================================
% Now let us have a look at th dt.cc file 
% First read in all the data
create_ttimefile('dt.cc',inputdatadir,outputdir,events_oi)


% =================================

function create_ttimefile(fname,inputdatadir,outputdir,events_oi)
fid = fopen([inputdatadir,fname]);
dt_data = textscan(fid, '%s', 'delimiter', '\n');
fclose(fid);

%Now lets go line by line and write out the data we want
fid = fopen([outputdir,fname],'w');
writerflag = 0;
for i = 1:length(dt_data{1})
    tmp = dt_data{1}(i);
    if strcmp(tmp{1}(1),'#') % we have hit a new event pair
        eventpair = str2num(tmp{1}(2:end));
        ind1 = find(eventpair(1) == events_oi);
        ind2 = find(eventpair(2) == events_oi);
        if ~isempty(ind1) & ~isempty(ind2) % we have a desire dpair
            writerflag = 1;  % set the writer flag to 1 to ensure the following station data is written
            fprintf(fid, '%s\n', tmp{1}); % write the evnt pair line to file
        else 
            writerflag = 0; % set th flag to 0 to ensure the following station data is ignored
        end
    else  % we are looking at station data
        if writerflag ==1
            fprintf(fid, '%s\n', tmp{1}); % write the evnt pair line to file
        end       
    end
end    
fclose(fid)

function create_phasefile(fname,inputdatadir,outputdir,events_oi)
fid = fopen([inputdatadir,fname]);
phase_data = textscan(fid, '%s', 'delimiter', '\n');
fclose(fid);

%Now lets go line by line and write out the data we want
fid = fopen([outputdir,fname],'w');
writerflag = 0;
for i = 1:length(phase_data{1})
    tmp = phase_data{1}(i);
    if strcmp(tmp{1}(1),'#') % we have hit a new event pair
        tmp2 = str2num(tmp{1}(2:end));
        IDinLine = tmp2(end);
        ind = find(IDinLine == events_oi);
        if isempty(ind)
            writerflag = 0;
        elseif ~isempty(ind)
            writerflag = 1;  % set the writer flag to 1 to ensure the following station data is written
            fprintf(fid, '%s\n', tmp{1}); % write the evnt line to file
        end
    else
        if writerflag ==1
            fprintf(fid, '%s\n', tmp{1}); % write the evnt pair line to file
        end       
    end
end    
fclose(fid)
   
        
        
        
   

