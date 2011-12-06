function [] = subsample_hypoDDdata_to_stations(outputdir, number_of_stats)
%
% This function is used to subsample the hypoDD data to desired stations
% (i.e. number_of_stats stations are sub-sampled from 
% ../../CWI/importantstats.txt). A new set of subsampled hypoDD input files
% are saved to outputdir. This function assumes that the input files have 
% already been sub-sampled to the deisred events (for details see
% subsample_hypoDDdata_to_69eq.m).
%

inputdatadir = './orig_68eqonly/';

%=============================================================
% Read in the station data
s = importdata('../../CWI/importantstats.txt');
importantstats_tmp = s.textdata(6:end);
importantstats_data_tmp = s.data;
importantstats = importantstats_tmp(1:number_of_stats);
importantstats_data = importantstats_data_tmp(1:number_of_stats,:);


% ============================================================
% Sub-sample the station file
fid = fopen([inputdatadir,'station.dat']);
station_data = textscan(fid, '%s', 'delimiter', '\n');
fclose(fid);
fid = fopen([outputdir,filesep,'station.dat'],'w');
for i = 1:length(station_data{1})
    tmp = station_data{1}(i);
    statname1 = tmp{1}(1:6);
    for j = 1:length(importantstats)
        statname2 = ['NC',importantstats{j},' '];
        if strcmp(statname1,statname2)
            fprintf(fid, '%s \n', tmp{1});
        end
    end
 %   if ~isempty(ind1)
        %fprintf(fid, '%8.0f %8.0f %7.4f %9.4f %5.3f %3.1f %4.2f %4.2f %4.2f %f \n', tmp2); % write the evnt pair line to file
  %      fprintf(fid, '%s \n', tmp{1})
  %  end
end
fclose(fid);

% ============================================================
% Sub-sample the dt.cc file
events_rep1 = create_ttimefile('dt.cc',inputdatadir,outputdir,importantstats);

%===================================================
% Now let us have a look at the Calaveras.pha file used by ph2dt
% First read in all the data
events_rep2 = create_phasefile('Calaveras.pha',inputdatadir,outputdir,importantstats)

% ============================================================
% Now lets analyse the events we actually have constraints for
% and create the event file
events_oi = union(events_rep1,events_rep2);
whos events_rep1 events_rep2 events_oi

fid = fopen([inputdatadir,'event.dat']);
event_data = textscan(fid, '%s', 'delimiter', '\n');
fclose(fid);
fid = fopen([outputdir,filesep,'event.dat'],'w');
for i = 1:length(event_data{1})
    tmp = event_data{1}(i);
    tmp2 = str2num(tmp{1}); 
    ind1 = find(tmp2(end) == events_oi);
    if ~isempty(ind1)
        %fprintf(fid, '%8.0f %8.0f %7.4f %9.4f %5.3f %3.1f %4.2f %4.2f %4.2f %f \n', tmp2); % write the evnt pair line to file
        fprintf(fid, '%s \n', tmp{1});
    end
end
fclose(fid);



function [events_rep] = create_ttimefile(fname,inputdatadir,outputdir,importantstats)
fid = fopen([inputdatadir,fname]);
dt_data = textscan(fid, '%s', 'delimiter', '\n');
fclose(fid);

%Now lets go line by line and write out the data we want
fid = fopen([outputdir,filesep,fname],'w');
writerflag = 1;
tmp_cell = {};
events_rep = [];
for i = 1:length(dt_data{1})
    tmp = dt_data{1}(i);
    if strcmp(tmp{1}(1),'#') % we have hit a new event pair
        %First check to see if we have already got some data in the
        %temporary cell and if so write it to file
        if length(tmp_cell) >1  % i.e. if we have a station hit
            tmp_events = str2num(tmp_cell{1}(2:end));
            events_rep = [events_rep, tmp_events(1), tmp_events(2)];
            for k=1:length(tmp_cell)
                fprintf(fid, '%s \n', tmp_cell{k});
            end
        end
        clear tmp_cell  % clear the slate to start a new temporary cell array
        tmp_cell{1} = tmp{1}; % store line into a temporary cell
    else
        % Now lets look for the stations
        statname1 = tmp{1}(1:6);
        for j = 1:length(importantstats)
            statname2 = ['NC',importantstats{j},' '];
            if strcmp(statname1,statname2)
                tmp_cell = [tmp_cell, tmp{1}];
                %fprintf(fid, '%s \n', tmp{1});
            end  % end if on station name compare
        end % end for loop over station names
    end % end if block that sorts between event pairs and station data
end  % end of main loop over all data

% Now we must check the last tmp_cell array
if length(tmp_cell) >1  % i.e. if we have a station hit
    for k=1:length(tmp_cell)
        fprintf(fid, '%s \n', tmp_cell{k});
    end
end
fclose(fid);

events_rep = unique(events_rep);
%===============================================
function events_rep = create_phasefile(fname,inputdatadir,outputdir,importantstats)
fid = fopen([inputdatadir,fname]);
phase_data = textscan(fid, '%s', 'delimiter', '\n');
fclose(fid);

% Now lets go line by line and write out the data we want
fid = fopen([outputdir,filesep,fname],'w');
writerflag = 0;
events_rep = [];

for i = 1:length(phase_data{1})
    tmp = phase_data{1}(i);
    if strcmp(tmp{1}(1),'#') % we have hit a new event 
        fprintf(fid, '%s\n', tmp{1}); % write the evnt line to file        
        tmp2 = str2num(tmp{1}(2:end));
        IDinLine = tmp2(end);
        events_rep = [events_rep, IDinLine];
    else % we have station data 
        % Now lets look for the stations
        statname1 = tmp{1}(1:6);
        for j = 1:length(importantstats)
            statname2 = ['NC',importantstats{j},' '];
            if strcmp(statname1,statname2)
                fprintf(fid, '%s\n', tmp{1}); % write the evnt line to file        
            end  % end if on station name compar
        end
    end
end       
  
fclose(fid)

events_rep = unique(events_rep);
   
        
        
        
   




