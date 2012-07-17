function [waves] = grab2Calaveras_waveforms(eventid1,eventid2,station,channel,fbounds, qa_switch)
% Get's two waveforms from the Calaveras data set and returns them based on
% properties defined by wavetype
% INPUTS: 
% eventid1      [str] eventid for firts event
% eventid2      [str] eventid for second event
% station       [str] desired station
% channel       [str] desired channel
% fbounds       [double 1x2] bounds for the bandpass filter [Hcut,Lcut]
% qa_switch     [integer] diagnostic switch
%                   1 => create the plot
%                   0 => do not do any plotting 
%
% OUTPUTS:

%%setup
Hcut = fbounds(1);
Lcut = fbounds(2);


%% get the filenames
wave_store = '/export/storage/davidr/play/HYPODD_examples/Calaveras/';
supportdir = '/export/storage/davidr/play/HYPODD_examples/example2';
wavedir1 = [wave_store,eventid1,filesep];
wavedir2 = [wave_store,eventid2,filesep];
[filename1] = LOC_findthewaveform(wavedir1,station,channel);
[filename2] = LOC_findthewaveform(wavedir2,station,channel);

%% Get the pahse data - needed for the alignments
[A, B, C, D, E, F, G, H, I , J, K, L, M, N, O] = textread([supportdir,'/Calaveras.pha'],'%s%f%f%s%f%f%f%f%f%f%f%f%f%f%f');
ind2phasepicks = strmatch('#',A);

ind2phasepicks = [ind2phasepicks;length(A)];
for g = 1:length(ind2phasepicks)-1
    phasedata.(['event', num2str(O(ind2phasepicks(g)))]).year = B(ind2phasepicks(g));
    phasedata.(['event', num2str(O(ind2phasepicks(g)))]).month = C(ind2phasepicks(g));
    phasedata.(['event', num2str(O(ind2phasepicks(g)))]).day = str2double(D(ind2phasepicks(g)));
    phasedata.(['event', num2str(O(ind2phasepicks(g)))]).hour = E(ind2phasepicks(g));
    phasedata.(['event', num2str(O(ind2phasepicks(g)))]).minute = F(ind2phasepicks(g));
    phasedata.(['event', num2str(O(ind2phasepicks(g)))]).second = G(ind2phasepicks(g));
    for h = ind2phasepicks(g)+1:ind2phasepicks(g+1)-1
        phasedata.(['event', num2str(O(ind2phasepicks(g)))]).(A{h}(3:end)) = B(h);
    end
end

%% Load the waveforms
[event1]= sacread([wavedir1,filename1],'b');
[month1 day1] = nzday2monthday(event1.nzjday,event1.nzyear);
                wavestarttime1 = [event1.nzyear, ...
                    month1, ...
                    day1, ...
                    event1.nzhour, ...
                    event1.nzmin, ...
                    event1.nzsec+event1.nzmsec/1000];
eventorigintime1 = [phasedata.(['event', num2str(eventid1)]).year, ...
                    phasedata.(['event', num2str(eventid1)]).month, ...
                    phasedata.(['event', num2str(eventid1)]).day, ...
                    phasedata.(['event', num2str(eventid1)]).hour, ...
                    phasedata.(['event', num2str(eventid1)]).minute, ...
                    phasedata.(['event', num2str(eventid1)]).second];
timediff1 = etime(wavestarttime1,eventorigintime1);

[event2]= sacread([wavedir2,filename2],'b');
[month2 day2] = nzday2monthday(event2.nzjday, event2.nzyear);
wavestarttime2 = [event2.nzyear, ...
                    month2, ...
                    day2, ...
                    event2.nzhour, ...
                    event2.nzmin, ...
                    event2.nzsec+event2.nzmsec/1000];
eventorigintime2 = [phasedata.(['event', num2str(eventid2)]).year, ...
                    phasedata.(['event', num2str(eventid2)]).month, ...
                    phasedata.(['event', num2str(eventid2)]).day, ...
                    phasedata.(['event', num2str(eventid2)]).hour, ...
                    phasedata.(['event', num2str(eventid2)]).minute, ...
                    phasedata.(['event', num2str(eventid2)]).second];
timediff2 = etime(wavestarttime2,eventorigintime2);
intime1 = timediff1 + event1.b+[0:event1.delta:(length(event1.data)-1)*event1.delta];
inwave1 = event1.data;
intime2 = timediff2 + event2.b+[0:event2.delta:(length(event2.data)-1)*event2.delta];
inwave2 = event2.data;

%=================================================
% ALIGN WAVEFORMS

 alignval1 = phasedata.(['event',eventid1]).(station);%(3:end));
alignval2 = phasedata.(['event',eventid2]).(station);% (3:end));
[time1,alignwave1,time2,alignwave2,align_details] = ...
    align_waveforms(intime1,inwave1,intime2,inwave2,alignval1,alignval2,'v');
time1 = time1-align_details.alignval1_out;
time2 = time2-align_details.alignval2_out;

%%FILTER WAVEFORMS
[filt_wave2, filt_wave1,filt_extras] = filter_waveform(alignwave1,time1, alignwave2,time2,Lcut,Hcut,0);


%% Now let's put it all together
waves.(['event',eventid1]).rawtime = intime1;
waves.(['event',eventid1]).rawdata = inwave1;
waves.(['event',eventid1]).rawtp = alignval1;
waves.(['event',eventid1]).aligntime = time1;
waves.(['event',eventid1]).aligndata = alignwave1;
waves.(['event',eventid1]).aligntp = 0;
waves.(['event',eventid1]).filttime = time1;
waves.(['event',eventid1]).filtdata = filt_wave1;
waves.(['event',eventid1]).filttp = 0;

waves.(['event',eventid2]).rawtime = intime2;
waves.(['event',eventid2]).rawdata = inwave2;
waves.(['event',eventid2]).rawtp = alignval2;
waves.(['event',eventid2]).aligntime = time2;
waves.(['event',eventid2]).aligndata = alignwave2;
waves.(['event',eventid2]).aligntp = 0;
waves.(['event',eventid2]).filttime = time2;
waves.(['event',eventid2]).filtdata = filt_wave2;
waves.(['event',eventid2]).filttp = 0;


%% Diagnostic plotting
if qa_switch ==1
    figure
    subplot(3,1,1)
    l1 = plot(waves.(['event',eventid1]).rawtime,waves.(['event',eventid1]).rawdata,'b');
    hold on
    l2 = plot(waves.(['event',eventid2]).rawtime,waves.(['event',eventid2]).rawdata,'r');   
    legend([l1,l2], {['event',eventid1],['event',eventid2]})
    title('raw waveforms')
    
    subplot(3,1,2)
    l1 = plot(waves.(['event',eventid1]).aligntime,waves.(['event',eventid1]).aligndata,'b');
    hold on
    l2 = plot(waves.(['event',eventid2]).aligntime,waves.(['event',eventid2]).aligndata,'r');
    title('aligned waveforms')
    
    subplot(3,1,3)
    l1 = plot(waves.(['event',eventid1]).filttime,waves.(['event',eventid1]).filtdata,'b');
    hold on
    l2 = plot(waves.(['event',eventid2]).filttime,waves.(['event',eventid2]).filtdata,'r');
    title('filtered waveforms')
    xlabel('time (s)')
end

function [filename] = LOC_findthewaveform(wavedir,station,channel)
    tmp = dir(wavedir);
    for i = 1:length(tmp)
        if strfind(tmp(i).name,station) & strfind(tmp(i).name,channel)
            filename = tmp(i).name;
        end
    end
            
        


