%% IMPORTANT SWITCHES ================================
clear all

waveform_figmode = 0;   % 0=> waveform figures designed for portrait half page
% 1=> waveform figures designed for landscape full page
waveform_cliptime = 1;  % 0=> show entire waveforms in plots
% 1=> clip waveforms to tp-tstart<=t<=tp-tstart+twidth
waveform_rotate = 0;    % 0=> Use original waveforms N,E,Z
% 1=> Use rotated waveforms T,R,Z

tstart.MORW = 5;
twidth.MORW = 75;
tstart.BLDU = 5;
twidth.BLDU = 40;
tstart.KLBR = 5;
twidth.KLBR = 75;

tsurface.MORW = 27;
tsurface.BLDU = 11;
tsurface.KLBR = 30;

tshear.MORW = 19;
tshear.BLDU = 8;
tshear.KLBR = 20;

if waveform_rotate ==0
    channel_list = {'N','E','Z'}
elseif waveform_rotate ==1
    channel_list = {'R','T','Z'}
end


%% MORW ==============================================
Kal_CWI.Event1.MORW.(channel_list{1}).fname = '2005264224434.00.MORW.SHN';
Kal_CWI.Event1.MORW.(channel_list{2}).fname = '2005264224442.85.MORW.SHE';
Kal_CWI.Event1.MORW.(channel_list{3}).fname = '2005264224443.41.MORW.SHZ';


Kal_CWI.Event2.MORW.(channel_list{2}).fname = '2005264225837.21.MORW.SHE';
Kal_CWI.Event2.MORW.(channel_list{1}).fname = '2005264225846.30.MORW.SHN';
Kal_CWI.Event2.MORW.Z.fname = '2005264225847.10.MORW.SHZ';

Kal_CWI.Event3.MORW.(channel_list{1}).fname = '2005265035136.96.MORW.SHN';
Kal_CWI.Event3.MORW.(channel_list{2}).fname = '2005265035147.80.MORW.SHE';
Kal_CWI.Event3.MORW.Z.fname = '2005265035149.30.MORW.SHZ';

Kal_CWI.Event4.MORW.(channel_list{1}).fname = '2005265183329.96.MORW.SHN';
Kal_CWI.Event4.MORW.(channel_list{2}).fname = '2005265183333.10.MORW.SHE';
Kal_CWI.Event4.MORW.Z.fname = '2005265183334.00.MORW.SHZ';

%% BLDU =============================================
Kal_CWI.Event1.BLDU.(channel_list{1}).fname = '2005264224442.09.BLDU.BHN';
Kal_CWI.Event1.BLDU.(channel_list{2}).fname = '2005264224445.95.BLDU.BHE';
Kal_CWI.Event1.BLDU.Z.fname = '2005264224448.75.BLDU.BHZ';

Kal_CWI.Event2.BLDU.(channel_list{2}).fname = '2005264225845.65.BLDU.BHE';
Kal_CWI.Event2.BLDU.(channel_list{1}).fname = '2005264225857.45.BLDU.BHN';
Kal_CWI.Event2.BLDU.Z.fname = '2005264225859.50.BLDU.BHZ';

Kal_CWI.Event3.BLDU.(channel_list{2}).fname = '2005265035143.84.BLDU.BHE';
Kal_CWI.Event3.BLDU.Z.fname = '2005265035149.54.BLDU.BHZ';
Kal_CWI.Event3.BLDU.(channel_list{1}).fname = '2005265035158.54.BLDU.BHN';

Kal_CWI.Event4.BLDU.Z.fname = '2005265183348.84.BLDU.BHZ';
Kal_CWI.Event4.BLDU.(channel_list{2}).fname = '2005265183358.00.BLDU.BHE';
Kal_CWI.Event4.BLDU.(channel_list{1}).fname = '2005265183359.00.BLDU.BHN';

% KLBR ==============================================
Kal_CWI.Event1.KLBR.Z.fname = '2005264224436.84.KLBR.SHZ';
Kal_CWI.Event1.KLBR.(channel_list{1}).fname = '2005264224440.79.KLBR.SHN';
Kal_CWI.Event1.KLBR.(channel_list{2}).fname = '2005264224445.25.KLBR.SHE';

Kal_CWI.Event2.KLBR.Z.fname = '2005264225854.45.KLBR.SHZ';
Kal_CWI.Event2.KLBR.(channel_list{2}).fname = '2005264225857.40.KLBR.SHE';
Kal_CWI.Event2.KLBR.(channel_list{1}).fname = '2005264225858.00.KLBR.SHN';

Kal_CWI.Event3.KLBR.(channel_list{1}).fname = '2005265035135.65.KLBR.SHN';
Kal_CWI.Event3.KLBR.Z.fname = '2005265035137.65.KLBR.SHZ';
Kal_CWI.Event3.KLBR.(channel_list{2}).fname = '2005265035144.15.KLBR.SHE';

Kal_CWI.Event4.KLBR.Z.fname = '2005265183341.59.KLBR.SHZ';
Kal_CWI.Event4.KLBR.(channel_list{1}).fname = '2005265183344.25.KLBR.SHN';
Kal_CWI.Event4.KLBR.(channel_list{2}).fname = '2005265183349.95.KLBR.SHE';

% MEEK ==============================================
Kal_CWI.Event1.MEEK.(channel_list{2}).fname = '2005264224442.53.MEEK.SHE';
Kal_CWI.Event1.MEEK.(channel_list{1}).fname = '2005264224448.69.MEEK.SHN';
Kal_CWI.Event1.MEEK.Z.fname = '2005264224459.08.MEEK.SHZ';

Kal_CWI.Event2.MEEK.Z.fname = '2005264225839.39.MEEK.SHZ';
Kal_CWI.Event2.MEEK.(channel_list{1}).fname = '2005264225849.33.MEEK.SHN';
Kal_CWI.Event2.MEEK.(channel_list{2}).fname = '2005264225858.83.MEEK.SHE';

Kal_CWI.Event3.MEEK.Z.fname = '2005265035131.83.MEEK.SHZ';
Kal_CWI.Event3.MEEK.(channel_list{1}).fname = '2005265035136.28.MEEK.SHN';
Kal_CWI.Event3.MEEK.(channel_list{2}).fname = '2005265035155.78.MEEK.SHE';

Kal_CWI.Event4.MEEK.(channel_list{1}).fname = '2005265183334.14.MEEK.SHN';
Kal_CWI.Event4.MEEK.Z.fname = '2005265183356.78.MEEK.SHZ';
Kal_CWI.Event4.MEEK.(channel_list{2}).fname = '2005265183358.19.MEEK.SHE';

% KMBL ==============================================
Kal_CWI.Event1.KMBL.(channel_list{1}).fname = '2005264224443.20.KMBL.BHN';
Kal_CWI.Event1.KMBL.(channel_list{2}).fname = '2005264224452.75.KMBL.BHE';
Kal_CWI.Event1.KMBL.Z.fname = '2005264224455.85.KMBL.BHZ';

Kal_CWI.Event2.KMBL.Z.fname = '2005264225841.85.KMBL.BHZ';
Kal_CWI.Event2.KMBL.(channel_list{2}).fname = '2005264225845.20.KMBL.BHE';
Kal_CWI.Event2.KMBL.(channel_list{1}).fname = '2005264225846.65.KMBL.BHN';

Kal_CWI.Event3.KMBL.(channel_list{1}).fname = '2005265035143.45.KMBL.BHN';
Kal_CWI.Event3.KMBL.Z.fname = '2005265035151.65.KMBL.BHZ';
Kal_CWI.Event3.KMBL.(channel_list{2}).fname = '2005265035158.20.KMBL.BHE';

Kal_CWI.Event4.KMBL.(channel_list{2}).fname = '2005265183346.10.KMBL.BHE';
Kal_CWI.Event4.KMBL.Z.fname = '2005265183351.50.KMBL.BHZ';
Kal_CWI.Event4.KMBL.(channel_list{1}).fname = '2005265183353.75.KMBL.BHN';


if strcmp(filesep,'/')  % in linux
    if waveform_rotate ==0
        wavedir = '/export/storage/davidr/play/Kalannie/sac/';
    else
        wavedir = '/export/storage/davidr/play/Kalannie/sac/tmp/';
    end
else  % in windows
    if waveform_rotate ==0
        wavedir = 'C:\workstuff\Kalannie\sac\';
    else
        wavedir = 'C:\workstuff\Kalannie\sac\tmp';
    end
end
statnames = {'MORW', 'BLDU', 'KLBR', 'MEEK', 'KMBL'};
statnames = {'MORW', 'BLDU', 'KLBR'};
statlocations = [   -29.068300   116.038800; ... % MORW
    -30.614700  116.709100; ... % BLDU
    -31.591500   117.754600; ... %KLBR
    -26.638000   118.615000; ... %MEEK
    -31.366900   121.882100]; %KMBL

%local_store = './Kalannie_figs/';  % directory in which to store pictures
plotcounter = 1;
Hcut = 1;
Lcut = 2;

tp.event1 = [153.8,144.4,144.95,...% row 1
    130.8,124.2,126.9,... % row 2
    147.75+0.25,151.8+0.1,143.5+0.1,... % row3
    171,160.6,178.5,... % row 4
    181.5,168.8,172.2];
tp.event2 = [65.4-0.25,64.6-0.25,74.5-0.25,... % row1
    39.4-0.4,37.2-0.2,51-0.25,... %row2
    54.3,58-0.2,55.5-0.1-0.45,...%row3
    94.5,104.4,85.65,...%row4
    102.3,106.5,102.85];%row5
tp.event3 = [109-0.1,96.65-0.05,98.1,...%row1
    72.65-0.2,81.7-0.2,87.3-0.2,...%row2
    111.25,109.1+0.3-0.2,102.6+0.7-0.5,...%row3
    142,146-0.1, 123,...%row4
    139.7,130.9,124.4]; %row5
tp.event4 = [89.8,85.75,86.65,...%row1
    46.1-0.15,56.2-0.05,47.05-0.15,...%row2
    76.54+0.15, 79.28, 70.91+0.2, ...%76.5+1,79.2+1,70.8+1,...%row3
    118,94.6,94.4,...%row4
    103.5,105.2,110.6];%row5

% landscape dimensions
hor_size = 6;
vert_size = 2.5;
vert_space = 1;
hor_space = 1;
land_dim.counter7 = [hor_space, vert_space, hor_size, vert_size];
land_dim.counter8 = [2*hor_space+hor_size, vert_space, hor_size, vert_size];
land_dim.counter9 = [3*hor_space+2*hor_size, vert_space, hor_size, vert_size];

land_dim.counter4 = [hor_space, 2*vert_space + vert_size, hor_size, vert_size];
land_dim.counter5 = [2*hor_space+hor_size, 2*vert_space + vert_size, hor_size, vert_size];
land_dim.counter6 = [3*hor_space+2*hor_size, 2*vert_space + vert_size, hor_size, vert_size];

land_dim.counter1 = [hor_space, 3*vert_space + 2*vert_size, hor_size, vert_size];
land_dim.counter2 = [2*hor_space+hor_size, 3*vert_space + 2*vert_size, hor_size, vert_size];
land_dim.counter3 = [3*hor_space+2*hor_size, 3*vert_space + 2*vert_size, hor_size, vert_size];

figure
for i = 1:2  % loop over the events of interest (usually only 2)
    counter = 0;
    for j = 1:1 % length(statnames)
        for k = 1:3
            counter = counter+1;
            if plotcounter ==1
                if waveform_rotate ==0
                    channel = 'N';
                elseif waveform_rotate ==1
                    channel = 'R';
                else
                    error('invalid value for waveform_rotate')
                end
                plotcounter = plotcounter +1;  % convert it to value 2
            elseif plotcounter ==2
                channel = 'Z';
                plotcounter = plotcounter +1;  % convert to value 3
            elseif plotcounter == 3
                if waveform_rotate ==0
                    channel = 'E';
                elseif waveform_rotate ==1
                    channel = 'T';
                else
                    error('invalid value for waveform_rotate')
                end
                plotcounter = 1;
            else
                error('you should never get here')
            end

            fname = Kal_CWI.(['Event',num2str(i)]).(statnames{j}).(channel).fname;
            [tmp101] = sacread([wavedir,fname],'b');
            Kal_CWI.(['Event',num2str(i)]).(statnames{j}).(channel).event = tmp101;
            Kal_CWI.(['Event',num2str(i)]).(statnames{j}).(channel).event.tp = tp.(['event',num2str(i)])(counter);
            timevec = 0:tmp101.delta:tmp101.delta*length(tmp101.data)-tmp101.delta;
            Kal_CWI.(['Event',num2str(i)]).(statnames{j}).(channel).timevec = timevec;


            [filt_wave2, filt_wave1,filt_extras] = filter_waveform(tmp101.data,timevec, tmp101.data,timevec,Lcut,Hcut,0);
            Kal_CWI.(['Event',num2str(i)]).(statnames{j}).(channel).fdom = filt_extras.end.u.fdom;

            ind = find(timevec<=400);


            if i ==1 % only plot the waveforms for event 1
                subplotnum.(['num',num2str(counter)]) = subplot(3,3,counter);

                if waveform_cliptime ==0
                    plot(timevec(ind), filt_wave2(ind)/max(abs(filt_wave2)) )%tmp101.data(ind)/max(abs(tmp101.data)));
                    hold on
                    plot(tp.(['event',num2str(i)])(counter)*[1 1],[-1 1],'r')
                elseif waveform_cliptime ==1
                    plot(timevec(ind)-tp.(['event',num2str(i)])(counter), filt_wave2(ind)/max(abs(filt_wave2)))
                    hold on
                    plot([0 0],[-1 1],'r')
                    plot(tsurface.(statnames{j})*[1 1], [-1 1],'r')
                    plot(tshear.(statnames{j})*[1 1], [-1 1],'r')

                    set(gca,'xlim', [-tstart.(statnames{j}),-tstart.(statnames{j})+twidth.(statnames{j})])
                end

                if counter ==1 |  counter ==2 |  counter ==3

                    title(channel)

                end
                %if counter == 3*length(statnames) | counter == 3*length(statnames)-1 | counter == 3*length(statnames)-2
                xlabel('t (s)')
                %end


            end

        end


    end


    % convert the figure to landscape mode if desired
    if waveform_figmode ==1
        if length(statnames) ~=3
            error('can only convert to landscape for three stations')
        else
            counter = 1;
            for i = 1:9
                axes(eval(['subplotnum.num',num2str(counter)]))
                set(gcf,'units','centimeters')
                set(gcf,'PaperPosition',[0.6, 6.3,24,15])
                set(gca,'units','centimeters')
                posydata = eval(['land_dim.counter',num2str(counter)]);
                set(gca,'position', posydata)
                counter = counter +1;
            end
        end
    end

end % END lopp over the four events

%error('stop here while coding')
ylimits = [-0.4, 0.4; -0.4,0.4; -0.4,0.4; -0.05, 0.05; -0.05,0.05; -0.05,0.05; -0.4, 0.4; -0.4,0.4;-0.4,0.4];

plotcounter = 1;
for i = 1:1 %1:3   % loop over the first event
    for m= 2:2 %i+1:4 % loop over the second event

        counter = 3;
        for j = 1:1 % length(statnames)  % loop over the stations

            for k = 1:3  % loop over the channels
                counter = counter+1;
                if plotcounter ==1

                    if waveform_rotate ==0
                        channel = 'N';
                    elseif waveform_rotate ==1
                        channel = 'R';
                    end
                    plotcounter = plotcounter +1;  % convert it to value 2
                elseif plotcounter ==2
                    channel = 'Z';
                    plotcounter = plotcounter +1;  % convert to value 3
                elseif plotcounter == 3
                    if waveform_rotate ==0
                        channel = 'E';
                    elseif waveform_rotate ==1
                        channel = 'T';
                    end
                    plotcounter = 1;
                else
                    error('you should never get here')
                end


                %=================================================
                % ALIGN WAVEFORMS
                intime1 = Kal_CWI.(['Event',num2str(i)]).(statnames{j}).(channel).timevec;
                inwave1 = Kal_CWI.(['Event',num2str(i)]).(statnames{j}).(channel).event.data;
                intime2 = Kal_CWI.(['Event',num2str(m)]).(statnames{j}).(channel).timevec;
                inwave2 = Kal_CWI.(['Event',num2str(m)]).(statnames{j}).(channel).event.data;
                alignval1 = Kal_CWI.(['Event',num2str(i)]).(statnames{j}).(channel).event.tp;
                alignval2 = Kal_CWI.(['Event',num2str(m)]).(statnames{j}).(channel).event.tp;
                [time1,alignwave1,time2,alignwave2,tester101] = ...
                    align_waveforms(intime1,inwave1,intime2,inwave2,alignval1,alignval2,'v');


                % FILTER WAVEFORMS
                [filt_wave2, filt_wave1,filt_extras] = filter_waveform(alignwave1,time1, alignwave2,time2,Lcut,Hcut,0);

                time1 = time1-tester101.alignval1_out;
                time2 = time2 - tester101.alignval1_out;
                arrival_lims = [-0.2,3];
                ind = find(time1>arrival_lims(1) & time1<arrival_lims(2));


%                 %figure(fh1)
%                 subplot(3,3,counter)
%                 plot(time1(ind),filt_wave1(ind)/max(abs(filt_wave1)),'b')
%                 hold on
%                 plot(time2(ind),filt_wave2(ind)/max(abs(filt_wave2)),'r')
%                 set(gca,'xlim',arrival_lims,'ylim',ylimits(counter,:))
%                 if counter == 1 |  counter == 4  | counter == 7
%                     ylabel(statnames{j})
%                 end
%                 if counter ==1 |  counter ==2 |  counter ==3
%                     title(channel)
%                 end
%                 if counter == 3*length(statnames) | counter == 3*length(statnames)-1 | counter == 3*length(statnames)-2
%                     xlabel('t (s)')
%                 end
                 pairname = ['pair',num2str(i),'_',num2str(m)];


                %figure(fh2)
                subplot(3,3,counter)
                plot(time1(ind),filt_wave1(ind)/max(abs(filt_wave1)),'b')
                hold on
                plot(time2(ind),filt_wave2(ind)/max(abs(filt_wave2)),'r')
                %set(gca,'xlim',arrival_lims,'ylim',ylimits(counter,:))
                set(gca,'xlim', arrival_lims,'ylim', [-0.4, 0.4])

                plot(time1(ind),filt_wave1(ind)/max(abs(filt_wave1)),'b.')
                plot(time2(ind),filt_wave2(ind)/max(abs(filt_wave2)),'r.')
                %if counter == 3*length(statnames) | counter == 3*length(statnames)-1 | counter == 3*length(statnames)-2
                xlabel('t (s)')
                %end



            end
        end

    end
end

%%=======================================
% plot the CWI Data 
%Hcut = 1;
%Lcut = 2;
tw = 5;
vs = 2728.3;

tstart.MORW = 5;
twidth.MORW = 22; % 30
tstart.BLDU = 5;
twidth.BLDU = 7; % 15
tstart.KLBR = 5;
twidth.KLBR = 25;  % 30

datastore = '../../../thesis_version2/diags/eq_location_simple/';

if waveform_rotate==0 & waveform_cliptime ==0
    load([datastore,'Kalannie_CWI_',num2str(Hcut),'to',num2str(Lcut),'Hz_tw',num2str(tw),'.mat'],'Kalpairs','Kal_CWI')
elseif  waveform_rotate==1   &  waveform_cliptime==1
    load([datastore,'Kalannie_CWI_rotated_timeclip_',num2str(Hcut),'to',num2str(Lcut),'Hz_tw',num2str(tw),'.mat'],'Kalpairs','Kal_CWI')
elseif waveform_rotate==0 & waveform_cliptime ==1
    load([datastore,'Kalannie_CWI_timeclip_',num2str(Hcut),'to',num2str(Lcut),'Hz_tw',num2str(tw),'.mat'],'Kalpairs','Kal_CWI')

else
    error('INVALID combination of waveform_rotate and waveform_cliptime')
end    

subplot(3,3,[7:9])
for k = 1:3 % loop over the channels
    tmp_sep = Kalpairs.pair1_2.sep.MORW.(channel_list{k});
    plot(tmp_sep(:,1),tmp_sep(:,3),'color','k')
    hold on
end
set(gca,'ylim',[0 1000],'xlim',[-30,150])
plot(tstart.MORW*[1 1], get(gca,'ylim'), 'color','k','linewidth',2)%[0.8 0.8 0.8])
plot((tstart.MORW+twidth.MORW)*[1 1], get(gca,'ylim'), 'color','k','linewidth',2)%[0.8 0.8 0.8])
xlabel('t (s)')
ylabel('{\it\delta_{CWIN}}  (m)')

print -depsc waves_and_CWIest.eps