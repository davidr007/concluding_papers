function [meanR, meanS, omega2] = plot_flinders(fname_u, fname_p, sample_rate)

% fname_u   [string] full path file name to the wave from the the unperturbed
%           source
% fname_p   [string] full path file name to the wave from the perturbed
%           source. Must be recorded on the same instrument as fname_u
% sample_rate   [float] sample rate of instrument. 

if nargin==0
    %fname_u = '/export/storage/davidr/play/flinders/2004/117/FR01.HHZ.2004:117:22:22:19';
    %fname_p = '/export/storage/davidr/play/flinders/2004/127/FR01.HHZ.2004:127:01:12:45';  %plot_par.trans2 = 0.6
    fname_u = '/export/storage/davidr/play/flinders/2004/220/FR01.HHZ.2004:220:11:31:19';
    fname_p = '/export/storage/davidr/play/flinders/2004/131/FR01.HHZ.2004:131:01:44:42';  %plot_par.trans2 = 0.6
    fname_u = 'FR01.HHZ.2004_220_11_31_19';
    fname_p = 'FR01.HHZ.2004_131_01_44_42';  %plot_par.trans2 = 0.6
    sample_rate = 200
end


wave_u_full = load(fname_u);
wave_p_full = load(fname_p);


src_type = 2 %(=2 => double couple source in a 3D elastic medium)

%Now setup the plotting structure for PLOT_SEIS_XCORR
plot_par.trans1 = 0.05;    % horizontal translation of event1(:,1) [time]
plot_par.trans2 =  0%0.585;% 0.006;         % horizontal translation of event2(:,1) [time]
plot_par.start_arrivwin = 29.7;
plot_par.end_arrivwin = 32.7;
plot_par.start_codawin = 40;
plot_par.end_codawin = 43;
plot_par.event1name = 'event1';
plot_par.event2name = 'event2';
plot_par.parrival = 29.5;
plot_par.twindow = 10;

timevector_u = 0:1/sample_rate: (length(wave_u_full)-1)*(1/sample_rate);
timevector_p = 0:1/sample_rate: (length(wave_p_full)-1)*(1/sample_rate);

% % Now we need to mean shift both waveforms so that they are around 0
% wave_u_full = wave_u_full - median(wave_u_full);
% wave_p_full = wave_p_full - median(wave_p_full);
% 
% % Now we need High Pass filter the signals to remove the very low freq
% wave_u_full = apply_butter(wave_u_full, 5,'high', 'time',1, 200);
% wave_p_full = apply_butter(wave_p_full, 5,'high', 'time',1, 200);
% 
% % Now we need to mean shift both waveforms so that they are around 0
% wave_u_full = wave_u_full - median(wave_u_full);
% wave_p_full = wave_p_full - median(wave_p_full);
% 
% %Now we low pass filter the signals to remove all of the very high freq
% wave_u_full = apply_butter(wave_u_full, 5,'low', 'time',3, 200);
% wave_p_full = apply_butter(wave_p_full, 5,'low', 'time',3, 200);
% 
% % Now we need to mean shift both waveforms so that they are around 0
% wave_u_full = wave_u_full - median(wave_u_full);
% wave_p_full = wave_p_full - median(wave_p_full);
% 
% % 
% % Finally, we scale both waveforms
% wave_u_full = wave_u_full/max(abs(wave_u_full));
% wave_p_full = wave_p_full/max(abs(wave_p_full));


Lcut=[3 4 5 6 7 8 9 10];
Hcut=[1 2 3 4 5 6 7 8];

Lcut = 5;
Hcut = 1;

for i=1:length(Lcut);
    % Now we need to mean shift both waveforms so that they are around 0
    wave_u_full = wave_u_full - median(wave_u_full);
    wave_p_full = wave_p_full - median(wave_p_full);
    generate_plots(wave_u_full, timevector_u, wave_p_full, timevector_p, ...
        'original waveforms (step 1)', [0,2],[-2500/200,2500/200] )

    %Now we should taper them to zero at the end
    u_taper = ztaper_tsignal(length(wave_u_full), 4001, 1);
    p_taper = ztaper_tsignal(length(wave_p_full), 4001, 1);
    wave_u_full = wave_u_full.*u_taper;
    wave_p_full = wave_p_full.*p_taper;
    generate_plots(wave_u_full, timevector_u, wave_p_full, timevector_p, ...
        'Tapered to zero at each end (step 2)', [0,2],[-2500/200,2500/200] )


    % Now we need High Pass filter the signals to remove the very low freq
    wave_u_full = apply_butter(wave_u_full, 5,'high', 'time',Hcut(i), 200);
    wave_p_full = apply_butter(wave_p_full, 5,'high', 'time',Hcut(i), 200);
    generate_plots(wave_u_full, timevector_u, wave_p_full, timevector_p, ...
        'High Pass 1Hz (step 3)', [0,40],[-0.25, 0.25] )



    % Now we need to mean shift both waveforms so that they are around 0
    wave_u_full = wave_u_full - median(wave_u_full);
    wave_p_full = wave_p_full - median(wave_p_full);
    generate_plots(wave_u_full, timevector_u, wave_p_full, timevector_p, ...
        'DC Removal of mean (step 4)',[0,40],[-0.25, 0.25] )

      %Now we should taper them to zero at the end - we need to do this again
  %because of the second DC shift
    u_taper = ztaper_tsignal(length(wave_u_full), 232, 1);
    p_taper = ztaper_tsignal(length(wave_p_full), 232, 1);
    wave_u_full = wave_u_full.*u_taper;
    wave_p_full = wave_p_full.*p_taper;
    generate_plots(wave_u_full, timevector_u, wave_p_full, timevector_p, ...
        'Tapered to zero at each end (step 2)', [0,2],[-2500/200,2500/200] )
        
        
    %Now we low pass filter the signals to remove all of the very high freq

    wave_u_full = apply_butter(wave_u_full, 5,'low', 'time',Lcut(i), 200);
    wave_p_full = apply_butter(wave_p_full, 5,'low', 'time',Lcut(i),200);
    generate_plots(wave_u_full, timevector_u, wave_p_full, timevector_p, ...
        ['Low Pass ',num2str(Lcut(i)),'Hz (step 5)'],[0,40],[-0.25, 0.25] )



    % Now we need to mean shift both waveforms so that they are around 0
    wave_u_full = wave_u_full - median(wave_u_full);
    wave_p_full = wave_p_full - median(wave_p_full);
    generate_plots(wave_u_full, timevector_u, wave_p_full, timevector_p, ...
        'DC Removal of mean (step 6)', [0,40],[-0.25, 0.25] )

    %
    % Finally, we scale both waveforms
    wave_u_full = wave_u_full/max(abs(wave_u_full));
    wave_p_full = wave_p_full/max(abs(wave_p_full));
    omega2(i) = generate_plots(wave_u_full, timevector_u, wave_p_full, timevector_p, ...
        'Scale Both Waveforms (step 7)', [0,40],[-0.25, 0.25] )


    [norm_maxcorr,sep] = plot_seis_xcorr(plot_par,[timevector_u(:),wave_u_full(:)],[timevector_p(:),wave_p_full(:)],src_type);
    
    start1 = length(sep) - floor(length(sep)/2);
    end1 = length(sep);
    meanR(i) = mean(norm_maxcorr(start1:end1,2));
    meanS(i) = mean(sep(start1:end1,2));
    close all

end