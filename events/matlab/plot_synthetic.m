function [meanR, meanS, omega2] = plot_flinders(fname_u, fname_p, diag_switch)

% fname_u   [string] full path file name to the wave from the the unperturbed
%           source
% fname_p   [string] full path file name to the wave from the perturbed
%           source. Must be recorded on the same instrument as fname_u
% sample_rate   [float] sample rate of instrument. 

if nargin==0
    %fname_u = '/export/storage/davidr/play/flinders/2004/117/FR01.HHZ.2004:117:22:22:19';
    %fname_p = '/export/storage/davidr/play/flinders/2004/127/FR01.HHZ.2004:127:01:12:45';  %plot_par.trans2 = 0.6
    fname_u = 'receiver000036.txt';
    fname_p = 'receiver000039.txt';  %plot_par.trans2 = 0.6
    %fname_u = 'run0021/station208.txt';
    %fname_p = 'run0022/station208.txt';  
    %sample_rate = 1250
    diag_switch = 1
end
known_sep = sqrt(2*(6*20)^2);
%known_sep = 0.0228

wave_u_full_matrix = load(fname_u);
wave_p_full_matrix = load(fname_p);

wave_u_full = wave_u_full_matrix(:,11);
wave_p_full = wave_p_full_matrix(:,11);

src_type = 1    %(1 => double couple source in a 2D elastic medium
                %(2 => double couple source in a 3D elastic medium
                %(3 => 3D double couple source [delta source variation])
%Now setup the plotting structure for PLOT_SEIS_XCORR
plot_par.trans1 = 0;    % horizontal translation of event1(:,1) [time]
plot_par.trans2 =  0%0.585;% 0.006;         % horizontal translation of event2(:,1) [time]
 plot_par.start_arrivwin = 1;%30.5;
    plot_par.end_arrivwin = 2;%34.5;
    plot_par.start_codawin = 4;%40;
    plot_par.end_codawin = 5;%44;
    plot_par.parrival = 1.1;
plot_par.event1name = 'event1';
plot_par.event2name = 'event2';
plot_par.twindow = 0.75;
plot_par.lbnd4nenergy = 1;
plot_par.ubnd4nenergy = 2;

timevector_u = wave_u_full_matrix(:,1); %0:1/sample_rate: (length(wave_u_full)-1)*(1/sample_rate);
timevector_p = wave_p_full_matrix(:,1); %0:1/sample_rate: (length(wave_p_full)-1)*(1/sample_rate);

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


Lcut=[2.5 3 3.5 4 4.5 5];
Hcut=[1 1.5 2 2.5 3 3.5];

Lcut = 2;  %cutoff used with the low pass filter (i.e. (Hcut<Lcut)
Hcut = 1;  %cutoff used with the high pass filter (i.e. Hcut<Lcut)

% Lcut=[3 4 5 6 7 8 9 10];
% Hcut=[1 2 3 4 5 6 7 8];
% 
% Lcut=[2.5 5];
% Hcut=[1 2.5];

%Lcut = 10;
%Hcut = 1;
for i=1:length(Lcut);
   
    [wave_u_full, wave_p_full] = filter_waveform(wave_u_full, timevector_u,wave_p_full, timevector_p,Lcut(i),Hcut(i),diag_switch); 
    [norm_maxcorr,sep] = plot_seis_xcorr(plot_par,[timevector_u(:),wave_u_full(:)],[timevector_p(:),wave_p_full(:)],src_type,diag_switch);
    
    if diag_switch ==1 
        x_lim = get(gca,'xlim');
        plot(x_lim, known_sep*[1 1],'--')
    end
   
    start1 = length(sep) - floor(length(sep)/2);
    end1 = length(sep);
    meanR(i) = mean(norm_maxcorr(start1:end1,2));
    meanS(i) = mean(sep(start1:end1,2));   
   
end