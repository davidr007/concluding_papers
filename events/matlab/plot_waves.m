function [] = plot_waves(fname_u, fname_p, sample_rate)

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
 
    sample_rate = 200
end


wave_u_full = load(fname_u);
wave_p_full = load(fname_p);




%Now setup the plotting structure for PLOT_SEIS_XCORR
plot_par.trans1 = 0.05;    % horizontal translation of event1(:,1) [time]
plot_par.trans2 =  0%0.585;% 0.006;         % horizontal translation of event2(:,1) [time]
plot_par.start_arrivwin = 29.7;
plot_par.end_arrivwin = 32.7;
plot_par.start_codawin = 40;
plot_par.end_codawin = 43;
plot_par.event1name = 'event1';
plot_par.event2name = 'event2';

timevector_u = 0:1/sample_rate: (length(wave_u_full)-1)*(1/sample_rate);
timevector_p = 0:1/sample_rate: (length(wave_p_full)-1)*(1/sample_rate);

figure
subplot(2,1,1)
h1 = plot(timevector_u, wave_u_full, 'b');
xlabel('time(sec)')
title('Event 1')
subplot(2,1,2)
h2 = plot(timevector_p, wave_p_full, 'r' );
xlabel('time(sec)')
title('Event 1')


% Now we need to mean shift both waveforms so that they are around 0
wave_u_full = wave_u_full - median(wave_u_full);
wave_p_full = wave_p_full - median(wave_p_full);
generate_plots(wave_u_full, timevector_u, wave_p_full, timevector_p, ...
    'original waveforms (step 1)', [0,2],[-2500/200,2500/200] )



% Now we need High Pass filter the signals to remove the very low freq
wave_u_full = apply_butter(wave_u_full, 5,'high', 'time',1, 200);
wave_p_full = apply_butter(wave_p_full, 5,'high', 'time',1, 200);
generate_plots(wave_u_full, timevector_u, wave_p_full, timevector_p, ...
    'High Pass 1Hz (step 2)', [0,40],[-0.25, 0.25] )



% Now we need to mean shift both waveforms so that they are around 0
wave_u_full = wave_u_full - median(wave_u_full);
wave_p_full = wave_p_full - median(wave_p_full);
generate_plots(wave_u_full, timevector_u, wave_p_full, timevector_p, ...
    'DC Removal of mean (step 3)',[0,40],[-0.25, 0.25] )

%Now we low pass filter the signals to remove all of the very high freq
wave_u_full = apply_butter(wave_u_full, 5,'low', 'time',10, 200);
wave_p_full = apply_butter(wave_p_full, 5,'low', 'time',10, 200);
generate_plots(wave_u_full, timevector_u, wave_p_full, timevector_p, ...
    'Low Pass 10Hz (step 4)',[0,40],[-0.25, 0.25] )



% Now we need to mean shift both waveforms so that they are around 0
wave_u_full = wave_u_full - median(wave_u_full);
wave_p_full = wave_p_full - median(wave_p_full);
generate_plots(wave_u_full, timevector_u, wave_p_full, timevector_p, ...
'DC Removal of mean (step 5)', [0,40],[-0.25, 0.25] )

% 
% Finally, we scale both waveforms
wave_u_full = wave_u_full/max(abs(wave_u_full));
wave_p_full = wave_p_full/max(abs(wave_p_full));
generate_plots(wave_u_full, timevector_u, wave_p_full, timevector_p, ...
    'Scale Both Waveforms (step 6)', [0,40],[-0.25, 0.25] )

%plot_seis_xcorr(plot_par,[timevector_u(:),wave_u_full(:)],[timevector_p(:),wave_p_full(:)])

% function LOC_generate_plot(wave_u, time_u, wave_p, time_p,title1, pow_xlim, auto_xlim)
% figure
% sample_rate = 200;
% 
% % first lets plot the two waveforms on the same axes
% subplot(5,2,1:2)
% h2 = plot(time_p, wave_p, 'r' );
% hold on
% h1 = plot(time_u, wave_u, 'b'); 
% legend([h1, h2], {'Event 1', 'Event 2'});
% title(title1)
% 
% % Now lets plot the power spectra (to the Nyquist)
% 
% 
% 
% nu = length(wave_u);
% nu2 = nextpow2(nu);
% fft_wave_u=fft(wave_u,2^nu2);
% subplot(5,2,3)
% Pyy_u = fft_wave_u.* conj(fft_wave_u) / 2^nu2;
% f_u = 200*(0:2^nu2/2)/2^nu2;
% plot(f_u,Pyy_u(1:2^nu2/2+1),'b')
% xlabel('Frequency (Hz)')
% ylabel('Power')
% grad_wave_u =gradient(wave_u,1/sample_rate);
% omega2_numerator_u = compute_energy(grad_wave_u,time_u);
% omega2_denominator_u= compute_energy(wave_u,time_u);
% f_dom_u = sqrt(omega2_numerator_u/(4*pi^2*omega2_denominator_u));
% title(['computed dominant freq = ', num2str(f_dom_u)])
% disp('    ')
% disp('=======================================================')
% disp('   ')
% disp(title1)
% disp('   ')
% disp(['omega2_numerator_u = ', num2str(omega2_numerator_u)])
% disp(['omega2_denominator_u = ', num2str(omega2_denominator_u)])
% disp(['omega2 = ', num2str(omega2_numerator_u/omega2_denominator_u)])
% disp(['f_dom_u = ', num2str(f_dom_u)])
% 
% 
% np = length(wave_p);
% np2 = nextpow2(np);
% fft_wave_p = fft(wave_p,2^np2);
% subplot(5,2,4)
% Pyy_p = fft_wave_p.* conj(fft_wave_p) / 2^np2;
% f_p = 200*(0:2^np2/2)/2^np2;
% plot(f_p,Pyy_p(1:2^np2/2+1),'r')
% xlabel('Frequency (Hz)')
% ylabel('Power')
% grad_wave_p =gradient(wave_p,1/sample_rate);
% omega2_numerator_p = compute_energy(grad_wave_p,time_p);
% omega2_denominator_p= compute_energy(wave_p,time_p);
% f_dom_p = sqrt(omega2_numerator_p/(4*pi^2*omega2_denominator_p));
% title(['computed dominant freq = ', num2str(f_dom_p)])
% disp('    ')
% disp(['omega2_numerator_p = ', num2str(omega2_numerator_p)])
% disp(['omega2_denominator_p = ', num2str(omega2_denominator_p)])
% disp(['omega2 = ', num2str(omega2_numerator_p/omega2_denominator_p)])
% disp(['f_dom_p = ', num2str(f_dom_p)])
% disp('    ')
% 
% 
% % No we plot the power spectra again to some freq<Nyquist
% subplot(5,2,5)
% plot(f_u,Pyy_u(1:2^nu2/2+1),'b')
% xlabel('Frequency (Hz)')
% ylabel('Power')
% set(gca,'xlim',pow_xlim)
% 
% subplot(5,2,6)
% plot(f_p,Pyy_p(1:2^nu2/2+1),'r')
% xlabel('Frequency (Hz)')
% ylabel('Power')
% set(gca,'xlim',pow_xlim)
% 
% % Now lets have a look at the auto-correlation
% [autocorr_u,lags_u] = xcorr(wave_u);
% subplot(5,2,7)
% plot(lags_u/200,autocorr_u,'b');
% xlabel('Time (sec)')
% ylabel('A-xcorr')
% 
% [autocorr_p, lags_p] = xcorr(wave_p);
% subplot(5,2,8)
% plot(lags_p/200,autocorr_p,'r');
% xlabel('Time (sec)')
% ylabel('A-xcorr')
% 
% % Now lets zoom into the auto-correlation
% subplot(5,2,9)
% plot(lags_u/200,autocorr_u,'b');
% set(gca,'xlim',auto_xlim)
% xlabel('Time (sec)')
% ylabel('A-xcorr')
% 
% subplot(5,2,10)
% plot(lags_p/200,autocorr_p,'r');
% set(gca,'xlim',auto_xlim)
% xlabel('Time (sec)')
% ylabel('A-xcorr')