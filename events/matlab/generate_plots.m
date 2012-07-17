function [hf, omega2,genplot]=generate_plots(wave_u, time_u, wave_p, time_p,title1, pow_xlim, auto_xlim, diag_switch)

%OUTPUTS:
% hf        [scalar] figure handle
% omega2    [scalar] omega2 value ??
% genplot   [structure] extra information for the waves
%               genplot.u.fdom = dominant frequency of unperturbed wave
%               genplot.u.power = power spectra [frequency,power] of unperturbed
%               genplot.p.fdom = dominant frequency of perturbed wave
%               genplot.p.power = power spectra [frequency,power] of perturbed
%
% diag_switch   [scalar]
%                   diag_switch ==1 => make plots and disp output to screen
%                   diag_switch ==0 => no plots and not output to screen




%sample_rate = 1/0.0017;  %% No this stuffs it up - you must compute!!!! 
sample_rate = 1/(time_u(2)-time_u(1)); % sample rate of perturbed wave is the same
hf = [];

if diag_switch ==1
    hf = figure
    % first lets plot the two waveforms on the same axes
    subplot(5,2,1:2)
    h2 = plot(time_p, wave_p, 'r' );
    hold on
    h1 = plot(time_u, wave_u, 'b');
    legend([h1, h2], {'Event 1', 'Event 2'});
    title(title1)
end

% Now lets plot the power spectra (to the Nyquist)
nu = length(wave_u);
nu2 = nextpow2(nu);
fft_wave_u=fft(wave_u,2^nu2);
Pyy_u = fft_wave_u.* conj(fft_wave_u) / 2^nu2;
f_u = sample_rate*(0:2^nu2/2)/2^nu2;
grad_wave_u =gradient(wave_u,1/sample_rate);
omega2_numerator_u = compute_energy(grad_wave_u,time_u);
omega2_denominator_u= compute_energy(wave_u,time_u);
f_dom_u = sqrt(omega2_numerator_u/(4*pi^2*omega2_denominator_u));
% Build output structure
genplot.u.fdom = f_dom_u;
if diag_switch ==1
    subplot(5,2,3)
    semilogx(f_u,Pyy_u(1:2^nu2/2+1),'b')
    xlabel('Frequency (Hz)')
    ylabel('Power')
    title(['computed dominant freq = ', num2str(f_dom_u)])
    disp('    ')
    disp('=======================================================')
    disp('   ')
    disp(title1)
    disp('   ')
    disp(['omega2_numerator_u = ', num2str(omega2_numerator_u)])
    disp(['omega2_denominator_u = ', num2str(omega2_denominator_u)])
    disp(['omega2 = ', num2str(omega2_numerator_u/omega2_denominator_u)])
    disp(['f_dom_u = ', num2str(f_dom_u)])
end


np = length(wave_p);
np2 = nextpow2(np);
fft_wave_p = fft(wave_p,2^np2);
Pyy_p = fft_wave_p.* conj(fft_wave_p) / 2^np2;
%f_p = 200*(0:2^np2/2)/2^np2;
f_p = sample_rate*(0:2^np2/2)/2^np2;
grad_wave_p =gradient(wave_p,1/sample_rate);
omega2_numerator_p = compute_energy(grad_wave_p,time_p);
omega2_denominator_p= compute_energy(wave_p,time_p);
f_dom_p = sqrt(omega2_numerator_p/(4*pi^2*omega2_denominator_p));
% Build output structure
genplot.p.fdom = f_dom_p;
omega2 = omega2_numerator_p/omega2_denominator_p;
if diag_switch ==1
    subplot(5,2,4)
    semilogx(f_p,Pyy_p(1:2^np2/2+1),'r')
    xlabel('Frequency (Hz)')
    ylabel('Power')
    title(['computed dominant freq = ', num2str(f_dom_p)])
    disp('    ')
    disp(['omega2_numerator_p = ', num2str(omega2_numerator_p)])
    disp(['omega2_denominator_p = ', num2str(omega2_denominator_p)])
    disp(['omega2 = ', num2str(omega2_numerator_p/omega2_denominator_p)])
    disp(['f_dom_p = ', num2str(f_dom_p)])
    disp('    ')
end

% Build output structure
genplot.u.power = [f_u',Pyy_u(1:2^nu2/2+1)];
% No we plot the power spectra again to some freq<Nyquist
if diag_switch ==1
    subplot(5,2,5)
    plot(f_u,Pyy_u(1:2^nu2/2+1),'b')
    xlabel('Frequency (Hz)')
    ylabel('Power')
    set(gca,'xlim',pow_xlim)
end

% Build output structure
genplot.p.power = [f_p',Pyy_p(1:2^np2/2+1)];
if diag_switch ==1
    subplot(5,2,6)
    plot(f_p,Pyy_p(1:2^nu2/2+1),'r')
    xlabel('Frequency (Hz)')
    ylabel('Power')
    set(gca,'xlim',pow_xlim)
end

% Now lets have a look at the auto-correlation
[autocorr_u,lags_u] = xcorr(wave_u);
if diag_switch ==1
    subplot(5,2,7)
    plot(lags_u/200,autocorr_u,'b');
    xlabel('Time (sec)')
    ylabel('A-xcorr')
end

[autocorr_p, lags_p] = xcorr(wave_p);
if diag_switch ==1
    subplot(5,2,8)
    plot(lags_p/200,autocorr_p,'r');
    xlabel('Time (sec)')
    ylabel('A-xcorr')

    % Now lets zoom into the auto-correlation
    subplot(5,2,9)
    plot(lags_u/200,autocorr_u,'b');
    set(gca,'xlim',auto_xlim)
    xlabel('Time (sec)')
    ylabel('A-xcorr')

    subplot(5,2,10)
    plot(lags_p/200,autocorr_p,'r');
    set(gca,'xlim',auto_xlim)
    xlabel('Time (sec)')
    ylabel('A-xcorr')

end