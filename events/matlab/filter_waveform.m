function [waveout_p, waveout_u,extras] = filter_waveform(wavein_u, time_u,wavein_p, time_p,Lcut,Hcut,diag_switch )
% INPUTS:
% wave_in = wave to be filtered
% time = time vector for wave_in
% diag_switch   [scalar]
%                   1 => draw plots
%                   0 => do not draw plots (compute values only)
%
% Outputs:
% waveout_p     [vector mx1]: the perturbed wave after filtering (note time data not
%               included: - this must be re-built from sample_rate
% waveout_u     [vector mx1]: the un-perturbed or reference wave after filtering 
%               (note time data not included: - this must be re-built from
%               sample_rate)
% extras        [structure] contains optional extra outputs
%                 **** For the initial (input) waves
%                   extras.init.u.fdom = dominant frequency of unperturbed wave
%                   extras.init.u.power = power spectra [frequency,power] of unperturbed
%                   extras.init.p.fdom = dominant frequency of perturbed wave
%                   extras.init.p.power = power spectra [frequency,power] of perturbed
%                 **** For the end (after filtering) waves
%                   extras.end.u.fdom = dominant frequency of unperturbed wave
%                   extras.end.u.power = power spectra [frequency,power] of unperturbed
%                   extras.end.p.fdom = dominant frequency of perturbed wave
%                   extras.end.p.power = power spectra [frequency,power] of perturbed

sample_rate = 1/(time_u(2)-time_u(1));

    % Now we need to mean shift the waveform so it is around 0
    waveout_u = wavein_u - median(wavein_u);
    waveout_p = wavein_p - median(wavein_p);
    [hf_loc_init, omega2_init, genplots_init] = generate_plots(wavein_u, time_u, wavein_p, time_p, ...
            'original waveforms (step 1)', [0,2],[-2500/200,2500/200], 0);
    extras.init = genplots_init;
    if diag_switch ~=1
        close(hf_loc_init)
    end
%     
    %Now we should taper it to zero at the end
    if length(wavein_p) > 1000
        extras.taperlength = 201;
    else
        extras.taperlength = ceil(length(wavein_p)/20);
    end
    taper_u = ztaper_tsignal(length(wavein_u),extras.taperlength, 1);
    waveout_u = waveout_u.*taper_u;      
    taper_p = ztaper_tsignal(length(wavein_p), extras.taperlength, 1);
        
    try, waveout_p = waveout_p.*taper_p; catch, end %try catch so that function can be used to filter a single waveform
    
    if diag_switch ==1
        generate_plots(waveout_u, time_u, waveout_p, time_p, ...
            'Tapered to zero at each end (step 2)', [0,2],[-2500/200,2500/200], diag_switch )
    end

    % Now we need High Pass filter the signal to remove the very low freq
    waveout_u = apply_butter(waveout_u, 5,'high', 'time',Hcut, sample_rate);
    waveout_p = apply_butter(waveout_p, 5,'high', 'time',Hcut, sample_rate);
    
    if diag_switch ==1
        generate_plots(waveout_u, time_u, waveout_p, time_p, ...
            ['High Pass ',num2str(Hcut),'Hz (step 3)'], [0,40],[-0.25, 0.25],  diag_switch)
    end


    % Now we need to mean shift the waveform so that it is around 0
    waveout_u = waveout_u - median(waveout_u);
    waveout_p = waveout_p - median(waveout_p);
    
    if diag_switch ==1
        generate_plots(waveout_u, time_u, waveout_p, time_p, ...
            'DC Removal of mean (step 4)',[0,40],[-0.25, 0.25], diag_switch )
    end
    
      %Now we should taper it to zero at the end - we need to do this again
  %because of the second DC shift
    taper2_u = ztaper_tsignal(length(waveout_u), 50, 1);
    waveout_u = waveout_u.*taper2_u;
    taper2_p = ztaper_tsignal(length(waveout_p), 50, 1);
    try, waveout_p = waveout_p.*taper2_p; catch, end
    
    if diag_switch ==1
        generate_plots(waveout_u, time_u, waveout_p, time_p, ...
            'Tapered to zero at each end (step 2)', [0,2],[-2500/200,2500/200], diag_switch );
    end
        
    %Now we low pass filter the signal to remove all of the very high freq

    waveout_u = apply_butter(waveout_u, 5,'low', 'time',Lcut, sample_rate);
    waveout_p = apply_butter(waveout_p, 5,'low', 'time',Lcut, sample_rate);
    if diag_switch ==1
        generate_plots(waveout_u, time_u, waveout_p, time_p, ...
            ['Low Pass ',num2str(Lcut),'Hz (step 5)'],[0,40],[-0.25, 0.25], diag_switch )
    end


    % Now we need to mean shift the waveform so that it is around 0
    waveout_u = waveout_u - median(waveout_u);
    waveout_p = waveout_p - median(waveout_p);
    
    if diag_switch ==1
        generate_plots(waveout_u, time_u, waveout_p, time_p, ...
            'DC Removal of mean (step 6)', [0,40],[-0.25, 0.25], diag_switch )
    end
%     
    %
    % Finally, we scale both waveforms
    waveout_u = waveout_u/max(abs(waveout_u));
    waveout_p = waveout_p/max(abs(waveout_p));
    [hf_loc_end, omega2_end, genplots_end] = generate_plots(waveout_u, time_u, waveout_p, time_p, ...
            'Scale Both Waveforms (step 7)', [0,40],[-0.25, 0.25],0 );
    extras.end = genplots_end;
    if diag_switch ~=1
        close(hf_loc_end)
    end