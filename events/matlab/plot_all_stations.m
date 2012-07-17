function [xcorrStruc, h1, sepStruc, h2,extras, omega2Struc] = plot_all_stations(fname_u,fname_p,known_sep, Lcut, Hcut,diag_vec,plot_par)

% PLOT_ALL_STATIONS is used to plot the CWI source separation results 
% from the 11 surface stations. Two files are provided representing  
% the waveforms for two events as produced by ac2d_for.f90. Two figures are
% returned showing (1) the cross correlations and (2) CWI computed and known
% separationss at each of the stations. The cross correlations and
% separations are returned in XCORRMAT and SEPMAT.
%
% CAUTION: - this function is for the 2D acoustic modelled waveforms only.
% It will also break if more than 11 stations are used. 
%
%Inputs:
% fname_u   [string] waveforms for the unperturbed event. This file 
%           must have 12 columns. The first column is time, the second
%           to twelth columns are waveforms measured at different receivers
%           (i.e. the file usually comes from ac2d_for.f90)
% fname_p   [string] waveforms for the perturbed event. Same format as 
%           fname_u.
% known_sep [scalar] the known separation between events
% Lcut      [scalar] cutoff used with the low pass filter (i.e. Hcut<Lcut)
% Hcut      [scalar] cutoff used with the high pass filter (i.e. Hcut<Lcut)
% diag_vec  [vector] vector containg diag_switch values. Note that 1 =>
%           diagnostic plots produced and 0 => plots not produced.
%               diag_vec(1) = diag_switch for FILTER_WAVEFORM
%               diag_vec(2) = diag_switch for PLOT_SEIS_XCORR
%               diag_vec(3) = diag_swicth for parent (PLOT_ALL_STATIONS)
%               diag_vec(4) = diagnostic writing on PLOT_ALL_STATIONS
%                             figure if diag_switch(3)==1
%
% OUTPUTS
%xcorrStruc     [Structure] structure containing normalised maximum cross
%               correlation matrices as produced by plot_seis_xcorr. There
%               is a field for each station.
% h1            [Scalar] figure handle for Max-xcorr plots
%sepMat         [Structure] structure containing separation matrices
%               as produced by plot_seis_xcorr. There is a field for each 
%               station.
% h2            [Scalar] figure handle for Sep plots
% extras        [Structure] contains extra information about waveforms (see 
%               FILTER_WAVEFORM for a complete description.
% omega2Struc   [Structure] structure containing the omega^2 estimates used
%               to compute separation. There is a field for each station.
%
% David Robinson
% 28 March 2006


h1 = [];
h2 = [];

%load the files
wave_u_full_matrix_tmp = load(fname_u);
wave_p_full_matrix_tmp = load(fname_p);
fname_u_short = fname_u(end-17:end-4);
fname_p_short = fname_p(end-17:end-4);

%check that both files have the required number of columns
[n m] = size(wave_u_full_matrix_tmp);
if m==12
    wave_u_full_matrix = wave_u_full_matrix_tmp;
elseif m==1 % we must re-arrange to a 12 column matrix due to problem with output from ac2d_for (i.e. data in 1 column vector
    wave_u_full_matrix = zeros(length(wave_u_full_matrix_tmp)/12,12);
    count =1;
    for i = 1:length(wave_u_full_matrix_tmp)/12
        wave_u_full_matrix(i,1:12) = wave_u_full_matrix_tmp(count:count+11)';
        count = count+12;
    end
else
    error('Incorrect number of columns in fname_u')
end


[n m] = size(wave_p_full_matrix_tmp);
if m==12
    wave_p_full_matrix = wave_p_full_matrix_tmp;
elseif m==1 % we must re-arrange to a 12 column matrix due to problem with output from ac2d_for (i.e. data in 1 column vector
    wave_p_full_matrix = zeros(length(wave_p_full_matrix_tmp)/12,12);
    count =1;
    for i = 1:length(wave_p_full_matrix_tmp)/12
        wave_p_full_matrix(i,1:12) = wave_p_full_matrix_tmp(count:count+11)';
        count = count+12;
    end
else
    error('Incorrect number of columns in fname_p') 
end



src_type = plot_par.src_type;    %(1 => double couple source in a 2D elastic medium
                %(2 => double couple source in a 3D elastic medium
                %(3 => 3D double couple source [delta source variation])
%Now setup the plotting structure for PLOT_SEIS_XCORR


timevector_u_tmp = wave_u_full_matrix(:,1); %0:1/sample_rate: (length(wave_u_full)-1)*(1/sample_rate);
timevector_p_tmp = wave_p_full_matrix(:,1); %0:1/sample_rate: (length(wave_p_full)-1)*(1/sample_rate);

% Initialise vectors to store ylimits for plots
ylim_Xcorr = [NaN,NaN];
ylim_Sep = [NaN,NaN];

for i=1:11;
   
    wave_u_full_tmp = wave_u_full_matrix(:,i+1);
    wave_p_full_tmp = wave_p_full_matrix(:,i+1);
    alignval_u = plot_par.firstarrival.(fname_u_short).(['station',num2str(i)]);
    alignval_p = plot_par.firstarrival.(fname_p_short).(['station',num2str(i)]);
    
    [timevector_u,wave_u_full,timevector_p,wave_p_full,aligndetails] = align_waveforms(timevector_u_tmp,wave_u_full_tmp,timevector_p_tmp,wave_p_full_tmp,alignval_u,alignval_p);
    [wave_u_full, wave_p_full,extras_tmp] = filter_waveform(wave_u_full, timevector_u,wave_p_full, timevector_p,Lcut,Hcut,diag_vec(1)); 
   % figure, plot(timevector_u(:),wave_u_full(:),'r'), hold on, plot(timevector_p(:),wave_p_full(:))
    [norm_maxcorr,sep, omega2] = plot_seis_xcorr(plot_par,[timevector_u(:),wave_u_full(:)],[timevector_p(:),wave_p_full(:)],diag_vec(2));
    if diag_vec(2)==1
        x_lim = get(gca,'xlim');
        plot(x_lim, known_sep*[1 1],'--')
    end
    
    xcorrStruc.(['station', num2str(i)]) = norm_maxcorr;
    sepStruc.(['station', num2str(i)]) = sep;
    omega2Struc.(['station', num2str(i)]) = omega2;
    extras.(['station', num2str(i)]) = extras_tmp;
    
    % build ylimits for plotting as you go
    ylim_Xcorr(1) = min([ylim_Xcorr(1), min(xcorrStruc.(['station', num2str(i)])(:,2))]);
    ylim_Xcorr(2) = max([ylim_Xcorr(2), max(xcorrStruc.(['station', num2str(i)])(:,2))]);
    ylim_Sep(1) = min([known_sep-20, ylim_Sep(1), min(sepStruc.(['station', num2str(i)])(:,2))]);
    ylim_Sep(2) = max([known_sep+20, ylim_Sep(2), max(sepStruc.(['station', num2str(i)])(:,2))]);

end

if diag_vec(3) ==1

    for i=1:11
        % create plots for the cross correlations
        h1 = figure; % figure for cross correlations
        s1(i) = subplot(6,2,i)
        plot(xcorrStruc.(['station', num2str(i)])(:,1), xcorrStruc.(['station', num2str(i)])(:,2))
        %xlabel('time (sec)')
        %ylabel('X-corr max')
        %title(['Station ', num2str(i)])
        hleg1 = legend(['Station ', num2str(i)])
        set(hleg1,'Location','SouthEast', 'box','off')
        tmp = get(hleg1,'children')
        set(tmp(2), 'lineStyle','none')
        % set uniform ylimits for all sub-plots
        set(s1(i), 'ylim', ylim_Xcorr)
               
        % create plots for the separations
        h2 = figure;
        s2(i)=subplot(6,2,i)
        plot(sepStruc.(['station', num2str(i)])(:,1), sepStruc.(['station', num2str(i)])(:,2))
        %xlabel('time (sec)')
        %ylabel('separation (m)')
        %title(['Station ', num2str(i)])
        hleg2 = legend(['Station ', num2str(i)])
        set(hleg2, 'box','off')
        tmp = get(hleg2,'children')
        set(tmp(2), 'lineStyle','none')
        
        x_lim = get(gca,'xlim');
        hold on, plot(x_lim, known_sep*[1 1],'--')
        % set uniform ylimits for all sub-plots
        set(s2(i), 'ylim', ylim_Sep);
    end
    
    if diag_vec(4)==1 & diag_vec(3)==1
        % Add details for the figure
        figure(h1)
        s1(12) = subplot(6,2,12)
        axis off
        text(0.1, 0.6, {'Maximum X-corr:'; [fname_u(end-17:end-4), ' and ', ...
                            fname_p(end-17:end-4)]; [num2str(Hcut), ...
                            ' to ', num2str(Lcut), ' Hz']; ['known separation = ', ...
                            num2str(known_sep), ' m '] });

        % Add details for the figure
        figure(h2)
        s2(12) = subplot(6,2,12)
        axis off
        text(0.2, 0.7, {'Separation:'; [fname_u(end-17:end-4), ' and ', ...
                            fname_p(end-17:end-4)]; [num2str(Hcut), ...
                            ' to ', num2str(Lcut), ' Hz']; ['known separation = ', ...
                            num2str(known_sep), ' m ']});
    end
    
end