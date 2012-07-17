function [xcorrStruc, h1, sepStruc, h2] = plot_all_stations3D(dname_u,dname_p, ...
    sample_rate,known_sep, Lcut, Hcut,diag_vec,plot_par)

% PLOT_ALL_STATIONS3D is used to plot the CWI source separation results 
% from the three channels of one station from a 3D simulation.
% Note that unlike the 2D version of this function we do not do all
% stations together (with teh 2D synthetics all stations are in the single
% file, whereas for the 3D synthetics there is a separate file for each
% station. 
% Two files are provided representing  
% the waveforms for two events as produced by pmcl3d. Two figures are
% returned showing (1) the cross correlations and (2) CWI computed and known
% separations at each of the stations. The cross correlations and
% separations are returned in XCORRMAT and SEPMAT.
%
% CAUTION: - this function is for the 2D acoustic modelled waveforms only.
% It will also break if more than 11 stations are used. 
%
%Inputs:
% dname_u   [string] directory for location of station files for the unperturbed event 
%           of interest. Typically, the directory will contain 36 station
%           files. Each file must have 4 columns. Note that there is not 
%           always a space between the columns (see load technique in code).
%           Also note that the first column does not contain time - it contains
%           the reading number. A time vector must be created using the sample
%           rate. 
% dname_p   [string] station directory for the perturbed event. Same format as 
%           dname_u.
% sample_rate [scalar] the sample rate (IN SECONDS) used to create the timevector
%               e.g. sample_rate=0.0017 is used for many sims
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
%xcorrStruc [Structure] structure containing normalised maximum cross
%           correlation matrices as produced by plot_seis_xcorr. There
%           is a field for each station.
% h1        [Scalar] figure handle for Max-xcorr plots
%sepMat     [Structure] structure containing separation matrices
%           as produced by plot_seis_xcorr. There is a field for each 
%           station.
% h2        [Scalar] figure handle for Sep plots
%
% David Robinson
% 28 March 2006

stations = [101, 102, 103, 104, 105, 106, 107, 108, 109,...
            201, 202, 203, 204, 205, 206, 207, 208, 209, ...
            301, 302, 303, 304, 305, 306, 307, 308, 309, ...
            401, 402, 403, 404, 405, 406, 407, 408, 409]; 
% 

for i =1:length(stations) % loop over all of the stations computing CWI estimates
    [xcorrStruc.(['station',num2str(stations(i))]), ...
        sepStruc.(['station',num2str(stations(i))])] = ...
                LOC_plot4allchannels(stations(i),dname_u,dname_p, ...
                    sample_rate,known_sep, Lcut, Hcut,diag_vec,plot_par);
end


h1 = [];
h2 = [];

%% Local function to do CWI for each channel of a single station
function [xcorrStruc,sepStruc] = LOC_plot4allchannels(station_name,dname_u,dname_p, ...
    sample_rate,known_sep, Lcut, Hcut,diag_vec,plot_par)

stat_str = ['station',num2str(station_name),'.txt'];

fid = fopen([dname_u,'/',stat_str]);
event1 = fscanf(fid,'%i %e %e %e', [4 Inf]);
event1 = event1';
fclose(fid);
wave_u_full = event1(:,2:4); %

fid = fopen([dname_p,'/',stat_str]);
event2 = fscanf(fid,'%i %e %e %e', [4 Inf]);
event2 = event2';
fclose(fid);
wave_p_full = event2(:,2:4);

% src_type = plot_par.src_type;    %(1 => double couple source in a 2D elastic medium
% %(2 => double couple source in a 3D elastic medium
% %(3 => 3D double couple source [delta source variation])
timevector_u = 0:sample_rate: (length(wave_u_full)-1)*sample_rate;
timevector_p = 0:sample_rate: (length(wave_p_full)-1)*sample_rate;

% Do the work
for i =1:3
    [wave_u_full(:,i), wave_p_full(:,i)] = filter_waveform(wave_u_full(:,i), timevector_u,wave_p_full(:,i), timevector_p,Lcut,Hcut,diag_vec(1));
    [xcorrStruc.(['chan',num2str(i)]),sepStruc.(['chan',num2str(i)])] = plot_seis_xcorr(plot_par,[timevector_u(:),wave_u_full(:,i)],[timevector_p(:),wave_p_full(:,i)],diag_vec(2));
end

%     if diag_vec(2)==1
%         x_lim = get(gca,'xlim');
%         plot(x_lim, known_sep*[1 1],'--')
%     end










%
% % Initialise vectors to store ylimits for plots
% ylim_Xcorr = [NaN,NaN];
% ylim_Sep = [NaN,NaN];
%
% for i=1:11;
%
%     wave_u_full = wave_u_full_matrix(:,i+1);
%     wave_p_full = wave_p_full_matrix(:,i+1);
%
%     [wave_u_full, wave_p_full] = filter_waveform(wave_u_full, timevector_u,wave_p_full, timevector_p,Lcut,Hcut,diag_vec(1));
%     [norm_maxcorr,sep] = plot_seis_xcorr(plot_par,[timevector_u(:),wave_u_full(:)],[timevector_p(:),wave_p_full(:)],diag_vec(2));
%     if diag_vec(2)==1
%         x_lim = get(gca,'xlim');
%         plot(x_lim, known_sep*[1 1],'--')
%     end
%
%     xcorrStruc.(['station', num2str(i)]) = norm_maxcorr;
%     sepStruc.(['station', num2str(i)]) = sep;
%
%     % build ylimits for plotting as you go
%     ylim_Xcorr(1) = min([ylim_Xcorr(1), min(xcorrStruc.(['station', num2str(i)])(:,2))]);
%     ylim_Xcorr(2) = max([ylim_Xcorr(2), max(xcorrStruc.(['station', num2str(i)])(:,2))]);
%     ylim_Sep(1) = min([known_sep-20, ylim_Sep(1), min(sepStruc.(['station', num2str(i)])(:,2))]);
%     ylim_Sep(2) = max([known_sep+20, ylim_Sep(2), max(sepStruc.(['station', num2str(i)])(:,2))]);
%
% end
%
% if diag_vec(3) ==1
%
%     for i=1:11
%         % create plots for the cross correlations
%         h1 = figure; % figure for cross correlations
%         s1(i) = subplot(6,2,i)
%         plot(xcorrStruc.(['station', num2str(i)])(:,1), xcorrStruc.(['station', num2str(i)])(:,2))
%         %xlabel('time (sec)')
%         %ylabel('X-corr max')
%         %title(['Station ', num2str(i)])
%         hleg1 = legend(['Station ', num2str(i)])
%         set(hleg1,'Location','SouthEast', 'box','off')
%         tmp = get(hleg1,'children')
%         set(tmp(2), 'lineStyle','none')
%         % set uniform ylimits for all sub-plots
%         set(s1(i), 'ylim', ylim_Xcorr)
%
%         % create plots for the separations
%         h2 = figure;
%         s2(i)=subplot(6,2,i)
%         plot(sepStruc.(['station', num2str(i)])(:,1), sepStruc.(['station', num2str(i)])(:,2))
%         %xlabel('time (sec)')
%         %ylabel('separation (m)')
%         %title(['Station ', num2str(i)])
%         hleg2 = legend(['Station ', num2str(i)])
%         set(hleg2, 'box','off')
%         tmp = get(hleg2,'children')
%         set(tmp(2), 'lineStyle','none')
%
%         x_lim = get(gca,'xlim');
%         hold on, plot(x_lim, known_sep*[1 1],'--')
%         % set uniform ylimits for all sub-plots
%         set(s2(i), 'ylim', ylim_Sep);
%     end
%
%     if diag_vec(4)==1 & diag_vec(3)==1
%         % Add details for the figure
%         figure(h1)
%         s1(12) = subplot(6,2,12)
%         axis off
%         text(0.1, 0.6, {'Maximum X-corr:'; [fname_u(end-17:end-4), ' and ', ...
%                             fname_p(end-17:end-4)]; [num2str(Hcut), ...
%                             ' to ', num2str(Lcut), ' Hz']; ['known separation = ', ...
%                             num2str(known_sep), ' m '] });
%
%         % Add details for the figure
%         figure(h2)
%         s2(12) = subplot(6,2,12)
%         axis off
%         text(0.2, 0.7, {'Separation:'; [fname_u(end-17:end-4), ' and ', ...
%                             fname_p(end-17:end-4)]; [num2str(Hcut), ...
%                             ' to ', num2str(Lcut), ' Hz']; ['known separation = ', ...
%                             num2str(known_sep), ' m ']});
%     end
%
% end