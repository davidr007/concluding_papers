% This script is used to check the St Helens Data with my CWI software
% Note that we work with stacks09(:,1:24) as this is the same as what
% was used in Roel's in press paper


plot_par.trans1 = 0;    % horizontal translation of event1(:,1) [time]
plot_par.trans2 = 0;         % horizontal translation of event2(:,1) [time]
plot_par.start_arrivwin = -0.05;
plot_par.end_arrivwin = 1.5;
plot_par.start_codawin = 2.2;
plot_par.end_codawin = 3.7;
plot_par.parrival = 0;
plot_par.event1name = '0 hour';
plot_par.event2name = 'later hour';
plot_par.save = 1;
plot_par.u_name = '0 hour';
plot_par.p_name = '1 hour';
plot_par.twindow = 5;
plot_par.lbnd4nenergy = 0.02;
plot_par.ubnd4nenergy = 0.14;






Lcut = [2 4 6];
Hcut = [4 6 10];

load data_V6.mat
%  loads:   abs_peaks08 abs_peaks09 abs_peaks10
%           stacks08    stacks09    stacks10    time
clear norm_maxcorr  sep
ind = 1; % used sometimes to cut the waveform at start


f1 = figure
f2 = figure
for j = 1:length(Lcut)
    for i = 2:24
        time_u = time(ind:end)';
        wave_u_tmp = stacks09(ind:end,1);
        time_p = time(ind:end)';
        wave_p_tmp = stacks09(ind:end,i);
        
        [wave_u, wave_p] = filter_waveform(wave_u_tmp, time_u,wave_p_tmp, time_p,Lcut(j),Hcut(j),0);
        
        %[wave_u] = filter_waveform(wave_u_tmp, time_u,[],[],Lcut(j),Hcut(j),0);
        %[wave_p] = filter_waveform(wave_p_tmp, time_u,[],[],Lcut(j),Hcut(j),0);
       
        [tmp_norm_maxcorr,tmp_sep, tmp_sep2, tmp_xcorr2] = seis_xcorr(wave_u, wave_p, time_u, time_p,plot_par.twindow,'i',2);
        %[tmp_norm_maxcorr,tmp_sep] = plot_seis_xcorr(plot_par,[time_u wave_u],[time_u wave_u],2)
        if i ==1
            [n m] = length(tmp_norm_maxcorr);
            norm_maxcorr = zeros(n,24);
            sep = zeros(n,24);
            norm_maxcorr(:,1) = tmp_norm_maxcorr(:,1);
            sep(:,1) = tmp_sep(:,1);
        end
        norm_maxcorr(:,i) =  tmp_norm_maxcorr(:,2);
        sep(:,i) = tmp_sep(:,2);
        % close all
    end

    figure(f1)
    subplot(1,3,j)
    l1 = plot([1:23], norm_maxcorr(1,2:end),'k','linewidth',2), hold on
    l2 = plot([1:23], norm_maxcorr(2,2:end),'k--','linewidth',2)
    l3 = plot([1:23], norm_maxcorr(3,2:end),'-.k','linewidth',2)
    l4 = plot([1:23], norm_maxcorr(4,2:end),'k:','linewidth',2)
    legend([l1,l2,l3,l4],{'0-5 s','5-10 s','10-15 s','15-20 s'},'location','SouthWest')
    xlabel('age (h)')
    ylabel('max xcorr')
    grid on
    title([num2str(Lcut(j)), ' to ', num2str(Hcut(j)), ' Hz'])
    set(gca,'ylim',[0.4 1])

    figure(f2)
    subplot(1,3,j)
    l1 = plot([1:23], sep(1,2:end),'k','linewidth',2), hold on
    l2 = plot([1:23], sep(2,2:end),'k--','linewidth',2)
    l3 = plot([1:23], sep(3,2:end),'-.k','linewidth',2)
    l4 = plot([1:23], sep(4,2:end),'k:','linewidth',2)
    legend([l1,l2,l3,l4],{'0-5 s','5-10 s','10-15 s','15-20 s'},'location','SouthEast')
    xlabel('age (h)')
    ylabel('separation (m)')
    grid on
    title([num2str(Lcut(j)), ' to ', num2str(Hcut(j)), ' Hz'])
    set(gca,'ylim',[0 120])
end




