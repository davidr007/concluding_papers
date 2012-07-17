function [norm_maxcorr,sep, omega2] = plot_seis_xcorr(plot_par,event_u,event_p,diag_switch)

% INPUTS:
% diag_switch   [scalar]
%                   1 => draw plots
%                   0 => do not draw plots (compute values only)


if isempty(plot_par)
% important plotting parameters
    plot_par.trans1 = 0;    % horizontal translation of event1(:,1) [time]
    plot_par.trans2 = 0;% 0.006;         % horizontal translation of event2(:,1) [time]
    plot_par.start_arrivwin = 30.5;
    plot_par.end_arrivwin = 34.5;
    plot_par.start_codawin = 40;
    plot_par.end_codawin = 44;
    plot_par.parrival = 29.5;
    plot_par.event1name = 'event1';
    plot_par.event2name = 'event2';
    plot_par.save = 1;
    plot_par.u_name = 'receiver100001.txt';
    plot_par.p_name = 'receiver100002.txt';
    plot_par.twindow = 10;
    plot_par.lbnd4nenergy = -1;
    plot_par.ubnd4nenergy = -0.1;
end
trans1=plot_par.trans1;    % horizontal translation of event1(:,1) [time]
trans2=plot_par.trans2;% 0.006;         % horizontal translation of event2(:,1) [time]
start_arrivwin=plot_par.start_arrivwin;
end_arrivwin=plot_par.end_arrivwin;
start_codawin=plot_par.start_codawin;
end_codawin=plot_par.end_codawin;
event1name=plot_par.event1name;
event2name=plot_par.event2name;
twindow = plot_par.twindow;
%src_type = plot_par.src_type;

% Now we need to make sure the waves are lined up and that they are the same length
[indcut_u,closest,mind] = find_closest(event_u(:,1),event_u(1,1)+ trans1,'euclidean');
[indcut_p,closest,mind] = find_closest(event_u(:,1),event_u(1,1)+ trans2,'euclidean');
wave_u = [zeros(size(event_u(indcut_p:end,2))), event_u(indcut_p:end,2), event_u(indcut_p:end,2)];
wave_p = [zeros(size(event_p(indcut_u:end,2))), event_p(indcut_u:end,2), event_p(indcut_u:end,2)];
ndiff = abs(size(wave_u,1)-size(wave_p,1));
cutlength_lower = floor(ndiff/2);
cutlength_upper = ceil(ndiff/2);
np = size(wave_p,1);
nu = size(wave_u,1);
if size(wave_u,1)>size(wave_p,1)  % we need to shorten wave_u
    %%%%
    %wave_u_final = wave_u(cutlength_lower:nu-cutlength_upper,:);
    %the above line is wrong - we lined them up and then by cutting
    % from the start they are no longer aligned.
    %Instead we must cut from the end (see below).
    %%%%
    wave_u_final = wave_u(1:nu-ndiff,:);
    wave_p_final = wave_p;

    %Now add the time information
    timestep = event_p(2,1)-event_p(1,1);
    timevector_u = 0:timestep:(size(wave_u_final,1)-1)*timestep;
    timevector_p = 0:timestep:(size(wave_p_final,1)-1)*timestep;
    wave_u_final(:,1) = timevector_u';
    wave_p_final(:,1) = timevector_p';

elseif size(wave_u,1)<size(wave_p,1)  % we need to shorten wave_p
    %wave_p_final = wave_p(cutlength_lower:np-cutlength_upper,:);
    wave_p_final = wave_p(1:np-ndiff,:);
    wave_u_final = wave_u;

    %Now add the time information
    timestep = event_p(2,1)-event_p(1,1);
    timevector_u = 0:timestep:(size(wave_u_final,1)-1)*timestep;
    timevector_p = 0:timestep:(size(wave_p_final,1)-1)*timestep;
    wave_u_final(:,1) = timevector_u';
    wave_p_final(:,1) = timevector_p';

else % neither needs to be changed
    wave_u_final = wave_u;
    wave_p_final = wave_p;
    wave_u_final(:,1) = event_u(:,1);
    wave_p_final(:,1) = event_p(:,1);
end

%% CHECK LATER
if plot_par.rem_noise == 1 % remove noise if requested
    % Now we must compute the noise to signal ratio of both waveforms.
    % Lets do the uperturbed wave first
    [ind_un_start,closest,mind] = find_closest(wave_u_final(:,1),plot_par.lbnd4nenergy,'euclidean');
    [ind_un_end,closest,mind] = find_closest(wave_u_final(:,1),plot_par.ubnd4nenergy,'euclidean');
    %[ind_us_start,closest,mind] = find_closest(wave_u_final(:,1),30,'euclidean');
    %[ind_us_end,closest,mind] = find_closest(wave_u_final(:,1),100,'euclidean');
    nenergy_u = compute_energy(wave_u_final(ind_un_start:ind_un_end,2)-median(wave_u_final(ind_un_start:ind_un_end,2)),wave_u_final(ind_un_start:ind_un_end,1))
    %senergy_u = compute_energy(wave_u_final(ind_us_start:ind_us_end,2),wave_u_final(ind_us_start:ind_us_end,1))
    %ns2ratio_u = nenergy_u/senergy_u;

    % Now we need to do the perturbed wave
    [ind_pn_start,closest,mind] = find_closest(wave_p_final(:,1),plot_par.lbnd4nenergy,'euclidean');
    [ind_pn_end,closest,mind] = find_closest(wave_p_final(:,1),plot_par.ubnd4nenergy,'euclidean');
    %[ind_ps_start,closest,mind] = find_closest(wave_p_final(:,1),30,'euclidean');
    %[ind_ps_end,closest,mind] = find_closest(wave_p_final(:,1),100,'euclidean');
    nenergy_p = compute_energy(wave_p_final(ind_pn_start:ind_pn_end,2)-median(wave_p_final(ind_pn_start:ind_pn_end,2)),wave_p_final(ind_pn_start:ind_pn_end,1))
    %senergy_p = compute_energy(wave_p_final(ind_ps_start:ind_ps_end,2),wave_p_final(ind_ps_start:ind_ps_end,1))
    %ns2ratio_p = nenergy_p/senergy_p;
end

% Now we do 1 final cut so that the waveforms start near the first arrival
[ind_final_cut,closest,mind] = find_closest(wave_p_final(:,1),plot_par.parrival,'euclidean');
wave_p_final2= wave_p_final(ind_final_cut:end,:);
wave_u_final2= wave_u_final(ind_final_cut:end,:);
% and foix up the time
wave_u_final2(:,1) = wave_u_final2(:,1)-wave_u_final2(1,1);
wave_p_final2(:,1) = wave_p_final2(:,1)-wave_p_final2(1,1);


[norm_maxcorr,sep, tmp_sep, tmp_xcorr, omega2] = seis_xcorr_tmp(wave_u_final2(:,2), wave_p_final2(:,2), wave_u_final2(:,1), wave_p_final2(:,1),plot_par,'i',diag_switch);

    
if diag_switch ==1
    f1 = figure
    hr1 = subplot(4,2,1:2);
    h1 = plot(event_u(:,1)+trans1, event_u(:,2),'b');
    hold on;
    h2 = plot(event_p(:,1)+trans2, event_p(:,2),'r');
    htopYlims = get(hr1,'Ylim');
    arriv_box = LOC_makebox(htopYlims,start_arrivwin,end_arrivwin);
    plot(arriv_box(:,1),arriv_box(:,2),'k');
    coda_box = LOC_makebox(htopYlims,start_codawin,end_codawin);
    plot(coda_box(:,1),coda_box(:,2),'k');
    legend([h1,h2],{event1name,event2name});
    Xlims = get(hr1,'Xlim')
    twindowline = LOC_make_twindowline(htopYlims,Xlims,twindow);
    plot(twindowline(:,1),twindowline(:,2))

    hr2l = subplot(4,2,3);
    [indstart1,closest,mind] = find_closest(event_u(:,1),start_arrivwin-trans1,'euclidean');
    [indend1,closest,mind] = find_closest(event_u(:,1),end_arrivwin-trans1,'euclidean');
    h1 = plot(event_u(indstart1:indend1,1)+trans1, event_u(indstart1:indend1,2),'b');
    hold on;
    [indstart2,closest,mind] = find_closest(event_p(:,1),start_arrivwin-trans2,'euclidean');
    [indend2,closest,mind] = find_closest(event_p(:,1),end_arrivwin-trans2,'euclidean');
    h2 = plot(event_p(indstart2:indend2,1)+trans2, event_p(indstart2:indend2,2),'r');
    set(hr2l,'XTickLabel',{''},'YTickLabel',{''},'xlim',[max([event_u(indstart1),event_p(indstart2)]) ...
        min([event_u(indend1),event_p(indend2)])]);
    title('Left window - Direct arrival');

    hr2r=subplot(4,2,4);
    [indstart1,closest,mind] = find_closest(event_u(:,1),start_codawin-trans1,'euclidean');
    [indend1,closest,mind] = find_closest(event_u(:,1),end_codawin-trans1,'euclidean');
    h1 = plot(event_u(indstart1:indend1,1)+trans1, event_u(indstart1:indend1,2),'b');
    hold on;
    [indstart2,closest,mind] = find_closest(event_p(:,1),start_codawin-trans2,'euclidean');
    [indend2,closest,mind] = find_closest(event_p(:,1),end_codawin-trans2,'euclidean');
    h2 = plot(event_p(indstart2:indend2,1)+trans2, event_p(indstart2:indend2,2),'r');
    set(hr2r,'xTickLabel',{''},'YTickLabel',{''},'xlim',[max([event_u(indstart1),event_p(indstart2)]) ...
        min([event_u(indend1),event_p(indend2)])]);
    title('Right window - Section of coda');  
 


    hr3 = subplot(4,2,5:6)
    plot(norm_maxcorr(:,1),norm_maxcorr(:,2),'g')
    ylabel('X-xorr max')

    hr4 = subplot(4,2,7:8)
    plot(sep(:,1),sep(:,2),'g')
    hold on
    plot(sep(:,1), sep(:,3), 'r')
    plot(sep(:,1), sep(:,4), 'b')
    if src_type == 1| src_type==2
        ylabel('Separation (m)')
    elseif src_type ==3
        ylabel('r_t')
    end
    xlabel('time (s)')

end

% ind1 = find_closest(event1(:,1),trans1,'euclidean');
% ind2 = find_closest(event2(:,1),trans2,'euclidean');
% [maxcorr,sep,tmp_sep, tmp_xcorr ] = seis_xcorr(event1(ind1:end,2), event2(ind2:end,2), event1(ind1:end,1)-trans1, event2(ind2:end,1)-trans2,twindow,'ni');
% 
% 
% hr3 = subplot(4,2,5:6);
% %plot(maxcorr(:,1),maxcorr(:,2)/max(maxcorr(:,2)));
% hold on
% plot(maxcorr(:,1),maxcorr(:,2));
% 
% 
% hr3 = subplot(4,2,7:8);
% %plot(maxcorr(:,1),maxcorr(:,2)/max(maxcorr(:,2)));
% plot(sep(:,1),sep(:,2));
% xlims = get(hr3,'Xlim')
% actual_sep = sqrt(2*40^2);
% hold on
% plot(xlims, actual_sep*ones(1,2),'--')
% 
% 
% 
function box = LOC_makebox(Ylims,startwin,endwin)
% box   is a two coulum vector time coordinate in first column and amplitude
%       coordinate in the second column
sep = (Ylims(2)-Ylims(1))*5/100;  % add 5% separation to keep away from edge
box(1,:) = [startwin,   Ylims(1)+sep];  % bottom left
box(2,:) = [startwin,   Ylims(2)-sep];  % top left
box(3,:) = [endwin,     Ylims(2)-sep];  % top right
box(4,:) = [endwin,     Ylims(1)+sep];  % bottom right
box(5,:) = [startwin,   Ylims(1)+sep];  % bottom left  (close box)

function twindowline = LOC_make_twindowline(Ylims,Xlims,twindow)
sep = (Ylims(2)-Ylims(1))*10/100;  % add 10% separation to keep away from edge
xstart = 2*(Xlims(2)-Xlims(1))/3;

twindowline(1,:) = [xstart, Ylims(2)-sep];
twindowline(2,:) = [xstart+twindow, Ylims(2)-sep];


% 
%                             