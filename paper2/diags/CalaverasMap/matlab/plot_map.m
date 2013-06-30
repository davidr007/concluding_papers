% The purpose of this script is to plot a map of the events and stations

%====================================
% Define desired figures
figure_switch = 1;
%   1 => compute small figure with subset of stations only (FAST)
%   2 => also compute large figure with all stations (SLOW)
%====================================

fsize = 10;


% Load the earthquakes
load hypoDD.reloc
figure
plot(hypoDD(:,3), hypoDD(:,2), 'ko','markersize',2,'markerfacecolor','k')
hold on

s = importdata('station_subset.dat');
[n m] =size(s.data)
for i = 1:n
    plot(s.data(i,2), s.data(i,1), 'k^','markersize',7,'markerfacecolor','k')
    tmp =s.textdata(i);
    rn = 11-i;  %remove number
    i
    if i>1
        text(s.data(i,2)-0.025, s.data(i,1)+0.015,[tmp{1}(3:end),'(',num2str(rn),')'],'fontsize',fsize)
    elseif i==1
        text(s.data(i,2)-0.025, s.data(i,1)+0.015,[tmp{1}(3:end)],'fontsize',fsize)
    end
    disp([num2str(i),':',tmp{1}(3:end)])
end
%set(gca,'ylim',[37 37.5],'xlim',[-121.9 -121.5])
set(gca,'xlim',[-121.9,-121.44],'ylim',[37,37.4])
set(gca,'ytick', [37:0.1:37.4],'xtick',[-121.9:0.1:-121.5])

ylabel('Latitude','fontsize',fsize)
xlabel('Longitude','fontsize',fsize)
set(gca,'fontsize',fsize)
%axis equal

% Let's put on a scale bar
[s, alpha12, alpha21] = vincenty_inverse(37.35,-121.65,37.35,-121.53705)
plot([-121.65 -121.53705],[37.35 37.35],'k','linewidth',2)
text(mean([-121.65 -121.53705])-0.015,37.36,'10km','fontsize',fsize)
arrow([-121.6169,37.3048],[-121.6563,37.2916],10)
text(-121.61,37.3180,'Calaveras')
text(-121.61,37.2964,'Events')


statwidth = 7.5;
statheight = 7.5;
statxstart = 1.5%4.6;
statystart = 1.5%1.95;
set(gca,'units','centimeters')
set(gca, 'position',[statxstart, statystart, statwidth, statheight])

%print -depsc Calaveras_substationmap.eps
print -depsc ../../../Figure7.eps

if figure_switch ==2
    figure
    plot(hypoDD(:,3), hypoDD(:,2), 'ko','markersize',2,'markerfacecolor','k')
    hold on

    s = importdata('station.dat');
    [n m] =size(s.data)
    for i = 1:n
        plot(s.data(i,2), s.data(i,1), 'r^','markersize',7,'markerfacecolor','r')
        tmp =s.textdata(i);
        text(s.data(i,2)-0.01, s.data(i,1)+0.01,tmp{1}(3:end))
    end
    %set(gca,'ylim',[37 37.5],'xlim',[-121.9 -121.5])
    axis equal
end
