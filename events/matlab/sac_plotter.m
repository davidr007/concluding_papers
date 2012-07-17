function event = sac_plotter(fname,entypcurrent)

% small function to wrap Erdinc's Sac utilities to plot the trace with
% picks if they exist.
%
% David Robinson
% November 2006


%[data,kstnm,delta,b,e,kcmnp,nzyear,nzjday,nzhour,nzmin,...
%    nzsec,stla,stlo,evla,evlo,cmpaz,cmpinc,tp,ts,t1,t2] = sacread(fname,entypcurrent)
event = sacread(fname,entypcurrent)
plot(event.b+[0:event.delta:(length(event.data)-1)*event.delta],event.data,'k');
ylimits = get(gca,'ylim');
hold on
% plot precise tp pick
if event.tp ~= -12345
    plot([event.tp event.tp],ylimits,'r--')
end
% plot precise ts pick
if event.ts ~=-12345
    plot([event.ts event.ts],ylimits,'r--')
end
