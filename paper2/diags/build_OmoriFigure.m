function [] = build_OmoriFigure()

% The purpose of this function is to build the Omori figure...

tstart = 4;
tdays = 150;
tend = tstart+tdays; 
ttotal = 1000;

Nstart =  Cumm_Omori(tstart)
Nend = Cumm_Omori(tend)
Ntotal = Cumm_Omori(ttotal)

prop_recorded = (Nend-Nstart)/Ntotal

t1 = [0:tstart];
N1 = Cumm_Omori(t1);
plot(t1,N1,'k--','linewidth',4)

t2 = [tstart:tend];
N2 = Cumm_Omori(t2);
hold on
plot(t2,N2,'k','linewidth',4)

t3 = [tend:200];
N3 = Cumm_Omori(t3);
plot(t3,N3,'k--','linewidth',4)

fsize = 16;
set(gca,'fontsize',fsize)
xlabel('time (days)','fontsize',fsize)
ylabel('Cummulative number of aftershocks','fontsize',fsize)

print -depsc OmoriFigure.eps

function [N] = Cumm_Omori(t)
% Assume the following values
K = 906.5;
c = 1.433; 
p = 1.256;

N = K *( c^(1-p) - (t+c).^(1-p))/(p-1);
