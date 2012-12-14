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
l1 = plot(t1,N1,'linewidth',2,'color',[0.7,0.7,0.7])

t2 = [tstart:tend];
N2 = Cumm_Omori(t2);
hold on
l2 = plot(t2,N2,'k','linewidth',2)

tend2 = 1000;
t3 = [tend:tend2];
N3 = Cumm_Omori(t3);
l3 = plot(t3,N3,'linewidth',2,'color',[0.7,0.7,0.7])

fsize = 16;
%set(gca,'fontsize',fsize)
xlabel('t (days)')%,'fontsize',fsize)
%ylabel('Cummulative number of aftershocks','fontsize',fsize)

ylabel('N(t)')
set(gca,'xlim',[-25, tend2],'xtick',[0:200:tend2])

set(gca,'units','centimeters','position',[4.6,1.5,7.5,6])
%print -depsc OmoriFigure.eps
print -depsc ../Figure12_bw.eps

% Now create the colour version
set(l1,'color','r')
set(l2,'color','b')
set(l3,'color','r')
print -depsc ../Figure12_c.eps


function [N] = Cumm_Omori(t)
% Assume the following values
K = 906.5;
c = 1.433; 
p = 1.256;

N = K *( c^(1-p) - (t+c).^(1-p))/(p-1);
