datadir = 'C:\datafiles\workstuff\sandpit\davidr\thesis_version2\diags\eq_location_optimisation\example_2D_50eq_3\';

if strcmp(filesep,'/')
    error('You must be on Windows to get the Latex Interpreter for these figures')
end

opttype = 'fortran'; % 'fortran' or 'matlab'
load([datadir,'randsearch_',opttype,'.mat'])
[sorted_x2, ind]=sort(xend(2,:));

[nevnts, nrand] = size(xend);

% Have a quick look at the solution space
figure
colorlookup = {'b.','r.','g.','c.','m.','y.','k.','bo','ro','go'}
colorlookup2 = {'b^','r^','g^','c^','m^','y^','k^','bd','rd','gd'}
for i = 1:nevnts
    h1 = plot(xend(i,:), yend(i,:),'ko','markerfacecolor','k', 'markersize',10);
    hold on
    h2 = plot(E(i,1),E(i,2), 'k^', 'markerfacecolor','k','markersize',10);
    hold on
    plot([xend(i,:) E(i,1)], [yend(i,:) E(i,2)],'Color',[0.5 0.5 0.5],'linewidth',2)
end
plot([ 0 0], get(gca,'ylim'),'k')
plot(get(gca,'xlim'),[0 0],'k')
legend([h2,h1],{'Actual Location','Optimisation Soln'},'Location','NorthWest') 



figure
colorlookup = {'b+','bo','b^','bV','bd','bs','b*','bp','b>','b<'}
colorlookup2 = {'r+','ro','r^','rV','rd','rs','r*','rp','r>','r<'}
ind = find(f<-127.38);
for i = 1:nevnts
    h1 = plot(E(i,1), E(i,2), 'ko', 'markerfacecolor','w','markersize',10);
    hold on
    %plot(xend(i,ind(j)),-yend(i,ind(j)), colorlookup{i})
    for j = 1:length(ind)
        h2 = plot(xend(i,ind(j)),yend(i,ind(j)), 'k^','markerfacecolor','w','markersize',10)
        hold on
        plot([xend(i,ind(j)) E(i,1)], [yend(i,ind(j)) E(i,2)],'Color',[0.7 0.7 0.7],'linewidth',2,'markersize',10)
    end
end
set(gca,'xlim',[-70,70], 'ylim',[-70,70], 'fontsize',18)
xlabel('$\hat{x}$ (m)','Interpreter','Latex')
ylabel('$\hat{y}$ (m)','Interpreter','Latex') 
legend([h2,h1],{'Actual Location','Optimisation Soln'},'Location','NorthWest') 

print -depsc locs_2D_50eq_3.eps

[nevnts, nrand] = size(xend);
% compute error as distancenorm
for k = 1:nrand
    count = 1
   sq_seperror(k) = sum( sqrt( (E(:,1) - xend(:,k)).^2 + (E(:,2) - yend(:,k)).^2 ) );
end
[sortedf,indf] = sort(f);
figure
subplot(2,1,1)
plot(sortedf)
ylabel('Lstar')
xlabel('Run Number after sorting')
subplot(2,1,2)
plot(sq_seperror(indf)/(nevnts*2-2-1))
ylabel('sep error per free location parameter')
xlabel('Run Number after sorting')
disp(['sep per invertable parameter = ', num2str(sq_seperror(indf)/(nevnts*2-2-1))])

% 




