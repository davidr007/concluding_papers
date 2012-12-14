% This script wraps the plotting function plot_outcomes.m 
% assembles all of the data in a convenient format and then 
% draws some plots analysing the simulations


%if strcmp(filesep,'/')
%    error('ERROR: you must run this on the Windows notebook to get the Latex figure labels')
%end

dirlist = { 'const_010perc'; 'const_020perc'; 'const_030perc'; ...
            'const_040perc'; 'const_050perc'; ...
            'const_060perc'; 'const_070perc'; 'const_080perc'; ...
            'const_090perc'; 'const_100perc'};

maxvals_all = zeros(1,length(dirlist));        
meanvals_all = zeros(1,length(dirlist));
stdvals_all = zeros(1,length(dirlist));
minvals_all = zeros(1,length(dirlist));
maxvals_best = zeros(1,length(dirlist));        
meanvals_best = zeros(1,length(dirlist));
stdvals_best = zeros(1,length(dirlist));
minvals_best = zeros(1,length(dirlist));
for i = 1:length(dirlist)
    [res.(dirlist{i}),fh] = plot_outcomes(dirlist{i});
    maxvals_all(i) = res.(dirlist{i}).allsol.maxval;
    meanvals_all(i) = res.(dirlist{i}).allsol.meanval;
    stdvals_all(i) = res.(dirlist{i}).allsol.stdval;
    minvals_all(i) = res.(dirlist{i}).allsol.minval;
    maxvals_best(i) = res.(dirlist{i}).bestsol.maxval;
    meanvals_best(i) = res.(dirlist{i}).bestsol.meanval;
    stdvals_best(i) = res.(dirlist{i}).bestsol.stdval;
    minvals_best(i) = res.(dirlist{i}).bestsol.minval;
    ticklabelstr{i} = [num2str(str2num(dirlist{i}(7:9))),'%'];
    figure(fh.allsolutions)
    switchval = dirlist{i};
    switch switchval
        case {'const_040perc'; 'const_050perc'; 'const_060perc'; ...
              'const_070perc'; 'const_080perc'; 'const_090perc'; 'const_100perc'}
            set(gca,'xlim',[-70,70],'ylim',[-70,70])
            set(gca,'xtick',[-50,0 50],'ytick',[-50,0 50])
        case {'const_030perc'; 'const_020perc'}
            set(gca,'xlim',[-110,110],'ylim',[-110,110])
            set(gca,'xtick',[-80,0 80],'ytick',[-80,0 80])
        case {'const_010perc'}
            set(gca,'xlim',[-150,150],'ylim',[-150,150])
            set(gca,'xtick',[-80,0 80],'ytick',[-80,0 80])
    end
    set(gca,'xminortick','on','yminortick','on')
    plot([ 0 0], get(gca,'ylim'),'k')
    plot(get(gca,'xlim'),[0 0],'k')
    fsize = 22;
    set(gca,'fontsize',fsize)
    title(ticklabelstr{i},'fontsize',fsize)
    xlabel('$\hat{x}$\,(m)','Interpreter','LaTex','fontsize',fsize)
    ylabel('$\hat{y}$\,(m)','Interpreter','LaTex','fontsize',fsize)
    set(gcf,'units','centimeters')
    set(gca,'units','centimeters')
    set(gca,'position',[ 2.7760, 2.7, 8.6, 8.6])
    print('-depsc',['soln_space_2Dsynth50eq_',dirlist{i},'.eps'])
    figure(fh.bestsolution)
      switch switchval
        case {'const_040perc'; 'const_050perc'; 'const_060perc'; ...
              'const_070perc'; 'const_080perc'; 'const_090perc'; 'const_100perc'}
            set(gca,'xlim',[-70,70],'ylim',[-70,70])
            set(gca,'xtick',[-50,0 50],'ytick',[-50,0 50])
        case {'const_030perc'; 'const_020perc'}
            set(gca,'xlim',[-150,150],'ylim',[-150,150])
            set(gca,'xtick',[-100,0 100],'ytick',[-100,0 100])
        case {'const_010perc'}
            set(gca,'xlim',[-150,150],'ylim',[-150,150])
            set(gca,'xtick',[-100,0 100],'ytick',[-100,0 100])
    end
    set(gca,'xminortick','on','yminortick','on')
    plot([ 0 0], get(gca,'ylim'),'k')
    plot(get(gca,'xlim'),[0 0],'k')
    fsize = 22;
    set(gca,'fontsize',fsize)
    title(ticklabelstr{i},'fontsize',fsize)
    xlabel('$\hat{x}$\,(m)','Interpreter','LaTex','fontsize',fsize)
    ylabel('$\hat{y}$\,(m)','Interpreter','LaTex','fontsize',fsize)
    set(gcf,'units','centimeters')
    set(gca,'units','centimeters')
    set(gca,'position',[ 2.7760, 2.7, 8.6, 8.6])
    print('-depsc',['soln_space_2Dsynth50eq_bestsoln',dirlist{i},'.eps'])
end


fsize = 12
% Now let's do the plotting
fh.maxvals = figure;
sh(1) = subplot(3,1,1)
h1b = plot(maxvals_best,'color','k','linewidth',2) %  best
hold on
h1g = plot(maxvals_all,'color',[0.7,0.7,0.7],'linewidth',1); % all
%legend([h1,h2], {'all', 'best'})
set(gca,'xtick',1:10,'ytick',[0 50 100 150])
set(gca,'xticklabel',{}) %ticklabelstr)
set(gca,'fontsize',fsize)
set(gca,'ylim',[0 150],'xlim',[1 10],'units','centimeters')
lh = legend([h1b,h1g],{'best','all'},'Location','Northeast')
set(lh,'fontsize',10)
set(lh,'Box','off')
ylabel('$\Delta_{max}$\,(m)','fontsize',fsize,'Interpreter','LaTex')
%xlabel('Number of constraints','fontsize',fsize)

sh(2) = subplot(3,1,2)
h2b = plot(meanvals_best,'color','k','linewidth',2) % best
hold on
h2g = plot(meanvals_all,'color',[0.7,0.7,0.7],'linewidth',1); % all
hold on
set(gca,'xtick',1:10,'ytick',[0 20 40])
set(gca,'xticklabel',{})%ticklabelstr)
set(gca,'fontsize',fsize)
set(gca,'ylim',[0 40],'xlim',[1 10],'units','centimeters')
ylabel('$\Delta_\mu$\,(m)','fontsize',fsize,'Interpreter','LaTex')

% subplot(4,1,3)
% h1 = plot(stdvals_all,'b','linewidth',2);
% hold on
% h2 = plot(stdvals_best,'r','linewidth',1)
% set(gca,'xtick',1:10,'ytick',[0 20 40])
% set(gca,'xticklabel',ticklabelstr)
% set(gca,'fontsize',fsize)
% set(gca,'ylim',[0 40],'xlim',[1 10])
% ylabel('$\Delta_\sigma$\,(m)','fontsize',fsize,'Interpreter','LaTex')


% now let us draw the constraints 
% setup a colormap for the pathlengths

colormap(jet(8))
C = colormap;
C = [C;1 1 1]; % add a white line at the end (i.e. for unconnected paths

for i = 1:length(dirlist)
    figure
    load([dirlist{i},filesep,'CWI_stat.txt']);
    [nCWI mCWI] = size(CWI_stat);
    E = res.(dirlist{i}).otherstuff.E;
    count = 1;
    segments = [];
    for j = 1:nCWI
        if  CWI_stat(j,3) ~= -99999
           plot([E(CWI_stat(j,1),1),E(CWI_stat(j,2),1)], [E(CWI_stat(j,1),2),E(CWI_stat(j,2),2)],'Color',[0.5,0.5,0.5]) 
           hold on
           % build a segment list for use with dijkstra later
           segments = [segments; count,CWI_stat(j,1),CWI_stat(j,2)];
           count = count+1;
        end
    end
    plot(E(:,1),E(:,2), 'bo', 'markersize',7)
    set(gca,'xlim',[-60,60],'ylim',[-60,60])
    set(gca,'xtick',[-30,0 30],'ytick',[-30,0 30])
    set(gca,'xminortick','on','yminortick','on')
    plot([ 0 0], get(gca,'ylim'),'k')
    plot(get(gca,'xlim'),[0 0],'k')
    fsize = 22;
    set(gca,'fontsize',fsize)
    title(ticklabelstr{i},'fontsize',fsize)
    xlabel('$\hat{x}$\,(m)','Interpreter','LaTex','fontsize',fsize)
    ylabel('$\hat{y}$\,(m)','Interpreter','LaTex','fontsize',fsize)
    set(gcf,'units','centimeters')
    set(gca,'units','centimeters')
    set(gca,'position',[ 2.7760, 2.7, 8.6, 8.6])
    print('-depsc',['CWIconst_links_',dirlist{i},'.eps'])    
    
    
    % Now let's have a look at the node linkages.......
    % We will use the third-party dijkstra2 to do this
    nodes = zeros(50,3);
    nodes(:,1) = 1:50;
    nodes(:,2:3) = E;
    for j = 1:50 % loop over the events
        [dist,path] = dijkstra2(nodes,segments,j);
        linkage_stats.dist.(dirlist{i}).(['E',num2str(j)]) = dist;
        linkage_stats.path.(dirlist{i}).(['E',num2str(j)]) = path;
    end   
    
    pathlength = zeros(50,50);
    pathlength_vec_tmp = [];
    for j = 1:50
        for k = 1:50 
            tmpdist = linkage_stats.dist.(dirlist{i}).(['E',num2str(j)])(k);
            if isinf(tmpdist)
                pathlength(j,k) = Inf;
            else
                pathlength(j,k) = length(linkage_stats.path.(dirlist{i}).(['E',num2str(j)]){k})-1;       
            end
            if k >j
                pathlength_vec_tmp = [pathlength_vec_tmp, pathlength(j,k)];
            end
        end
    end
    linkage_stats.pathlength.(dirlist{i}) = pathlength;
    average_pathlength_vec(i) = mean(pathlength_vec_tmp(~isinf(pathlength_vec_tmp)));
    
    
    figure
    % re-assign repeated data to white;
    tmp_mat = Inf*ones(50,50);
    for j = 1:50
        for k = j:50
            tmp_mat(j,k) = linkage_stats.pathlength.(dirlist{i})(j,k);
        end
    end
    linkage_stats.pathlength.(dirlist{i}) = tmp_mat;     
    imagesc(linkage_stats.pathlength.(dirlist{i}),[-0.5 8.5])
    title(ticklabelstr{i},'fontsize',fsize)
    set(gca,'fontsize',fsize)
    colormap(C)
    h = colorbar
    set(h,'fontsize',16)
    set(h,'Ytick',[0:8])
    set(h,'YTickLabel',{0,1,2,3,4,5,6,7,'NA'})
    set(gcf,'units','centimeters')
    set(gca,'units','centimeters')
    set(gca,'position',[ 2.7760, 2.7, 8.6, 8.6])
    print('-depsc',['CWIconst_linkmatrices_',dirlist{i},'.eps'])    
    
end

fsize = 12
figure(fh.maxvals)
sh(3) = subplot(3,1,3)
h3b = plot(average_pathlength_vec,'k','linewidth',2);
set(gca,'xtick',1:10)
%set(gca,'ytick',[0 50 100 150])
ticklabelstr = {' ','20%',' ','40%',' ','60%',' ','80%',' ','100%'};
set(gca,'xticklabel',ticklabelstr)
set(gca,'fontsize',fsize)
%set(gca,'ylim',[0 150],
set(gca,'xlim',[1 10],'units', 'centimeters')
ylabel('$\mu_{links}$','fontsize',fsize,'Interpreter','LaTex')
hold on
plot(get(gca,'xlim'),[2 2], 'k--')
xlabel('Number of constraints','fontsize',fsize)

statwidth = 7.5;
statheight = 1.5;
statxstart = 4.6;
statystart = 1.95;
statygap = 0.7;

set(sh(3), 'position',[statxstart, statystart, statwidth, statheight])
set(sh(2), 'position',[statxstart, statystart+statheight+statygap, statwidth, statheight])
set(sh(1), 'position',[statxstart, statystart+2*statheight+2*statygap, statwidth, statheight])

%xlabel('Number of constraints','fontsize',fsize)
%print -depsc ressummary_2Dsynth50eq.eps
print -depsc ../../Figure3_bw.eps


% Now lets make the colour version
set(h1b,'color','b')
set(h1g,'color','r')
set(h2b,'color','b')
set(h2g,'color','r')
print -depsc ../../Figure3_c.eps
