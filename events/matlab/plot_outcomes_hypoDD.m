function [hypoDDstats,Locs,fh] = plot_outcomes_hypoDD(dirname,plot_struct,colorchoice,path2thesis)
% This function plots the hypoDD outcomes. It is typically called by 
% batch_doplots. 
%
% INPUTS: 
% dirname       [char]  director where hypoDD outputs are stored. Note that 
%               eps figures are saved here as well. 
% plot_struct   [structure] settings for the figures. 
% colorchoice 
% path2thesis   [string] path to where the thesis is stored on the macjhine
%               being used
%
% OUTPUTS: 
% hypoDDstats   [structure] contains statistics of interest
%                   hypoDDstats.minval = min(abs(delta_coord));
%                   hypoDDstats.maxval = max(abs(delta_coord));
%                   hypoDDstats.meanval = mean(abs(delta_coord));
%                   hypoDDstats.stdval = std(abs(delta_coord));
%                   hypoDDstats.nE = nE; 
%                   hypoDDstats.plane => info on best fitting plane;    
% Locs          [structure] contains the locations in desired coordinate system
%                   Locs.relocx = hypoDD_reloc(i,5)
%                   Locs.relocy = hypoDD_reloc(i,6)
%                   Locs.relocz = hypoDD_reloc(i,7)
%                   Locs.not68x = hypodd_not68(i,5)
%                   Locs.not68y = hypodd_not68(i,6)
%                   Locs.not68z = hypodd_not68(i,7)
% fh            double [1Xn] figure handles
%
% David Robinson
% 9 December 2009

if nargin ==2 | isempty(colorchoice)
    colorchoice = 'c';
end

switch colorchoice
    case {'c'}
        c1 = 'b';
        c2 = 'r';
    case {'bw'}
        c1 = [0.5,0.5,0.5];
        c2 = 'k';
    otherwise
        error('ERROR: invalid value for colorchoice')
end

if isempty(plot_struct)
    plot_struct.xlimits = [-500 500];
    plot_struct.ylimits = [-500 500];
    plot_struct.zlimits = [-2000 2000];
    plot_struct.axes_posy = [5,3,6,6];
    plot_struct.xtickspots = [-250 0 250];
    plot_struct.ytickspots = [-250 0 250];
    plot_struct.ztickspots = [-1000 0 1000];
    plot_struct.fsize = 16;
    plot_struct.msize = 4;
end



% unwrap the plot structure
xlimits = plot_struct.xlimits;
ylimits = plot_struct.ylimits;
zlimits = plot_struct.zlimits;
axes_posy = plot_struct.axes_posy;
xtickspots = plot_struct.xtickspots;
ytickspots = plot_struct.ytickspots;
ztickspots = plot_struct.ztickspots;
fsize = plot_struct.fsize;
msize = plot_struct.msize;

if strcmp(filesep,'/') % i.e. in Linus
    parentdir = [path2thesis,'/diags/eq_location_optimisation/'];
else 
    parentdir = [path2thesis,'/diags/eq_location_optimisation/'];
%    error('Not yet operation for the Windows machine')
end

%%%%% Now let's put them onto the calaveras streaks
hypodd = load([parentdir,'hypoDD.reloc']);  % loads the complete hypoDD relocations for the 308 events
%hypodd_68 = load('HYPODD/runs/orig_68eqonly/hypoDD.reloc');  % loads the best relocations for the 68 events
hypodd_68 = load([parentdir,'hypoDD_loc68.txt']);



[tmp,ind101] = setdiff(hypodd(:,1),hypodd_68(:,1));
hypodd_not68 = hypodd(ind101,:);

hypoDD_reloc = load([dirname,filesep,'hypoDD.reloc']);
[nE,m] = size(hypoDD_reloc);

% Now we must translate the latest hypoDD outputs to be in the correct
% place for plotting

ind = find(hypodd(:,1) == hypodd_68(1,1));  %Note that hypodd_68(1,1)= 17842
translation = hypoDD_reloc(1,5:7) - hypodd(ind,5:7);
hypoDD_reloc(:,5:7) = hypoDD_reloc(:,5:7) - repmat(translation,nE,1);
hypoDDstats.translation = translation;
if isempty(ind)
    error('I am bombing')
end
disp(['hypodd_68(1,1) = ', num2str(hypodd_68(1,1))])
disp(['hypodd(ind,1) = ', num2str(hypodd(ind,1))])
fh(1) = figure
plot(hypodd_not68(:,5), hypodd_not68(:,6),'o','markeredgecolor',c1,'markerfacecolor',c1,'markersize',msize)
hold on
% for i = 1:nE
%     plot(hypoDD_reloc(i,5), hypoDD_reloc(i,6),'ro','markeredgecolor',c2,'markerfacecolor',c2,'markersize',msize)
% end
plot(hypoDD_reloc(:,5), hypoDD_reloc(:,6),'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',msize)

%     plot(xend(1,ind), yend(1,ind),'go','markerfacecolor','g','markersize',msize)
%     plot(xend(2,ind), yend(2,ind),'yo','markerfacecolor','y','markersize',msize)
%     plot(xend(3,ind), yend(3,ind),'co','markerfacecolor','c','markersize',msize)
    

xlabel('x (m)','fontsize',fsize)
ylabel('y (m)','fontsize',fsize)
set(gca,'xlim',xlimits, 'ylim',ylimits,'units','centimeters')
set(gca,'position',axes_posy)
set(gca,'xtick', xtickspots, 'ytick',ytickspots)
set(gca,'fontsize',fsize)
print([dirname,filesep,'example_68eq_hypoDD_xy.eps'], '-depsc') 

% First plot the catalogue locations
fh(2) = figure
plot(hypodd_not68(:,5), hypodd_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c1,'markersize',msize)
hold on
% for i = 1:nE
%     plot(hypoDD_reloc(i,5), hypoDD_reloc(i,7),'ro','markeredgecolor',c2,'markerfacecolor',c2,'markersize',msize)
% end
plot(hypoDD_reloc(:,5), hypoDD_reloc(:,7),'ro','markeredgecolor',c2,'markerfacecolor',c2,'markersize',msize)

%     plot(xend(1,ind), zend(1,ind),'go','markerfacecolor','g','markersize',msize)
%     plot(xend(2,ind), zend(2,ind),'yo','markerfacecolor','y','markersize',msize)
%     plot(xend(3,ind), zend(3,ind),'co','markerfacecolor','c','markersize',msize)

xlabel('x (m)','fontsize',fsize)
ylabel('z (m)','fontsize',fsize)
set(gca,'xlim',xlimits, 'ylim',zlimits,'units','centimeters')
set(gca,'position',axes_posy)
set(gca,'xtick', xtickspots, 'ytick',ztickspots)
set(gca,'fontsize',fsize)
print([dirname,filesep,'example_68eq_hypoDD_xz.eps'], '-depsc') 

% First plot the catalogue locations
fh(3) = figure
plot(hypodd_not68(:,6), hypodd_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c1,'markersize',msize)
hold on
% for i = 1:nE
%     plot(hypoDD_reloc(i,6), hypoDD_reloc(i,7),'ro','markeredgecolor',c2,'markerfacecolor',c2,'markersize',msize)
% end
plot(hypoDD_reloc(:,6), hypoDD_reloc(:,7),'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',msize)

%     plot(yend(1,ind), zend(1,ind),'go','markerfacecolor','g','markersize',msize)
%     plot(yend(2,ind), zend(2,ind),'yo','markerfacecolor','y','markersize',msize)
%     plot(yend(3,ind), zend(3,ind),'co','markerfacecolor','c','markersize',msize)

xlabel('y (m)','fontsize',fsize)
ylabel('z (m)','fontsize',fsize)
set(gca,'xlim',xlimits, 'ylim',zlimits,'units','centimeters')
set(gca,'position',axes_posy)
set(gca,'xtick', ytickspots, 'ytick',ztickspots)
set(gca,'fontsize',fsize)
print([dirname,filesep,'example_68eq_hypoDD_yz.eps'], '-depsc') 


% Do some simple statistics
E = zeros(nE,3);
for i = 1:nE
   ind = find(hypodd_68(:,1) ==hypoDD_reloc(i,1));
   E(i,:) = hypodd_68(ind,5:7);
end

hypoDD_origx = E(:,1);
best_hypoDDx = hypoDD_reloc(:,5);
hypoDD_origy = E(:,2);
best_hypoDDy = hypoDD_reloc(:,6);
hypoDD_origz = E(:,3);
best_hypoDDz = hypoDD_reloc(:,7);

deltax = hypoDD_origx-best_hypoDDx;
deltay = hypoDD_origy-best_hypoDDy;
deltaz = hypoDD_origz-best_hypoDDz;
delta_coord = [deltax;deltay;deltaz];

hypoDDstats.minval = min(abs(delta_coord));
hypoDDstats.maxval = max(abs(delta_coord));
hypoDDstats.meanval = mean(abs(delta_coord));
hypoDDstats.stdval = std(abs(delta_coord));
hypoDDstats.nE = nE; 


% Now let us find the best fitting plane
points_in = [best_hypoDDx,best_hypoDDy,best_hypoDDz];
[x0, normal, resid, norm_resid] = lsplane(points_in);
hypoDDstats.plane.normal = normal;
hypoDDstats.plane.resid = resid;
hypoDDstats.plane.norm_resid = norm_resid;
hypoDDstats.plane.rms = norm_resid/sqrt(length(resid));

Locs.relocx = hypoDD_reloc(:,5); 
Locs.relocy = hypoDD_reloc(:,6);
Locs.relocz = hypoDD_reloc(:,7);
Locs.sigmax = hypoDD_reloc(:,8); 
Locs.sigmay = hypoDD_reloc(:,9);
Locs.sigmaz = hypoDD_reloc(:,10);
Locs.not68x = hypodd_not68(:,5);
Locs.not68y = hypodd_not68(:,6);
Locs.not68z = hypodd_not68(:,7);

