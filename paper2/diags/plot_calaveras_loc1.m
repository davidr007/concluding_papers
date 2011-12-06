% The purpose of this script is to plot the Calaveras earthquake locations
% IMPORTANT NOTE
% .... This version (i.e. plot_calaveras_loc1.m utilises hypoDD with
% .... LSQR method........

%% Machine check - must be on notebook for best figure
if strcmp(filesep,'/')
    error('ERROR: must be on notebook for best figure quality')
end

%% Figure setup

colorswitch = 'bw';  %'c'=>color OR 'bw'=>black/white
switch colorswitch
    case 'bw'
        c1 = [0.7, 0.7, 0.7];   % background event - edge colour
        c2 = 'k';               % event of interest - edge and face color
        c3 = 'w';               % background event - face colour
    case 'c'
        c1 = 'b';               % background event - edge colour
        c2 = 'r';               % event of interest - edge and face color
        c3 = 'b';               % background event - face colour
end



plot_struct.xlimits = [-500 500];
plot_struct.ylimits = [-500 500];
plot_struct.zlimits = [-2000 2000];
plot_struct.axes_posy = [5,3,6,6];
plot_struct.xtickspots = [-250 0 250];
plot_struct.ytickspots = [-250 0 250];
plot_struct.ztickspots = [-1000 0 1000];
plot_struct.zticklabels =[-1000, 0, 1000];
plot_struct.fsize = 10;
plot_struct.msize = 2;

awidth = 3; % width of axis in centimeters
aheight = 3; % height of axis in centimeters
hgap = 0.2; % horizontal gap in centimeters
vgap = 0.2; % vertical gap in centimeters
ax_start = 1.7;
ay_start = 0.7;


%% Load inputdata
events_oi = load('../../../thesis_version2/diags/eq_location_optimisation/events_oi.txt');

CatLoc = load('../../../thesis_version2/diags/eq_location_optimisation/Calaveras_cat_locations.txt'); 
[n m] = size(CatLoc); 
CatLoc_68 = [];
CatLoc_not68 = [];
for i = 1:n
    ind = find(CatLoc(i,1)==events_oi);
    if ~isempty(ind) 
        CatLoc_68 = [CatLoc_68; CatLoc(ind,:)];
    else
        CatLoc_not68 = [CatLoc_not68; CatLoc(i,:)];
    end
end

%% Figure 1 - Catalogue - CWI and Travel Time Locations
% 3x3 figures showing:
%       col1 => Catalogue Locations
%       col2 => CWI locations
%       col3 => HypoDD locations

% Lets get the results we need for this figure
if strcmp(filesep,'/')
    dirname = '/export/storage/davidr/sandpit/davidr/thesis_version2/diags/eq_location_optimisation/';
else
    dirname = 'c:/datafiles/workstuff/sandpit/davidr/thesis_version2/diags/eq_location_optimisation/';
end
tmphypoDD68 = load([dirname,'hypoDD_loc68.txt']);
Locs_hypoDDbest.relocx = tmphypoDD68(:,5);
Locs_hypoDDbest.relocy = tmphypoDD68(:,6);
Locs_hypoDDbest.relocz = tmphypoDD68(:,7);
tmphypoDDnot68 = load([dirname,'hypodd_not68.txt']);
Locs_hypoDDbest.not68x = tmphypoDDnot68(:,5);
Locs_hypoDDbest.not68y = tmphypoDDnot68(:,6);
Locs_hypoDDbest.not68z = tmphypoDDnot68(:,7);
if strcmp(filesep,'/')
    dirname = '/export/storage/davidr/sandpit/davidr/thesis_version2/diags/eq_location_optimisation/example_69eq_calaveras3';
else
    dirname = 'c:/datafiles/workstuff/sandpit/davidr/thesis_version2/diags/eq_location_optimisation/example_69eq_calaveras3';
end
%dirname = '/export/storage/davidr/sandpit/davidr/thesis_version2/diags/eq_location_optimisation/CalaverasMultiSim/CWI/stat10';
[CWIstats_best,Locs_CWIbest,fh] = plot_outcomes_cwi(dirname,plot_struct);
for i = 1:length(fh)
    close(fh(i))
end


figure
set(gcf,'units','centimeter')

% catalogue x vs y
subplot(3,3,1)
plot(CatLoc_not68(:,5), CatLoc_not68(:,6),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(CatLoc_68(:,5), CatLoc_68(:,6),'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
%xlabel('x (m)','fontsize',plot_struct.fsize)
%ylabel('y (m)','fontsize',plot_struct.fsize)
ylabel('x vs y (m)','fontsize',plot_struct.fsize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.ylimits,'units','centimeters')
%set(gca,'position',axes_posy)
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ytickspots)
set(gca,'xticklabel',[])
set(gca,'fontsize',plot_struct.fsize)
title('Catalogue','fontsize',plot_struct.fsize)
set(gca,'position',[ax_start, ay_start+2*aheight+2*vgap, awidth, aheight])

% CWI x vs y
subplot(3,3,2)
plot(CatLoc_not68(:,5), CatLoc_not68(:,6),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_CWIbest.relocx, Locs_CWIbest.relocy,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.ylimits,'units','centimeters')
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ytickspots)
set(gca,'xticklabel',[],'yticklabel',[])
set(gca,'fontsize',plot_struct.fsize)
title('Coda Waves','fontsize',plot_struct.fsize)
set(gca,'position',[ax_start+awidth+hgap, ay_start+2*aheight+2*vgap, awidth, aheight])

% hypoDD x vs y
subplot(3,3,3)
plot(CatLoc_not68(:,5), CatLoc_not68(:,6),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_hypoDDbest.relocx, Locs_hypoDDbest.relocy,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.ylimits,'units','centimeters')
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ytickspots)
set(gca,'xticklabel',[],'yticklabel',[])
set(gca,'fontsize',plot_struct.fsize)
title('\DeltaTT','fontsize',plot_struct.fsize)
set(gca,'position',[ax_start+2*awidth+2*hgap, ay_start+2*aheight+2*vgap, awidth, aheight])

% catalogue x vs z
subplot(3,3,4)
plot(CatLoc_not68(:,5), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(CatLoc_68(:,5), CatLoc_68(:,7),'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
%xlabel('x (m)','fontsize',plot_struct.fsize)
%ylabel('z (m)','fontsize',plot_struct.fsize)
ylabel('x vs z (m)','fontsize',plot_struct.fsize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
%set(gca,'position',axes_posy)
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ztickspots,'yticklabel',plot_struct.zticklabels)
set(gca,'xticklabel',[])
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
set(gca,'position',[ax_start, ay_start+aheight+vgap, awidth, aheight])

% CWI x vs z
subplot(3,3,5)
plot(CatLoc_not68(:,5), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_CWIbest.relocx, Locs_CWIbest.relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ztickspots,'yticklabel',plot_struct.zticklabels)
set(gca,'xticklabel',[],'yticklabel',[])
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
set(gca,'position',[ax_start+awidth+hgap, ay_start+aheight+vgap, awidth, aheight])

% hypoDD x vs z
subplot(3,3,6)
plot(CatLoc_not68(:,5), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_hypoDDbest.relocx, Locs_hypoDDbest.relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ztickspots,'yticklabel',plot_struct.zticklabels)
set(gca,'xticklabel',[],'yticklabel',[])
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
set(gca,'position',[ax_start+2*awidth+2*hgap, ay_start+aheight+vgap, awidth, aheight])

% catalogue y vs z
subplot(3,3,7)
plot(CatLoc_not68(:,6), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(CatLoc_68(:,6), CatLoc_68(:,7),'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
%xlabel('y (m)','fontsize',plot_struct.fsize)
%ylabel('z (m)','fontsize',plot_struct.fsize)
ylabel('y vs z (m)','fontsize',plot_struct.fsize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
%set(gca,'position',axes_posy)
set(gca,'xtick', plot_struct.ytickspots, 'ytick',plot_struct.ztickspots, 'yticklabel',plot_struct.zticklabels)
%set(gca,'yticklabel',[])
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
set(gca,'position',[ax_start, ay_start, awidth, aheight])

% CWI y vs z
subplot(3,3,8)
plot(CatLoc_not68(:,6), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_CWIbest.relocy, Locs_CWIbest.relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
set(gca,'xtick', plot_struct.ytickspots, 'ytick',plot_struct.ztickspots, 'yticklabel',plot_struct.zticklabels)
set(gca,'yticklabel',[])
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
set(gca,'position',[ax_start+awidth+hgap, ay_start, awidth, aheight])

% hypoDD y vs z
subplot(3,3,9)
plot(CatLoc_not68(:,6), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_hypoDDbest.relocy, Locs_hypoDDbest.relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
set(gca,'xtick', plot_struct.ytickspots, 'ytick',plot_struct.ztickspots, 'yticklabel',plot_struct.zticklabels)
set(gca,'yticklabel',[])
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
set(gca,'position',[ax_start+2*awidth+2*hgap, ay_start, awidth, aheight])

print -depsc CalaverasLoc1.eps


%% Figure 2 - CWI Reduces Stations
figure
% 3x6 figure showing:
%       col1 => CWI with 7 stations
%       col2 => CWI with 5 stations
%       col3 => CWI with 4 stations
%       col4 => CWI with 3 stations
%       col5 => CWI with 2 stations
%       col6 => CWI with 1 stations

% Lets get the results we need for this figure
statlist = {'stat7','stat5','stat4','stat3','stat2','stat1'};
for j = 1: length(statlist)
    if strcmp(filesep,'/')
        dirname = ['/export/storage/davidr/sandpit/davidr/thesis_version2/diags/eq_location_optimisation/CalaverasMultiSim/CWI/',statlist{j},'/'];
    else
        dirname = ['c:/datafiles/workstuff/sandpit/davidr/thesis_version2/diags/eq_location_optimisation/CalaverasMultiSim/CWI/',statlist{j},'/'];
    end
    [tmpstats,Locs_CWI_rs.(statlist{j}),fh] = plot_outcomes_cwi(dirname,plot_struct);
    for i = 1:length(fh)
        close(fh(i))
    end
end

for i = 1:18
    subplot(3,6,i)
end
sh = zeros(1,18);

% CWI x vs y (7 stations)
sh(1) = subplot(3,6,1);
plot(CatLoc_not68(:,5), CatLoc_not68(:,6),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_CWI_rs.(statlist{1}).relocx, Locs_CWI_rs.(statlist{1}).relocy,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
%xlabel('x (m)','fontsize',plot_struct.fsize)
%ylabel('y (m)','fontsize',plot_struct.fsize)
ylabel('x vs y (m)','fontsize',plot_struct.fsize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.ylimits,'units','centimeters')
%set(gca,'position',axes_posy)
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ytickspots)
set(gca,'fontsize',plot_struct.fsize)
title('7 Stations','fontsize',plot_struct.fsize)
set(gca,'xticklabel',[])

%set(gca,'position',[ax_start, ay_start+2*aheight+2*vgap, awidth, aheight])

% CWI x vs y (5 stations)
sh(2) = subplot(3,6,2);
plot(CatLoc_not68(:,5), CatLoc_not68(:,6),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_CWI_rs.(statlist{2}).relocx, Locs_CWI_rs.(statlist{2}).relocy,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.ylimits,'units','centimeters')
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ytickspots)
set(gca,'fontsize',plot_struct.fsize)
title('5 Stations','fontsize',plot_struct.fsize)
%set(gca,'position',[ax_start+awidth+hgap, ay_start+2*aheight+2*vgap, awidth, aheight])
set(gca,'xticklabel',[],'yticklabel',[])

% CWI x vs y (4 stations)
sh(3) = subplot(3,6,3);
plot(CatLoc_not68(:,5), CatLoc_not68(:,6),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_CWI_rs.(statlist{3}).relocx, Locs_CWI_rs.(statlist{3}).relocy,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.ylimits,'units','centimeters')
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ytickspots)
set(gca,'fontsize',plot_struct.fsize)
title('4 Stations','fontsize',plot_struct.fsize)
%set(gca,'position',[ax_start+2*awidth+2*hgap, ay_start+2*aheight+2*vgap, awidth, aheight])
set(gca,'xticklabel',[],'yticklabel',[])

% CWI x vs y (4 stations)
sh(4) = subplot(3,6,4);
plot(CatLoc_not68(:,5), CatLoc_not68(:,6),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_CWI_rs.(statlist{4}).relocx, Locs_CWI_rs.(statlist{4}).relocy,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.ylimits,'units','centimeters')
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ytickspots)
set(gca,'fontsize',plot_struct.fsize)
title('3 Stations','fontsize',plot_struct.fsize)
%set(gca,'position',[ax_start+2*awidth+2*hgap, ay_start+2*aheight+2*vgap, awidth, aheight])
set(gca,'xticklabel',[],'yticklabel',[])

% CWI x vs y (3 stations)
sh(5) = subplot(3,6,5);
plot(CatLoc_not68(:,5), CatLoc_not68(:,6),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_CWI_rs.(statlist{5}).relocx, Locs_CWI_rs.(statlist{5}).relocy,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.ylimits,'units','centimeters')
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ytickspots)
set(gca,'fontsize',plot_struct.fsize)
title('2 Stations','fontsize',plot_struct.fsize)
%set(gca,'position',[ax_start+2*awidth+2*hgap, ay_start+2*aheight+2*vgap, awidth, aheight])
set(gca,'xticklabel',[],'yticklabel',[])

% CWI x vs y (1 station)
sh(6) = subplot(3,6,6);
plot(CatLoc_not68(:,5), CatLoc_not68(:,6),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_CWI_rs.(statlist{6}).relocx, Locs_CWI_rs.(statlist{6}).relocy,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.ylimits,'units','centimeters')
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ytickspots)
set(gca,'fontsize',plot_struct.fsize)
title('1 Station','fontsize',plot_struct.fsize)
%set(gca,'position',[ax_start+2*awidth+2*hgap, ay_start+2*aheight+2*vgap, awidth, aheight])
set(gca,'xticklabel',[],'yticklabel',[])

% CWI x vs z (7 stations)
sh(7) = subplot(3,6,7);
plot(CatLoc_not68(:,5), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_CWI_rs.(statlist{1}).relocx, Locs_CWI_rs.(statlist{1}).relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
%xlabel('x (m)','fontsize',plot_struct.fsize)
%ylabel('z (m)','fontsize',plot_struct.fsize)
ylabel('x vs z (m)','fontsize',plot_struct.fsize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
%set(gca,'position',axes_posy)
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ztickspots,'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
%set(gca,'position',[ax_start, ay_start+aheight+vgap, awidth, aheight])
set(gca,'xticklabel',[])

% CWI x vs z (5 stations)
sh(8) = subplot(3,6,8);
plot(CatLoc_not68(:,5), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_CWI_rs.(statlist{2}).relocx, Locs_CWI_rs.(statlist{2}).relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ztickspots,'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
%set(gca,'position',[ax_start+awidth+hgap, ay_start+aheight+vgap, awidth, aheight])
set(gca,'xticklabel',[],'yticklabel',[])

% CWI x vs z (4 stations)
sh(9) = subplot(3,6,9);
plot(CatLoc_not68(:,5), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_CWI_rs.(statlist{3}).relocx, Locs_CWI_rs.(statlist{3}).relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ztickspots,'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
%set(gca,'position',[ax_start+2*awidth+2*hgap, ay_start+aheight+vgap, awidth, aheight])
set(gca,'xticklabel',[],'yticklabel',[])

% CWI x vs z (3 stations)
sh(10) = subplot(3,6,10);
plot(CatLoc_not68(:,5), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_CWI_rs.(statlist{4}).relocx, Locs_CWI_rs.(statlist{4}).relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ztickspots,'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
set(gca,'xticklabel',[],'yticklabel',[])

% CWI x vs z (2 stations)
sh(11) = subplot(3,6,11);
plot(CatLoc_not68(:,5), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_CWI_rs.(statlist{5}).relocx, Locs_CWI_rs.(statlist{5}).relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ztickspots,'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
set(gca,'xticklabel',[],'yticklabel',[])

% CWI x vs z (1 stations)
sh(12) = subplot(3,6,12);
plot(CatLoc_not68(:,5), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_CWI_rs.(statlist{6}).relocx, Locs_CWI_rs.(statlist{6}).relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ztickspots,'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
set(gca,'xticklabel',[],'yticklabel',[])

% CWI y vs z (7 stations)
sh(13) = subplot(3,6,13);
plot(CatLoc_not68(:,6), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_CWI_rs.(statlist{1}).relocy, Locs_CWI_rs.(statlist{1}).relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
%xlabel('y (m)','fontsize',plot_struct.fsize)
%ylabel('z (m)','fontsize',plot_struct.fsize)
ylabel('y vs z (m)','fontsize',plot_struct.fsize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
%set(gca,'position',axes_posy)
set(gca,'xtick', plot_struct.ytickspots, 'ytick',plot_struct.ztickspots, 'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
%set(gca,'position',[ax_start, ay_start, awidth, aheight])
%set(gca,'yticklabel',[])

% CWI y vs z (5 stations)
sh(14) = subplot(3,6,14);
plot(CatLoc_not68(:,6), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_CWI_rs.(statlist{2}).relocy, Locs_CWI_rs.(statlist{2}).relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
set(gca,'xtick', plot_struct.ytickspots, 'ytick',plot_struct.ztickspots, 'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
%set(gca,'position',[ax_start+awidth+hgap, ay_start, awidth, aheight])
set(gca,'yticklabel',[])

% CWI y vs z (4 stations)
sh(15) = subplot(3,6,15);
plot(CatLoc_not68(:,6), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_CWI_rs.(statlist{3}).relocy, Locs_CWI_rs.(statlist{3}).relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
set(gca,'xtick', plot_struct.ytickspots, 'ytick',plot_struct.ztickspots, 'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
%set(gca,'position',[ax_start+2*awidth+2*hgap, ay_start, awidth, aheight])
set(gca,'yticklabel',[])

% CWI y vs z (3 stations)
sh(16) = subplot(3,6,16);
plot(CatLoc_not68(:,6), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_CWI_rs.(statlist{4}).relocy, Locs_CWI_rs.(statlist{4}).relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
set(gca,'xtick', plot_struct.ytickspots, 'ytick',plot_struct.ztickspots, 'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
set(gca,'yticklabel',[])

% CWI y vs z (2 stations)
sh(17) = subplot(3,6,17);
plot(CatLoc_not68(:,6), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_CWI_rs.(statlist{5}).relocy, Locs_CWI_rs.(statlist{5}).relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
set(gca,'xtick', plot_struct.ytickspots, 'ytick',plot_struct.ztickspots, 'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
set(gca,'yticklabel',[])

% CWI y vs z (1 station)
sh(18) = subplot(3,6,18);
plot(CatLoc_not68(:,6), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_CWI_rs.(statlist{6}).relocy, Locs_CWI_rs.(statlist{6}).relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
set(gca,'xtick', plot_struct.ytickspots, 'ytick',plot_struct.ztickspots, 'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
set(gca,'yticklabel',[])


set(sh(1),'position',[ax_start, ay_start+2*aheight+2*vgap, awidth, aheight])
set(sh(2),'position',[ax_start+awidth+hgap, ay_start+2*aheight+2*vgap, awidth, aheight])
set(sh(3),'position',[ax_start+2*awidth+2*hgap, ay_start+2*aheight+2*vgap, awidth, aheight])
set(sh(4),'position',[ax_start+3*awidth+3*hgap, ay_start+2*aheight+2*vgap, awidth, aheight])
set(sh(5),'position',[ax_start+4*awidth+4*hgap, ay_start+2*aheight+2*vgap, awidth, aheight])
set(sh(6),'position',[ax_start+5*awidth+5*hgap, ay_start+2*aheight+2*vgap, awidth, aheight])

set(sh(7),'position',[ax_start, ay_start+aheight+vgap, awidth, aheight])
set(sh(8),'position',[ax_start+awidth+hgap, ay_start+aheight+vgap, awidth, aheight])
set(sh(9),'position',[ax_start+2*awidth+2*hgap, ay_start+aheight+vgap, awidth, aheight])
set(sh(10),'position',[ax_start+3*awidth+3*hgap, ay_start+aheight+vgap, awidth, aheight])
set(sh(11),'position',[ax_start+4*awidth+4*hgap, ay_start+aheight+vgap, awidth, aheight])
set(sh(12),'position',[ax_start+5*awidth+5*hgap, ay_start+aheight+vgap, awidth, aheight])

set(sh(13),'position',[ax_start, ay_start, awidth, aheight])
set(sh(14),'position',[ax_start+awidth+hgap, ay_start, awidth, aheight])
set(sh(15),'position',[ax_start+2*awidth+2*hgap, ay_start, awidth, aheight])
set(sh(16),'position',[ax_start+3*awidth+3*hgap, ay_start, awidth, aheight])
set(sh(17),'position',[ax_start+4*awidth+4*hgap, ay_start, awidth, aheight])
set(sh(18),'position',[ax_start+5*awidth+5*hgap, ay_start, awidth, aheight])

set(gcf,'units','centimeters')
set(gcf,'position',[5, 5, 26, 12.5],'PaperOrientation','Portrait')
set(gcf,'PaperPosition',[0.6345    6.3452   27   15.2284])

print -depsc CalaverasLoc2.eps


%% Figure 3 - HYPODD reduced stations
figure
% 3x6 figure showing:
%       col1 => CWI with 7 stations
%       col2 => CWI with 5 stations
%       col3 => CWI with 4 stations
%       col4 => CWI with 3 stations
%       col5 => CWI with 2 stations
%       col6 => CWI with 1 stations

% Lets get the results we need for this figure
statlist = {'stat7','stat5','stat4','stat3','stat2','stat1'};
for j = 1: length(statlist)
    if strcmp(filesep,'/')
        dirname = ['/export/storage/davidr/sandpit/davidr/thesis_version2/diags/eq_location_optimisation/CalaverasMultiSim/HYPODD2/runs/',statlist{j},'/'];
    else
        dirname = ['c:/datafiles/workstuff/sandpit/davidr/thesis_version2/diags/eq_location_optimisation/CalaverasMultiSim/HYPODD2/runs/',statlist{j},'/'];
    end
    [tmpstats,Locs_HYPODD_rs.(statlist{j}),fh] = plot_outcomes_hypoDD(dirname,plot_struct);
    for i = 1:length(fh)
        close(fh(i))
    end
end

for i = 1:18
    subplot(3,6,i)
end
sh = zeros(1,18);

% CWI x vs y (7 stations)
sh(1) = subplot(3,6,1);
plot(CatLoc_not68(:,5), CatLoc_not68(:,6),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_HYPODD_rs.(statlist{1}).relocx, Locs_HYPODD_rs.(statlist{1}).relocy,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
%xlabel('x (m)','fontsize',plot_struct.fsize)
%ylabel('y (m)','fontsize',plot_struct.fsize)
ylabel('x vs y (m)','fontsize',plot_struct.fsize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.ylimits,'units','centimeters')
%set(gca,'position',axes_posy)
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ytickspots)
set(gca,'fontsize',plot_struct.fsize)
title('7 Stations','fontsize',plot_struct.fsize)
%set(gca,'position',[ax_start, ay_start+2*aheight+2*vgap, awidth, aheight])
set(gca,'xticklabel',[])

% CWI x vs y (5 stations)
sh(2) = subplot(3,6,2);
plot(CatLoc_not68(:,5), CatLoc_not68(:,6),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_HYPODD_rs.(statlist{2}).relocx, Locs_HYPODD_rs.(statlist{2}).relocy,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.ylimits,'units','centimeters')
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ytickspots)
set(gca,'fontsize',plot_struct.fsize)
title('5 Stations','fontsize',plot_struct.fsize)
%set(gca,'position',[ax_start+awidth+hgap, ay_start+2*aheight+2*vgap, awidth, aheight])
set(gca,'xticklabel',[],'yticklabel',[])

% CWI x vs y (4 stations)
sh(3) = subplot(3,6,3);
plot(CatLoc_not68(:,5), CatLoc_not68(:,6),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_HYPODD_rs.(statlist{3}).relocx, Locs_HYPODD_rs.(statlist{3}).relocy,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.ylimits,'units','centimeters')
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ytickspots)
set(gca,'fontsize',plot_struct.fsize)
title('4 Stations','fontsize',plot_struct.fsize)
%set(gca,'position',[ax_start+2*awidth+2*hgap, ay_start+2*aheight+2*vgap, awidth, aheight])
set(gca,'xticklabel',[],'yticklabel',[])

% CWI x vs y (4 stations)
sh(4) = subplot(3,6,4);
plot(CatLoc_not68(:,5), CatLoc_not68(:,6),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_HYPODD_rs.(statlist{4}).relocx, Locs_HYPODD_rs.(statlist{4}).relocy,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.ylimits,'units','centimeters')
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ytickspots)
set(gca,'fontsize',plot_struct.fsize)
title('3 Stations','fontsize',plot_struct.fsize)
%set(gca,'position',[ax_start+2*awidth+2*hgap, ay_start+2*aheight+2*vgap, awidth, aheight])
set(gca,'xticklabel',[],'yticklabel',[])

% CWI x vs y (3 stations)
sh(5) = subplot(3,6,5);
plot(CatLoc_not68(:,5), CatLoc_not68(:,6),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_HYPODD_rs.(statlist{5}).relocx, Locs_HYPODD_rs.(statlist{5}).relocy,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.ylimits,'units','centimeters')
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ytickspots)
set(gca,'fontsize',plot_struct.fsize)
title('2 Stations','fontsize',plot_struct.fsize)
%set(gca,'position',[ax_start+2*awidth+2*hgap, ay_start+2*aheight+2*vgap, awidth, aheight])
set(gca,'xticklabel',[],'yticklabel',[])

% CWI x vs y (1 station)
sh(6) = subplot(3,6,6);
plot(CatLoc_not68(:,5), CatLoc_not68(:,6),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_HYPODD_rs.(statlist{6}).relocx, Locs_HYPODD_rs.(statlist{6}).relocy,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.ylimits,'units','centimeters')
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ytickspots)
set(gca,'fontsize',plot_struct.fsize)
title('1 Station','fontsize',plot_struct.fsize)
%set(gca,'position',[ax_start+2*awidth+2*hgap, ay_start+2*aheight+2*vgap, awidth, aheight])
set(gca,'xticklabel',[],'yticklabel',[])

% CWI x vs z (7 stations)
sh(7) = subplot(3,6,7);
plot(CatLoc_not68(:,5), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_HYPODD_rs.(statlist{1}).relocx, Locs_HYPODD_rs.(statlist{1}).relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
%xlabel('x (m)','fontsize',plot_struct.fsize)
%ylabel('z (m)','fontsize',plot_struct.fsize)
ylabel('x vs z (m)','fontsize',plot_struct.fsize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
%set(gca,'position',axes_posy)
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ztickspots,'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
%set(gca,'position',[ax_start, ay_start+aheight+vgap, awidth, aheight])
set(gca,'xticklabel',[])

% CWI x vs z (5 stations)
sh(8) = subplot(3,6,8);
plot(CatLoc_not68(:,5), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_HYPODD_rs.(statlist{2}).relocx, Locs_HYPODD_rs.(statlist{2}).relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ztickspots,'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
%set(gca,'position',[ax_start+awidth+hgap, ay_start+aheight+vgap, awidth, aheight])
set(gca,'xticklabel',[],'yticklabel',[])

% CWI x vs z (4 stations)
sh(9) = subplot(3,6,9);
plot(CatLoc_not68(:,5), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_HYPODD_rs.(statlist{3}).relocx, Locs_HYPODD_rs.(statlist{3}).relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ztickspots,'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
%set(gca,'position',[ax_start+2*awidth+2*hgap, ay_start+aheight+vgap, awidth, aheight])
set(gca,'xticklabel',[],'yticklabel',[])

% CWI x vs z (3 stations)
sh(10) = subplot(3,6,10);
plot(CatLoc_not68(:,5), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_HYPODD_rs.(statlist{4}).relocx, Locs_HYPODD_rs.(statlist{4}).relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ztickspots,'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
set(gca,'xticklabel',[],'yticklabel',[])

% CWI x vs z (2 stations)
sh(11) = subplot(3,6,11);
plot(CatLoc_not68(:,5), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_HYPODD_rs.(statlist{5}).relocx, Locs_HYPODD_rs.(statlist{5}).relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ztickspots,'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
set(gca,'xticklabel',[],'yticklabel',[])

% CWI x vs z (1 stations)
sh(12) = subplot(3,6,12);
plot(CatLoc_not68(:,5), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_HYPODD_rs.(statlist{6}).relocx, Locs_HYPODD_rs.(statlist{6}).relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ztickspots,'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
set(gca,'xticklabel',[],'yticklabel',[])

% CWI y vs z (7 stations)
sh(13) = subplot(3,6,13);
plot(CatLoc_not68(:,6), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_HYPODD_rs.(statlist{1}).relocy, Locs_HYPODD_rs.(statlist{1}).relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
%xlabel('y (m)','fontsize',plot_struct.fsize)
%ylabel('z (m)','fontsize',plot_struct.fsize)
ylabel('y vs z (m)','fontsize',plot_struct.fsize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
%set(gca,'position',axes_posy)
set(gca,'xtick', plot_struct.ytickspots, 'ytick',plot_struct.ztickspots, 'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
%set(gca,'position',[ax_start, ay_start, awidth, aheight])
%set(gca,'xticklabel',[],'yticklabel',[])

% CWI y vs z (5 stations)
sh(14) = subplot(3,6,14);
plot(CatLoc_not68(:,6), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_HYPODD_rs.(statlist{2}).relocy, Locs_HYPODD_rs.(statlist{2}).relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
set(gca,'xtick', plot_struct.ytickspots, 'ytick',plot_struct.ztickspots, 'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
%set(gca,'position',[ax_start+awidth+hgap, ay_start, awidth, aheight])
set(gca,'yticklabel',[])

% CWI y vs z (4 stations)
sh(15) = subplot(3,6,15);
plot(CatLoc_not68(:,6), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_HYPODD_rs.(statlist{3}).relocy, Locs_HYPODD_rs.(statlist{3}).relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
set(gca,'xtick', plot_struct.ytickspots, 'ytick',plot_struct.ztickspots, 'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
%set(gca,'position',[ax_start+2*awidth+2*hgap, ay_start, awidth, aheight])
set(gca,'yticklabel',[])

% CWI y vs z (3 stations)
sh(16) = subplot(3,6,16);
plot(CatLoc_not68(:,6), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_HYPODD_rs.(statlist{4}).relocy, Locs_HYPODD_rs.(statlist{4}).relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
set(gca,'xtick', plot_struct.ytickspots, 'ytick',plot_struct.ztickspots, 'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
set(gca,'yticklabel',[])

% CWI y vs z (2 stations)
sh(17) = subplot(3,6,17);
plot(CatLoc_not68(:,6), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_HYPODD_rs.(statlist{5}).relocy, Locs_HYPODD_rs.(statlist{5}).relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
set(gca,'xtick', plot_struct.ytickspots, 'ytick',plot_struct.ztickspots, 'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
set(gca,'yticklabel',[])

% CWI y vs z (1 station)
sh(18) = subplot(3,6,18);
plot(CatLoc_not68(:,6), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_HYPODD_rs.(statlist{6}).relocy, Locs_HYPODD_rs.(statlist{6}).relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
set(gca,'xtick', plot_struct.ytickspots, 'ytick',plot_struct.ztickspots, 'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
set(gca,'yticklabel',[])


set(sh(1),'position',[ax_start, ay_start+2*aheight+2*vgap, awidth, aheight])
set(sh(2),'position',[ax_start+awidth+hgap, ay_start+2*aheight+2*vgap, awidth, aheight])
set(sh(3),'position',[ax_start+2*awidth+2*hgap, ay_start+2*aheight+2*vgap, awidth, aheight])
set(sh(4),'position',[ax_start+3*awidth+3*hgap, ay_start+2*aheight+2*vgap, awidth, aheight])
set(sh(5),'position',[ax_start+4*awidth+4*hgap, ay_start+2*aheight+2*vgap, awidth, aheight])
set(sh(6),'position',[ax_start+5*awidth+5*hgap, ay_start+2*aheight+2*vgap, awidth, aheight])

set(sh(7),'position',[ax_start, ay_start+aheight+vgap, awidth, aheight])
set(sh(8),'position',[ax_start+awidth+hgap, ay_start+aheight+vgap, awidth, aheight])
set(sh(9),'position',[ax_start+2*awidth+2*hgap, ay_start+aheight+vgap, awidth, aheight])
set(sh(10),'position',[ax_start+3*awidth+3*hgap, ay_start+aheight+vgap, awidth, aheight])
set(sh(11),'position',[ax_start+4*awidth+4*hgap, ay_start+aheight+vgap, awidth, aheight])
set(sh(12),'position',[ax_start+5*awidth+5*hgap, ay_start+aheight+vgap, awidth, aheight])

set(sh(13),'position',[ax_start, ay_start, awidth, aheight])
set(sh(14),'position',[ax_start+awidth+hgap, ay_start, awidth, aheight])
set(sh(15),'position',[ax_start+2*awidth+2*hgap, ay_start, awidth, aheight])
set(sh(16),'position',[ax_start+3*awidth+3*hgap, ay_start, awidth, aheight])
set(sh(17),'position',[ax_start+4*awidth+4*hgap, ay_start, awidth, aheight])
set(sh(18),'position',[ax_start+5*awidth+5*hgap, ay_start, awidth, aheight])

set(gcf,'units','centimeters')
set(gcf,'position',[5, 5, 26, 12.5],'PaperOrientation','Portrait')
set(gcf,'PaperPosition',[0.6345    6.3452   27   15.2284])

print -depsc CalaverasLoc3.eps

%% Figure 4 - statistical comparison of CWI and Travel Times
station_counts = [10,9,8,7,6,5,4,3,2,1];
%station_counts = 1; 
if strcmp(filesep,'/')
    parentdir = ['/export/storage/davidr/sandpit/davidr/thesis_version2/diags/eq_location_optimisation/CalaverasMultiSim/'];
else
    parentdir = ['c:/datafiles/workstuff/sandpit/davidr/thesis_version2/diags/eq_location_optimisation/CalaverasMultiSim/'];
end

for i = 1: length(station_counts)
    dirname2 = [parentdir,filesep, 'CWI/stat',num2str(station_counts(i))];
    [CWIstats_tmp, tmp, fh] = plot_outcomes_cwi(dirname2,plot_struct);
    for k = 1:length(fh)
        close(fh(k))
    end
    CWIstats_minval(i) = CWIstats_tmp.minval;
    CWIstats_maxval(i) = CWIstats_tmp.maxval;
    CWIstats_meanval(i) = CWIstats_tmp.meanval;
    CWIstats_stdval(i) = CWIstats_tmp.stdval;
    CWI_stats_rms(i) = CWIstats_tmp.plane.rms;
    CWIstats_nE(i) = CWIstats_tmp.nE;
    %pause
    dirname3 = [parentdir,filesep, 'HYPODD2/runs/stat',num2str(station_counts(i))];
    [hypoDD2stats_tmp, tmp, fh] = plot_outcomes_hypoDD(dirname3,plot_struct);
    for k = 1:length(fh)
        close(fh(k))
    end
    hypoDD2stats_minval(i) = hypoDD2stats_tmp.minval;
    hypoDD2stats_maxval(i) = hypoDD2stats_tmp.maxval;
    hypoDD2stats_meanval(i) = hypoDD2stats_tmp.meanval;
    hypoDD2stats_stdval(i) = hypoDD2stats_tmp.stdval;
    hypoDD2stats_rms(i) = hypoDD2stats_tmp.plane.rms;
    hypoDD2stats_nE(i) = hypoDD2stats_tmp.nE;
    hypoDD2stats_translation(i) = sqrt(sum(hypoDD2stats_tmp.translation.^2));

    

end

% ======================================================
% Plot the statistics for the updated HypoDD2 runs
sh = [];
figure
sh(1) = subplot(3,1,1); 
hCWI = plot(station_counts,CWIstats_nE,'color',c2,'linewidth',2);
hold on
hhypoDD2 = plot(station_counts,hypoDD2stats_nE,'color', c1,'linewidth',2);
set(gca,'xdir','reverse')
%xlabel('Number of stations')
ylabel('nE')
set(gca,'units','centimeters','xlim',[1,10])

sh(2) = subplot(3,1,2); 
hCWI = plot(station_counts,CWIstats_meanval,'color',c2,'linewidth',2);
hold on
hhypoDD2 = plot(station_counts,hypoDD2stats_meanval,'color', c1,'linewidth',2);
set(gca,'xdir','reverse')
legend([hCWI,hhypoDD2],{'CWI','Travel Times'},'Location','NorthWest')
ylabel('mean')
set(gca,'units','centimeters','xlim',[1,10])

sh(3) = subplot(3,1,3);
hCWI = plot(station_counts,CWIstats_maxval,'color',c2,'linewidth',2);
hold on
hhypoDD2 = plot(station_counts,hypoDD2stats_maxval,'color', c1,'linewidth',2);
set(gca,'xdir','reverse')
xlabel('Number of stations')
ylabel('max')
set(gca,'units','centimeters','xlim',[1,10])

statwidth = 7.5;
statheight = 1.5;
statxstart = 4.6;
statystart = 1.95;
statygap = 0.7;

set(sh(3), 'position',[statxstart, statystart, statwidth, statheight])
set(sh(2), 'position',[statxstart, statystart+statheight+statygap, statwidth, statheight])
set(sh(1), 'position',[statxstart, statystart+2*statheight+2*statygap, statwidth, statheight])

print -depsc CalaverasLoc4.eps

%% Figure 5 - Combined CWI & Travel Time - Low/Med/High Uncert on tt
% 3x3 figures showing:
%       col1 => Low uncert
%       col2 => Med Uncert
%       col3 => High Uncert

% Lets get the results we need for this figure
if strcmp(filesep,'/')
    parentdir = '/export/storage/davidr/sandpit/davidr/thesis_version2/diags/eq_location_optimisation/CalaverasCartesian2/';
else
    parentdir = 'c:/datafiles/workstuff/sandpit/davidr/thesis_version2/diags/eq_location_optimisation/CalaverasCartesian2/';
end

[CWIstats_tmp,CombLocs_Low,fh] = plot_outcomes_cwi_cartesian([parentdir,'example_68eq_calaveras1/'],plot_struct)
for i = 1:length(fh)
    close(fh(i))
end
[CWIstats_tmp,CombLocs_Med,fh] = plot_outcomes_cwi_cartesian([parentdir,'example_68eq_calaveras2/'],plot_struct)
for i = 1:length(fh)
    close(fh(i))
end
[CWIstats_tmp,CombLocs_High,fh] = plot_outcomes_cwi_cartesian([parentdir,'example_68eq_calaveras3/'],plot_struct)
for i = 1:length(fh)
    close(fh(i))
end

figure
set(gcf,'units','centimeter')

% catalogue x vs y
subplot(3,3,1)
plot(Locs_hypoDDbest.not68x, Locs_hypoDDbest.not68y,'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(CombLocs_Low.relocx, CombLocs_Low.relocy,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
%xlabel('x (m)','fontsize',plot_struct.fsize)
%ylabel('y (m)','fontsize',plot_struct.fsize)
ylabel('x vs y (m)','fontsize',plot_struct.fsize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.ylimits,'units','centimeters')
%set(gca,'position',axes_posy)
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ytickspots)
set(gca,'fontsize',plot_struct.fsize)
title('Small','fontsize',plot_struct.fsize)
set(gca,'position',[ax_start, ay_start+2*aheight+2*vgap, awidth, aheight])
set(gca,'xticklabel',[])

% CWI x vs y
subplot(3,3,2)
plot(Locs_hypoDDbest.not68x, Locs_hypoDDbest.not68y,'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(CombLocs_Med.relocx, CombLocs_Med.relocy,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.ylimits,'units','centimeters')
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ytickspots)
set(gca,'fontsize',plot_struct.fsize)
title('Medium','fontsize',plot_struct.fsize)
set(gca,'position',[ax_start+awidth+hgap, ay_start+2*aheight+2*vgap, awidth, aheight])
set(gca,'xticklabel',[],'yticklabel',[])

% hypoDD x vs y
subplot(3,3,3)
plot(Locs_hypoDDbest.not68x, Locs_hypoDDbest.not68y,'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(CombLocs_High.relocx, CombLocs_High.relocy,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.ylimits,'units','centimeters')
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ytickspots)
set(gca,'fontsize',plot_struct.fsize)
title('Large','fontsize',plot_struct.fsize)
set(gca,'position',[ax_start+2*awidth+2*hgap, ay_start+2*aheight+2*vgap, awidth, aheight])
set(gca,'xticklabel',[],'yticklabel',[])

% catalogue x vs z
subplot(3,3,4)
plot(Locs_hypoDDbest.not68x, Locs_hypoDDbest.not68z,'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(CombLocs_Low.relocx, CombLocs_Low.relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
%xlabel('x (m)','fontsize',plot_struct.fsize)
%ylabel('z (m)','fontsize',plot_struct.fsize)
ylabel('x vs z (m)','fontsize',plot_struct.fsize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
%set(gca,'position',axes_posy)
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ztickspots,'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
set(gca,'position',[ax_start, ay_start+aheight+vgap, awidth, aheight])
set(gca,'xticklabel',[])

% CWI x vs z
subplot(3,3,5)
plot(Locs_hypoDDbest.not68x, Locs_hypoDDbest.not68z,'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(CombLocs_Med.relocx, CombLocs_Med.relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ztickspots,'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
set(gca,'position',[ax_start+awidth+hgap, ay_start+aheight+vgap, awidth, aheight])
set(gca,'xticklabel',[],'yticklabel',[])

% hypoDD x vs z
subplot(3,3,6)
plot(Locs_hypoDDbest.not68x, Locs_hypoDDbest.not68z,'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(CombLocs_High.relocx, CombLocs_High.relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ztickspots,'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
set(gca,'position',[ax_start+2*awidth+2*hgap, ay_start+aheight+vgap, awidth, aheight])
set(gca,'xticklabel',[],'yticklabel',[])

% catalogue y vs z
subplot(3,3,7)
plot(Locs_hypoDDbest.not68y, Locs_hypoDDbest.not68z,'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(CombLocs_Low.relocy, CombLocs_Low.relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
%xlabel('y (m)','fontsize',plot_struct.fsize)
%ylabel('z (m)','fontsize',plot_struct.fsize)
ylabel('y vs z (m)','fontsize',plot_struct.fsize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
%set(gca,'position',axes_posy)
set(gca,'xtick', plot_struct.ytickspots, 'ytick',plot_struct.ztickspots, 'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
set(gca,'position',[ax_start, ay_start, awidth, aheight])

% CWI y vs z
subplot(3,3,8)
plot(Locs_hypoDDbest.not68y, Locs_hypoDDbest.not68z,'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(CombLocs_Med.relocy, CombLocs_Med.relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
set(gca,'xtick', plot_struct.ytickspots, 'ytick',plot_struct.ztickspots, 'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
set(gca,'position',[ax_start+awidth+hgap, ay_start, awidth, aheight])
set(gca,'yticklabel',[])

% hypoDD y vs z
subplot(3,3,9)
plot(Locs_hypoDDbest.not68y, Locs_hypoDDbest.not68z,'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(CombLocs_High.relocy, CombLocs_High.relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
set(gca,'xtick', plot_struct.ytickspots, 'ytick',plot_struct.ztickspots, 'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
set(gca,'position',[ax_start+2*awidth+2*hgap, ay_start, awidth, aheight])
set(gca,'yticklabel',[])

print -depsc CalaverasLoc5.eps



%% Figure 6 - Combining deployment types (network and temporary)
%  2x3 figures showing:
%       col1 => Combined CWI & Travel Time (subset)- 
%       col2 => Travel time only 

% Lets get the results we need for this figure
if strcmp(filesep,'/')
    parentdir = '/export/storage/davidr/sandpit/davidr/thesis_version2/diags/eq_location_optimisation/CalaverasCartesian2/example_68eq_calaveras7';
else
    parentdir = 'c:/datafiles/workstuff/sandpit/davidr/thesis_version2/diags/eq_location_optimisation/CalaverasCartesian2/example_68eq_calaveras7';
end
[CWIstats_tmp,CombLocs_CWI_2deploy,fh] = plot_outcomes_cwi_cartesian(parentdir,plot_struct)
for i = 1:length(fh)
    close(fh(i))
end

parentdir='./extra_hypoDD/';
[HYPODDstats_tmp,CombLocs_HYPODD_2deploy,fh] = plot_outcomes_hypoDD(parentdir,plot_struct)
for i = 1:length(fh)
    close(fh(i))
end


figure
set(gcf,'units','centimeter')

% catalogue x vs y
subplot(3,2,1)
plot(Locs_hypoDDbest.not68x, Locs_hypoDDbest.not68y,'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(CombLocs_CWI_2deploy.relocx, CombLocs_CWI_2deploy.relocy,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
%xlabel('x (m)','fontsize',plot_struct.fsize)
%ylabel('y (m)','fontsize',plot_struct.fsize)
ylabel('x vs y (m)','fontsize',plot_struct.fsize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.ylimits,'units','centimeters')
%set(gca,'position',axes_posy)
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ytickspots)
set(gca,'fontsize',plot_struct.fsize)
title('CWI + TT','fontsize',plot_struct.fsize)
set(gca,'position',[ax_start, ay_start+2*aheight+2*vgap, awidth, aheight])
set(gca,'xticklabel',[])

% CWI x vs y
subplot(3,2,2)
plot(Locs_hypoDDbest.not68x, Locs_hypoDDbest.not68y,'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(CombLocs_HYPODD_2deploy.relocx, CombLocs_HYPODD_2deploy.relocy,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.ylimits,'units','centimeters')
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ytickspots)
set(gca,'fontsize',plot_struct.fsize)
title('TT only','fontsize',plot_struct.fsize)
set(gca,'position',[ax_start+awidth+hgap, ay_start+2*aheight+2*vgap, awidth, aheight])
set(gca,'xticklabel',[],'yticklabel',[])

% catalogue x vs z
subplot(3,2,3)
plot(Locs_hypoDDbest.not68x, Locs_hypoDDbest.not68z,'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(CombLocs_CWI_2deploy.relocx, CombLocs_CWI_2deploy.relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
%xlabel('x (m)','fontsize',plot_struct.fsize)
%ylabel('z (m)','fontsize',plot_struct.fsize)
ylabel('x vs z (m)','fontsize',plot_struct.fsize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
%set(gca,'position',axes_posy)
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ztickspots,'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
set(gca,'position',[ax_start, ay_start+aheight+vgap, awidth, aheight])
set(gca,'xticklabel',[])

% CWI x vs z
subplot(3,2,4)
plot(Locs_hypoDDbest.not68x, Locs_hypoDDbest.not68z,'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(CombLocs_HYPODD_2deploy.relocx, CombLocs_HYPODD_2deploy.relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ztickspots,'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
set(gca,'position',[ax_start+awidth+hgap, ay_start+aheight+vgap, awidth, aheight])
set(gca,'xticklabel',[],'yticklabel',[])

% catalogue y vs z
subplot(3,2,5)
plot(Locs_hypoDDbest.not68y, Locs_hypoDDbest.not68z,'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(CombLocs_CWI_2deploy.relocy, CombLocs_CWI_2deploy.relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
%xlabel('y (m)','fontsize',plot_struct.fsize)
%ylabel('z (m)','fontsize',plot_struct.fsize)
ylabel('y vs z (m)','fontsize',plot_struct.fsize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
%set(gca,'position',axes_posy)
set(gca,'xtick', plot_struct.ytickspots, 'ytick',plot_struct.ztickspots, 'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
set(gca,'position',[ax_start, ay_start, awidth, aheight])

% CWI y vs z
subplot(2,3,6)
plot(Locs_hypoDDbest.not68y, Locs_hypoDDbest.not68z,'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(CombLocs_HYPODD_2deploy.relocy, CombLocs_HYPODD_2deploy.relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
set(gca,'xtick', plot_struct.ytickspots, 'ytick',plot_struct.ztickspots, 'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
set(gca,'position',[ax_start+awidth+hgap, ay_start, awidth, aheight])
set(gca,'yticklabel',[])

print -depsc CalaverasLoc6.eps






% dirname = '/export/storage/davidr/sandpit/davidr/thesis_version2/diags/eq_location_optimisation/CalaverasMultiSim/HYPODD/runs/orig_example2';
% [hypoDDstats_best,Locs_hypoDDbest,fh] = plot_outcomes_hypoDD(dirname,plot_struct,colorswitch)
% close(fh(1))
% close(fh(2))
% close(fh(3))


