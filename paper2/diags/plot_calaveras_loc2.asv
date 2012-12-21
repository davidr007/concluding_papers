% The purpose of this script is to plot the Calaveras earthquake locations
% IMPORTANT NOTE
% .... This version (i.e. plot_calaveras_loc2.m utilises hypoDD with
% .... SVD ........

 
homeswitch =0; %1=> on home computer - 0=> on work PC
if homeswitch ==1
    path2thesis = '../../../thesis_version2/';
elseif homeswitch ==0
    path2thesis = '../../../../invert/davidr/thesis_version2/';
end

if homeswitch ==1
        thesisdiags = '../../../thesis_version2/';
elseif homeswitch ==0
        thesisdiags = '../../../../invert/davidr/thesis_version2/diags/eq_location_optimisation/';
end

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

events_oi = load([path2thesis,'/diags/eq_location_optimisation/events_oi.txt']);
CatLoc = load([path2thesis,'/diags/eq_location_optimisation/Calaveras_cat_locations.txt']); 
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
    dirname = [path2thesis,'/diags/eq_location_optimisation/'];
else
    dirname = [path2thesis,'/diags/eq_location_optimisation/'];
    %dirname = '../../../thesis_version2/diags/eq_location_optimisation/'
    %dirname_hypoDD_SVD = '..\..\..\concluding_papers\paper2\diags\runs2\Hypodd_SVD\HYPODD2\runs\orig_68eqonly\';
    dirname_hypoDD_SVD = '.\runs2\Hypodd_SVD\HYPODD2\runs\orig_68eqonly\';
end

% Get the original best hypoDD locations (i.e. those where we used LSQR)
tmphypoDD68 = load([dirname,'hypoDD_loc68.txt']);
Locs_hypoDDbest_LSQR.relocx = tmphypoDD68(:,5);
Locs_hypoDDbest_LSQR.relocy = tmphypoDD68(:,6);
Locs_hypoDDbest_LSQR.relocz = tmphypoDD68(:,7);

% Now let's get the version where we use hypoDD with SVD
% We translate these to the same 'local' coordinate systems as the LSQR example 
tmphypoDD68 = load([dirname_hypoDD_SVD,'hypoDD.reloc']);
Locs_hypoDDbest.relocx = tmphypoDD68(:,5)+(Locs_hypoDDbest_LSQR.relocx(1) -tmphypoDD68(1,5));
Locs_hypoDDbest.relocy = tmphypoDD68(:,6)+(Locs_hypoDDbest_LSQR.relocy(1) -tmphypoDD68(1,6));
Locs_hypoDDbest.relocz = tmphypoDD68(:,7)+(Locs_hypoDDbest_LSQR.relocz(1) -tmphypoDD68(1,7));
Locs_hypoDDbest.sigmax = tmphypoDD68(:,8);
Locs_hypoDDbest.sigmay = tmphypoDD68(:,9);
Locs_hypoDDbest.sigmaz = tmphypoDD68(:,10);

tmphypoDDnot68 = load([dirname,'hypodd_not68.txt']);
Locs_hypoDDbest.not68x = tmphypoDDnot68(:,5);
Locs_hypoDDbest.not68y = tmphypoDDnot68(:,6);
Locs_hypoDDbest.not68z = tmphypoDDnot68(:,7);
if strcmp(filesep,'/')
    dirname = [path2thesis,'/diags/eq_location_optimisation/example_69eq_calaveras3'];
else
    dirname = [path2thesis,'/diags/eq_location_optimisation/example_69eq_calaveras3'];
end
%dirname = '/export/storage/davidr/sandpit/davidr/thesis_version2/diags/eq_location_optimisation/CalaverasMultiSim/CWI/stat10';
[CWIstats_best,Locs_CWIbest,fh] = plot_outcomes_cwi(dirname,plot_struct,path2thesis);
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

%print -depsc CalaverasLoc1_hypoDD_SVD.eps
print -depsc ../Figure6_bw.eps

%% Figure 2 - CWI Reduces Stations
% NO NEED TO MAKE THIS AGAIN SINCE IT DOESN'T CHANGE
% 3x6 figure showing:
%       col1 => CWI with 7 stations
%       col2 => CWI with 5 stations
%       col3 => CWI with 4 stations
%       col4 => CWI with 3 stations
%       col5 => CWI with 2 stations
%       col6 => CWI with 1 stations

% Lets get the results we need for this figure

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
statlist = {'stat7','stat5','stat4'}%,'stat3','stat2','stat1'};
for j = 1: length(statlist)
    if strcmp(filesep,'/')
        dirname = ['/export/storage/davidr/sandpit/davidr/concluding_papers/paper2/diags/runs2/Hypodd_SVD/HYPODD2/runs/',statlist{j},'/'];
    else
        dirname = ['..\..\..\concluding_papers\paper2\diags\runs2\Hypodd_SVD\HYPODD2\runs\',statlist{j},'/'];
    end
    [tmpstats,Locs_HYPODD_rs.(statlist{j}),fh] = plot_outcomes_hypoDD(dirname,plot_struct,[],path2thesis);
    for i = 1:length(fh)
        close(fh(i))
    end
end

for i = 1:18
    subplot(3,6,i)
end
sh = zeros(1,18);

% CWI x vs y (7 stations)
sh(1) = subplot(3,3,1);
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
sh(2) = subplot(3,3,2);
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
sh(3) = subplot(3,3,3);
plot(CatLoc_not68(:,5), CatLoc_not68(:,6),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_HYPODD_rs.(statlist{3}).relocx, Locs_HYPODD_rs.(statlist{3}).relocy,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.ylimits,'units','centimeters')
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ytickspots)
set(gca,'fontsize',plot_struct.fsize)
title('4 Stations','fontsize',plot_struct.fsize)
%set(gca,'position',[ax_start+2*awidth+2*hgap, ay_start+2*aheight+2*vgap, awidth, aheight])
set(gca,'xticklabel',[],'yticklabel',[])



% CWI x vs z (7 stations)
sh(4) = subplot(3,3,4);
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
sh(5) = subplot(3,3,5);
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
sh(6) = subplot(3,3,6);
plot(CatLoc_not68(:,5), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_HYPODD_rs.(statlist{3}).relocx, Locs_HYPODD_rs.(statlist{3}).relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ztickspots,'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
%set(gca,'position',[ax_start+2*awidth+2*hgap, ay_start+aheight+vgap, awidth, aheight])
set(gca,'xticklabel',[],'yticklabel',[])


% CWI y vs z (7 stations)
sh(7) = subplot(3,3,7);
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
sh(8) = subplot(3,3,8);
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
sh(9) = subplot(3,3,9);
plot(CatLoc_not68(:,6), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(Locs_HYPODD_rs.(statlist{3}).relocy, Locs_HYPODD_rs.(statlist{3}).relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
set(gca,'xtick', plot_struct.ytickspots, 'ytick',plot_struct.ztickspots, 'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
%set(gca,'position',[ax_start+2*awidth+2*hgap, ay_start, awidth, aheight])
set(gca,'yticklabel',[])




set(sh(1),'position',[ax_start, ay_start+2*aheight+2*vgap, awidth, aheight])
set(sh(2),'position',[ax_start+awidth+hgap, ay_start+2*aheight+2*vgap, awidth, aheight])
set(sh(3),'position',[ax_start+2*awidth+2*hgap, ay_start+2*aheight+2*vgap, awidth, aheight])

set(sh(4),'position',[ax_start, ay_start+aheight+vgap, awidth, aheight])
set(sh(5),'position',[ax_start+awidth+hgap, ay_start+aheight+vgap, awidth, aheight])
set(sh(6),'position',[ax_start+2*awidth+2*hgap, ay_start+aheight+vgap, awidth, aheight])

set(sh(7),'position',[ax_start, ay_start, awidth, aheight])
set(sh(8),'position',[ax_start+awidth+hgap, ay_start, awidth, aheight])
set(sh(9),'position',[ax_start+2*awidth+2*hgap, ay_start, awidth, aheight])


set(gcf,'units','centimeters')
set(gcf,'position',[5, 5, 26, 12.5],'PaperOrientation','Portrait')
set(gcf,'PaperPosition',[0.6345    6.3452   27   15.2284])

%print -depsc CalaverasLoc3_hypoDD_SVD.eps
print -depsc ../Figure9_bw.eps


%% Figure 4 - statistical comparison of CWI and Travel Times
station_counts = [10,9,8,7,6,5,4,3,2,1];
%station_counts = 1; 



if strcmp(filesep,'/')
    parentdir_CWI = [path2thesis,'/diags/eq_location_optimisation/CalaverasMultiSim/'];
    parentdir_HYPODD = ['/export/storage/davidr/sandpit/davidr/concluding_papers/paper2/diags/runs2/Hypodd_SVD/HYPODD2/runs/'];
else
    parentdir_CWI = [path2thesis,'/diags/eq_location_optimisation/CalaverasMultiSim/'];
    parentdir_HYPODD = ['..\..\..\concluding_papers\paper2\diags\runs2\Hypodd_SVD\HYPODD2\runs\'];
end

for i = 1: length(station_counts)
    dirname2 = [parentdir_CWI,filesep, 'CWI/stat',num2str(station_counts(i))];
    [CWIstats_tmp, tmp, fh] = plot_outcomes_cwi(dirname2,plot_struct,path2thesis);
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
    dirname3 = [parentdir_HYPODD,filesep, 'stat',num2str(station_counts(i))];
    
    try, % if the hypodd.reloc exists
        [hypoDD2stats_tmp, hypoDD_DetailedLocs.(['stat',num2str(station_counts(i))]), fh] = plot_outcomes_hypoDD(dirname3,plot_struct,[],path2thesis);
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
    catch, % if hypodd.reloc does not exist
        hypoDD2stats_minval(i) = NaN;
        hypoDD2stats_maxval(i) = NaN;
        hypoDD2stats_meanval(i) = NaN;
        hypoDD2stats_stdval(i) = NaN;
        hypoDD2stats_rms(i) = NaN;
        hypoDD2stats_nE(i) = NaN;
        hypoDD2stats_translation(i) = NaN;
    end

end

% ======================================================
% Plot the statistics for the updated HypoDD2 runs
sh = [];
figure
sh(1) = subplot(3,1,1); 
hCWI_s1 = plot(station_counts,CWIstats_nE,'color',c2,'linewidth',2);
hold on
hhypoDD2_s1 = plot(station_counts,hypoDD2stats_nE,'color', c1,'linewidth',2);
set(gca,'xdir','reverse')
%xlabel('Number of stations')
ylabel('nE')
set(gca,'units','centimeters','xlim',[1,10])

sh(2) = subplot(3,1,2); 
hCWI_s2 = plot(station_counts,CWIstats_meanval,'color',c2,'linewidth',2);
hold on
hhypoDD2_s2 = plot(station_counts,hypoDD2stats_meanval,'color', c1,'linewidth',2);
set(gca,'xdir','reverse')
ylabel('$\Delta_\mu$\,(m)','Interpreter','LaTex')
set(gca,'units','centimeters','xlim',[1,10])

sh(3) = subplot(3,1,3);
hCWI_s3 = plot(station_counts,CWIstats_maxval,'color',c2,'linewidth',2);
hold on
hhypoDD2_s3 = plot(station_counts,hypoDD2stats_maxval,'color', c1,'linewidth',2);
set(gca,'xdir','reverse')
xlabel('Number of stations')
ylabel('$\Delta_{max}$\,(m)','Interpreter','LaTex')
legend([hCWI_s3,hhypoDD2_s3],{'CWI','TT'},'Location','NorthEast')
set(gca,'units','centimeters','xlim',[1,10],'ylim',[0 5000])

statwidth = 7.5;
statheight = 1.5;
statxstart = 4.6;
statystart = 1.95;
statygap = 0.7;

set(sh(3), 'position',[statxstart, statystart, statwidth, statheight])
set(sh(2), 'position',[statxstart, statystart+statheight+statygap, statwidth, statheight])
set(sh(1), 'position',[statxstart, statystart+2*statheight+2*statygap, statwidth, statheight])

%print -depsc CalaverasLoc4_hypoDD_SVD.eps
print -depsc ..\Figure10_bw.eps

set(hCWI_s1,'color','b')
set(hhypoDD2_s1,'color','r')
set(hCWI_s2,'color','b')
set(hhypoDD2_s2,'color','r')
set(hCWI_s3,'color','b')
set(hhypoDD2_s3,'color','r')
print -depsc ..\Figure10_c.eps

%%============
% The following was not working on 17 Dec so commented it out - DR
% Now lets have a look at the sigmas 
% figure
% for i = 1:7
%     statstr = ['stat', num2str(10-i+1)];
%     
%     subplot(7,3,(i-1)*3+1)
%     hist(hypoDD_DetailedLocs.(statstr).sigmax)
%     if i ==7
%        xlabel('\sigma_x (m)')
%     end
%     
%     subplot(7,3,(i-1)*3+2)
%     hist(hypoDD_DetailedLocs.(statstr).sigmay)
%     title(statstr)
%     if i ==7
%        xlabel('\sigma_y (m)')
%     end
%     
%     subplot(7,3,(i-1)*3+3)
%     hist(hypoDD_DetailedLocs.(statstr).sigmaz)
%     if i ==7
%        xlabel('\sigma_x (m)')
%     end
%     
% end
% print -depsc CalaverasLoc4b_hypoDD_SVD_sigmahists.eps



%% Figure 5 - Combined CWI & Travel Time - Low/Med/High Uncert on tt
% 3x3 figures showing:
%       col1 => Low uncert
%       col2 => Med Uncert
%       col3 => High Uncert

% Lets get the results we need for this figure

parentdir = 'runs2\Hypodd_SVD\HYPODD2\runs\';


[CWIstats_tmp,CombLocs_Low,fh] = plot_outcomes_cwi_cartesian([parentdir,'example_68eq_calaveras1/'],plot_struct,path2thesis)
for i = 1:length(fh)
    close(fh(i))
end
[CWIstats_tmp,CombLocs_Med,fh] = plot_outcomes_cwi_cartesian([parentdir,'example_68eq_calaveras2/'],plot_struct,path2thesis)
for i = 1:length(fh)
    close(fh(i))
end
[CWIstats_tmp,CombLocs_High,fh] = plot_outcomes_cwi_cartesian([parentdir,'example_68eq_calaveras2/'],plot_struct,path2thesis)
for i = 1:length(fh)
    close(fh(i))
end

figure
set(gcf,'units','centimeter')

% catalogue x vs y
subplot(3,2,1)
%plot(Locs_hypoDDbest.not68x, Locs_hypoDDbest.not68y,'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
plot(CatLoc_not68(:,5), CatLoc_not68(:,6),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(CombLocs_Low.relocx, CombLocs_Low.relocy,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
ylabel('x vs y (m)','fontsize',plot_struct.fsize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.ylimits,'units','centimeters')
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ytickspots)
set(gca,'fontsize',plot_struct.fsize)
title('All','fontsize',plot_struct.fsize)
set(gca,'position',[ax_start, ay_start+2*aheight+2*vgap, awidth, aheight])
set(gca,'xticklabel',[])

% CWI x vs y
subplot(3,2,2)
%plot(Locs_hypoDDbest.not68x, Locs_hypoDDbest.not68y,'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
plot(CatLoc_not68(:,5), CatLoc_not68(:,6),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(CombLocs_Med.relocx, CombLocs_Med.relocy,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.ylimits,'units','centimeters')
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ytickspots)
set(gca,'fontsize',plot_struct.fsize)
title('5 Stations','fontsize',plot_struct.fsize)
set(gca,'position',[ax_start+awidth+hgap, ay_start+2*aheight+2*vgap, awidth, aheight])
set(gca,'xticklabel',[],'yticklabel',[])

% hypoDD x vs y
% subplot(3,3,3)
% plot(Locs_hypoDDbest.not68x, Locs_hypoDDbest.not68y,'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
% hold on
% plot(CombLocs_High.relocx, CombLocs_High.relocy,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
% set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.ylimits,'units','centimeters')
% set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ytickspots)
% set(gca,'fontsize',plot_struct.fsize)
% title('7 Stations','fontsize',plot_struct.fsize)
% set(gca,'position',[ax_start+2*awidth+2*hgap, ay_start+2*aheight+2*vgap, awidth, aheight])
% set(gca,'xticklabel',[],'yticklabel',[])

% catalogue x vs z
subplot(3,2,3)
%plot(Locs_hypoDDbest.not68x, Locs_hypoDDbest.not68z,'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
plot(CatLoc_not68(:,5), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
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
subplot(3,2,4)
%plot(Locs_hypoDDbest.not68x, Locs_hypoDDbest.not68z,'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
plot(CatLoc_not68(:,5), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(CombLocs_Med.relocx, CombLocs_Med.relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ztickspots,'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
set(gca,'position',[ax_start+awidth+hgap, ay_start+aheight+vgap, awidth, aheight])
set(gca,'xticklabel',[],'yticklabel',[])

% hypoDD x vs z
% subplot(3,3,6)
% plot(Locs_hypoDDbest.not68x, Locs_hypoDDbest.not68z,'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
% hold on
% plot(CombLocs_High.relocx, CombLocs_High.relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
% set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
% set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ztickspots,'yticklabel',plot_struct.zticklabels)
% set(gca,'fontsize',plot_struct.fsize)
% set(gca,'units','centimeters')
% set(gca,'position',[ax_start+2*awidth+2*hgap, ay_start+aheight+vgap, awidth, aheight])
% set(gca,'xticklabel',[],'yticklabel',[])

% catalogue y vs z
subplot(3,2,5)
%plot(Locs_hypoDDbest.not68y, Locs_hypoDDbest.not68z,'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
plot(CatLoc_not68(:,6), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
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
subplot(3,2,6)
%plot(Locs_hypoDDbest.not68y, Locs_hypoDDbest.not68z,'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
plot(CatLoc_not68(:,6), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(CombLocs_Med.relocy, CombLocs_Med.relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
set(gca,'xtick', plot_struct.ytickspots, 'ytick',plot_struct.ztickspots, 'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
set(gca,'position',[ax_start+awidth+hgap, ay_start, awidth, aheight])
set(gca,'yticklabel',[])

% hypoDD y vs z
% subplot(3,3,9)
% plot(Locs_hypoDDbest.not68y, Locs_hypoDDbest.not68z,'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
% hold on
% plot(CombLocs_High.relocy, CombLocs_High.relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
% set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
% set(gca,'xtick', plot_struct.ytickspots, 'ytick',plot_struct.ztickspots, 'yticklabel',plot_struct.zticklabels)
% set(gca,'fontsize',plot_struct.fsize)
% set(gca,'units','centimeters')
% set(gca,'position',[ax_start+2*awidth+2*hgap, ay_start, awidth, aheight])
% set(gca,'yticklabel',[])

%print -depsc CalaverasLoc5_hypoDD_SVD.eps
print -depsc ../Figure11_bw.eps


%% Figure 6 - Combining deployment types (network and temporary)
%  2x3 figures showing:
%       col1 => Combined CWI & Travel Time (subset)- 
%       col2 => Travel time only 

% Lets get the results we need for this figure
if strcmp(filesep,'/')
    parentdir = './runs2\HYpodd_SVD\HyPODD2\runs\example_68eq_calaveras4';
else
    parentdir = './runs2\HYpodd_SVD\HyPODD2\runs\example_68eq_calaveras4';
end
[CWIstats_tmp,CombLocs_CWI_2deploy,fh] = plot_outcomes_cwi_cartesian(parentdir,plot_struct,path2thesis)
for i = 1:length(fh)
    close(fh(i))
end

parentdir='./extra_hypoDD3/';
[HYPODDstats_tmp,CombLocs_HYPODD_2deploy,fh] = plot_outcomes_hypoDD(parentdir,plot_struct,[],path2thesis)
for i = 1:length(fh)
    close(fh(i))
end


figure
set(gcf,'units','centimeter')
ay_start = 1.2
% catalogue x vs y
subplot(3,2,1)
%plot(Locs_hypoDDbest.not68x, Locs_hypoDDbest.not68y,'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
plot(CatLoc_not68(:,5), CatLoc_not68(:,6),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(CombLocs_CWI_2deploy.relocx, CombLocs_CWI_2deploy.relocy,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
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
%plot(Locs_hypoDDbest.not68x, Locs_hypoDDbest.not68y,'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
plot(CatLoc_not68(:,5), CatLoc_not68(:,6),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
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
%plot(Locs_hypoDDbest.not68x, Locs_hypoDDbest.not68z,'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
plot(CatLoc_not68(:,5), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(CombLocs_CWI_2deploy.relocx, CombLocs_CWI_2deploy.relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
ylabel('x vs z (m)','fontsize',plot_struct.fsize)
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
set(gca,'xtick', plot_struct.xtickspots, 'ytick',plot_struct.ztickspots,'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
set(gca,'position',[ax_start, ay_start+aheight+vgap, awidth, aheight])
set(gca,'xticklabel',[])

% CWI x vs z
subplot(3,2,4)
%plot(Locs_hypoDDbest.not68x, Locs_hypoDDbest.not68z,'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
plot(CatLoc_not68(:,5), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
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
%plot(Locs_hypoDDbest.not68y, Locs_hypoDDbest.not68z,'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
plot(CatLoc_not68(:,6), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(CombLocs_CWI_2deploy.relocy, CombLocs_CWI_2deploy.relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
%xlabel('y (m)','fontsize',plot_struct.fsize)
%ylabel('z (m)','fontsize',plot_struct.fsize)
ylabel('y vs z (m)','fontsize',plot_struct.fsize)
xlabel(['nE = ', num2str(length(CombLocs_CWI_2deploy.relocy))])
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
%set(gca,'position',axes_posy)
set(gca,'xtick', plot_struct.ytickspots, 'ytick',plot_struct.ztickspots, 'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
set(gca,'position',[ax_start, ay_start, awidth, aheight])

% CWI y vs z
subplot(2,3,6)
%plot(Locs_hypoDDbest.not68y, Locs_hypoDDbest.not68z,'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
plot(CatLoc_not68(:,6), CatLoc_not68(:,7),'o','markeredgecolor',c1,'markerfacecolor',c3,'markersize',plot_struct.msize)
hold on
plot(CombLocs_HYPODD_2deploy.relocy, CombLocs_HYPODD_2deploy.relocz,'o','markeredgecolor',c2,'markerfacecolor',c2,'markersize',plot_struct.msize)
xlabel(['nE = ', num2str(length(CombLocs_HYPODD_2deploy.relocy))])
set(gca,'xlim',plot_struct.xlimits, 'ylim',plot_struct.zlimits,'units','centimeters')
set(gca,'xtick', plot_struct.ytickspots, 'ytick',plot_struct.ztickspots, 'yticklabel',plot_struct.zticklabels)
set(gca,'fontsize',plot_struct.fsize)
set(gca,'units','centimeters')
%set(gca,'position',[ax_start+awidth+hgap, 2*hgap+ay_start, awidth, aheight])
set(gca,'position',[ax_start+awidth+hgap, ay_start, awidth, aheight])
set(gca,'yticklabel',[])

%print -depsc CalaverasLoc6_hypoDD_SVD.eps
print -depsc ../Figure13_bw.eps




% dirname = '/export/storage/davidr/sandpit/davidr/thesis_version2/diags/eq_location_optimisation/CalaverasMultiSim/HYPODD/runs/orig_example2';
% [hypoDDstats_best,Locs_hypoDDbest,fh] = plot_outcomes_hypoDD(dirname,plot_struct,colorswitch)
% close(fh(1))
% close(fh(2))
% close(fh(3))


