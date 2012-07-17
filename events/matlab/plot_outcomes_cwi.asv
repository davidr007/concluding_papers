function [CWIstats,Locs,fh] = plot_outcomes_cwi(dirname,plot_struct)
% This function plots the CWI outcomes. It is typically called by 
% batch_doplots. 
%
% INPUTS: 
% dirname       [char]  director where hypoDD outputs are stored. Note that 
%               eps figures are saved here as well. 
% plot_struct   [structure] settings for the figures. 
% 
% OUTPUTS: 
% CWIstats      [structure] contains statistics of interest
%                   CWIstats.minval = min(abs(delta_coord));
%                   CWIstats.maxval = max(abs(delta_coord));
%                   CWIstats.meanval = mean(abs(delta_coord));
%                   CWIstats.stdval = std(abs(delta_coord));
%                   CWIstats.nE = nE; 
%                   CWIstats.plane => info on best fitting plane;    
% Locs          [structure] contains the locations in desired coordinate system
%                   Locs.relocx = best_CWIx (xloc for optimum CWI inversion)
%                   Locs.relocy = best_CWIy (yloc for optimum CWI inversion)
%                   Locs.relocz = best_CWIz (zloc for optimum CWI inversion)
%                   Locs.not68x = hypodd_not68(i,5)
%                   Locs.not68y = hypodd_not68(i,6)
%                   Locs.not68z = hypodd_not68(i,7)
% fh            double [1Xn] figure handles

%
% David Robinson
% 9 December 2009

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
    parentdir = '/export/storage/davidr/sandpit/davidr/thesis_version2/diags/eq_location_optimisation/';
else
    parentdir = 'c:/datafiles/workstuff/sandpit/davidr/thesis_version2/diags/eq_location_optimisation/';
    %error('Not yet operation for the Windows machine')
end



opttype = 'fortran'; % 'fortran' or 'matlab'
load([dirname,filesep,'randsearch_',opttype,'.mat']);
[sorted_x2, ind]=sort(xend(2,:));

[nE,nrand] = size(xend);
npar = nE*3-3-2-1;

% the xend, yend and zend loaded above could be the old ones so let's
% reload and manipulate them here
xendvec_tmp = load([dirname,filesep,'xend.txt']);
xendvec = zeros(nrand, npar);
count = 3;
for i = 1:nrand
    xendvec(i,:) = xendvec_tmp(count:count+npar-1);
    count = count +npar;
end
ftmp = load([dirname,filesep,'fend.txt']);
f = ftmp(:,1);
 

% Unwrap the solution into separate matrices for x and y coordinates
xend_tmp = zeros(nE,nrand);
yend_tmp = zeros(nE,nrand);
zend_tmp = zeros(nE,nrand);
xend_tmp(2,:) = xendvec(:,1)';
xend_tmp(3,:) = xendvec(:,2)';
yend_tmp(3,:) = xendvec(:,3)';
count = 4;
for j = 4:nE
    xend_tmp(j,:) =  xendvec(:,count)';
    count = count+1;
    yend_tmp(j,:) =  xendvec(:,count)';
    count = count+1;
    zend_tmp(j,:) = xendvec(:,count)';
    count = count+1;
end

for i = 1:nrand
    % check for correct orientation of x variables (i.e. x2,x3,x4)
    if xend_tmp(2,i) < 0 % must flip the x variables
       xend_tmp(:,i) = -xend_tmp(:,i);
    end
    % check for correct orientation of y variables (i.e. y3,y4)
    if yend_tmp(3,i) < 0 % must flip the y variables
       yend_tmp(:,i) = -yend_tmp(:,i);
    end
    % check for correct orientation of z variables (i.e. z4)
    if zend_tmp(3,i) < 0  % must flip the z variables
       zend_tmp(:,i) = -zend_tmp(:,i);   % flip z4
    end
    
end



xend = zeros(nE,nrand);
yend = zeros(nE,nrand);
zend = zeros(nE,nrand);
for i = 1:nrand
    [cart_out_tmp] = convert_eventcoordinates(E(1,:),E(2,:),E(3,:), 'loc2cart', [xend_tmp(:,i),yend_tmp(:,i), zend_tmp(:,i)]);
    xend(:,i) = cart_out_tmp(:,1);
    yend(:,i) = cart_out_tmp(:,2);
    zend(:,i) = cart_out_tmp(:,3);
end

% % Have a quick look at the solution space
% figure
% colorlookup = {'b.','r.','g.','c.','m.','y.','k.','bo','ro','go'}
% colorlookup2 = {'b^','r^','g^','c^','m^','y^','k^','bd','rd','gd'}
% for i = 1:10
%     plot(xend(i,:), yend(i,:),colorlookup{i},'markersize',10)
%     hold on
%     plot(E(i,1),E(i,2), colorlookup2{i}, 'markersize',10)
% end
% plot([ 0 0], get(gca,'ylim'),'k')
% plot(get(gca,'xlim'),[0 0],'k')

[nevnts, nrand] = size(xend);

invert_stats = load([dirname,filesep,'fend.txt']);
ind2converged = [];
ind2notconverged = [];
% compute error as distancenorm
for k = 1:nrand
    count = 1;
   if invert_stats(k,2) ==1200  % i.e. inversion did not converge
       sq_seperror(k) = NaN;
       ind2notconverged = [ind2notconverged,k];
   else
       sq_seperror(k) = sum( sqrt( (E(:,1) - xend(:,k)).^2 + (E(:,2) - yend(:,k)).^2   + (E(:,3) - zend(:,k)).^2   ) );
       ind2converged = [ind2converged,k];
   end
end
f(ind2notconverged) = NaN;
[sortedf,indf] = sort(f);
fh(1) = figure;
subplot(2,1,1)
plot(sortedf)

ind = indf(1);
sq_seperrorvec =  sqrt( (E(:,1) - xend(:,ind)).^2 + (E(:,2) - yend(:,ind)).^2   + (E(:,3) - zend(:,ind)).^2);
ylabel('L^*')
xlabel('Run Number after sorting')
subplot(2,1,2)
plot(sq_seperror(indf)/(nevnts*2-3-2-1))
print([dirname,filesep,'sortedf_for_allnrand.eps'], '-depsc') 

ylabel('sep error per free location parameter')
xlabel('Run Number after sorting')
disp(['error per invertable parameter: ', num2str(sq_seperror(indf)/(nevnts*2-3-2-1))])
disp(['min coordinate sepp error: ', num2str(min(sq_seperrorvec))])
disp(['max coordinate sepp error: ', num2str(max(sq_seperrorvec))])

fh(2) = figure;
title('sep per parameter')
vec_seperror = sqrt( (E(:,1) - xend(:,k)).^2 + (E(:,2) - yend(:,k)).^2   + (E(:,3) - zend(:,k)).^2 );
hist(vec_seperror,50)





fh(3) = figure;
%[minf,ind]=min(f);


for i = 1:nE
    if i ==1
        plot3(E(i,1), E(i,2), E(i,3),'bo','MarkerSize',10)
        %plot3(E(i,1), E(i,2), E(i,3),'go','markerFaceColor','g')
    elseif i ==2
        plot3(E(i,1), E(i,2), E(i,3),'bo','MarkerSize',10)
        %plot3(E(i,1), E(i,2), E(i,3),'mo','markerFaceColor','m')
    elseif i ==3
        plot3(E(i,1), E(i,2), E(i,3),'bo','MarkerSize',10)
        %plot3(E(i,1), E(i,2), E(i,3),'co','markerFaceColor','c')
    else
        plot3(E(i,1), E(i,2), E(i,3),'bo','MarkerSize',10)
    end
    hold on
    %plot(xend(i,ind(j)),-yend(i,ind(j)), colorlookup{i})
    if ttconst(i,2) == -99999
        plot3(xend(i,ind),yend(i,ind), zend(i,ind), 'ro','MarkerSize',10)  % i.e. no travel time constraint added
    else
        plot3(xend(i,ind),yend(i,ind), zend(i,ind), 'ro','MarkerSize',10,'MarkerFaceColor','r')
    end
    plot3([E(i,1) xend(i,ind)], [E(i,2) yend(i,ind)], [E(i,3) zend(i,ind)],'color',[0.7,0.7,0.7],'Linewidth',2)
end
%set(gca,'xlim',[-70,70], 'ylim',[-70,70], 'fontsize',18)
set(gca,'fontsize',18)
xlabel('x (m)')
ylabel('y (m)') 
zlabel('z (m)')
print([dirname,filesep,'locs_68eq_Calaveras_3_3D.eps'], '-depsc') 
   
fh(4) = figure;
% [minf,ind]=min(f);
% ind = 14
for i = 1:nE
    if i ==1
        plot(E(i,1), E(i,2),'bo','Markersize',10)
        %plot(E(i,1), E(i,2),'go','markerFaceColor','g')
    elseif i ==2
        plot(E(i,1), E(i,2),'bo','Markersize',10)
       % plot(E(i,1), E(i,2),'mo','markerFaceColor','m')
    elseif i ==3
        plot(E(i,1), E(i,2),'bo','Markersize',10)
        %plot(E(i,1), E(i,2),'co','markerFaceColor','c')
    else
        plot(E(i,1), E(i,2),'bo','Markersize',10)
    end
    hold on
    %plot(xend(i,ind(j)),-yend(i,ind(j)), colorlookup{i})
    if ttconst(i,2) == -99999
        plot(xend(i,ind),yend(i,ind), 'ro','Markersize',10)  % i.e. no travel time constraint added
    else
        plot(xend(i,ind),yend(i,ind), 'ro','MarkerFaceColor','r','Markersize',10)
    end
    plot([E(i,1) xend(i,ind)], [E(i,2) yend(i,ind)], 'color',[0.7,0.7,0.7],'linewidth',2)
end
set(gca,'fontsize',18)
xlabel('x (m)')
ylabel('y (m)')
print([dirname,filesep,'locs_68eq_Calaveras_xy.eps'], '-depsc') 
  


fh(5) = figure;
% [minf,ind]=min(f);
% ind = 14
for i = 1:nE
    if i ==1
        plot(E(i,1), E(i,3),'bo','Markersize',10)
        %plot(E(i,1), E(i,3),'go','markerFaceColor','g')
    elseif i ==2
        plot(E(i,1), E(i,3),'bo','Markersize',10)
        %plot(E(i,1), E(i,3),'mo','markerFaceColor','m')
    elseif i ==3
        plot(E(i,1), E(i,3),'bo','Markersize',10)
        %plot(E(i,1), E(i,3),'co','markerFaceColor','c')
    else
        plot(E(i,1), E(i,3),'bo','Markersize',10)
    end
    hold on
    %plot(xend(i,ind(j)),-yend(i,ind(j)), colorlookup{i})
    if ttconst(i,2) == -99999
        plot(xend(i,ind),zend(i,ind), 'ro','Markersize',10)  % i.e. no travel time constraint added
    else
        plot(xend(i,ind),zend(i,ind), 'ro','MarkerFaceColor','r','Markersize',10)
    end
    plot([E(i,1) xend(i,ind)], [E(i,3) zend(i,ind)], 'color',[0.7,0.7,0.7],'linewidth',2)
end
set(gca,'fontsize',18)
xlabel('x (m)')
ylabel('z (m)')
print([dirname,filesep,'locs_68eq_Calaveras_xz.eps'], '-depsc') 


fh(6) = figure;
% [minf,ind]=min(f);
% ind = 14
for i = 1:nE
    if i ==1
        plot(E(i,2), E(i,3),'bo','Markersize',10)
        %plot(E(i,2), E(i,3),'go','markerFaceColor','g')
    elseif i ==2
        plot(E(i,2), E(i,3),'bo','Markersize',10)
        %plot(E(i,2), E(i,3),'mo','markerFaceColor','m')
    elseif i ==3
        plot(E(i,2), E(i,3),'bo','Markersize',10)
        %plot(E(i,2), E(i,3),'co','markerFaceColor','c')
    else
        plot(E(i,2), E(i,3),'bo','Markersize',10)
    end
    hold on
    %plot(xend(i,ind(j)),-yend(i,ind(j)), colorlookup{i})
    if ttconst(i,2) == -99999
        plot(yend(i,ind),zend(i,ind), 'ro','Markersize',10)  % i.e. no travel time constraint added
    else
        plot(yend(i,ind),zend(i,ind), 'ro','MarkerFaceColor','r','Markersize',10)
    end
    plot([E(i,2) yend(i,ind)], [E(i,3) zend(i,ind)], 'color',[0.7,0.7,0.7],'linewidth',2)
end
set(gca,'fontsize',18)
xlabel('y (m)')
ylabel('z (m)')
print([dirname,filesep,'locs_68eq_Calaveras_yz.eps'], '-depsc') 

%%%%% Now let's put them onto the calaveras streaks
hypodd = load([parentdir,'/hypoDD.reloc']);  % loads the complete hypoDD relocations for the 308 events
%hypodd_68 = load('HYPODD/runs/orig_68eqonly/hypoDD.reloc');  % loads the best relocations for the 68 events
hypodd_68 = load([parentdir,'/hypoDD_loc68.txt']);
%cat_loc = load('Calaveras_cat_locations.txt');
%cat_loc_69 = load('Calaveras_cat_locations_69.txt');

[tmp,ind101] = setdiff(hypodd(:,1),hypodd_68(:,1));
hypodd_not68 = hypodd(ind101,:);


fh(7) = figure;
plot(hypodd_not68(:,5), hypodd_not68(:,6),'bo','markerfacecolor','b','markersize',msize)
hold on
for i = 1:nE
    plot(xend(i,ind), yend(i,ind),'ro','markerfacecolor','r','markersize',msize)
end
%     plot(xend(1,ind), yend(1,ind),'go','markerfacecolor','g','markersize',msize)
%     plot(xend(2,ind), yend(2,ind),'yo','markerfacecolor','y','markersize',msize)
%     plot(xend(3,ind), yend(3,ind),'co','markerfacecolor','c','markersize',msize)
    

xlabel('x (m)','fontsize',fsize)
ylabel('y (m)','fontsize',fsize)
set(gca,'xlim',xlimits, 'ylim',ylimits,'units','centimeters')
set(gca,'position',axes_posy)
set(gca,'xtick', xtickspots, 'ytick',ytickspots)
set(gca,'fontsize',fsize)
print([dirname,filesep,'example_68eq_calaveras_loc_xy.eps'], '-depsc') 

% First plot the catalogue locations
fh(8) = figure;
plot(hypodd_not68(:,5), hypodd_not68(:,7),'bo','markerfacecolor','b','markersize',msize)
hold on
for i = 1:nE
    plot(xend(i,ind), zend(i,ind),'ro','markerfacecolor','r','markersize',msize)
end
%     plot(xend(1,ind), zend(1,ind),'go','markerfacecolor','g','markersize',msize)
%     plot(xend(2,ind), zend(2,ind),'yo','markerfacecolor','y','markersize',msize)
%     plot(xend(3,ind), zend(3,ind),'co','markerfacecolor','c','markersize',msize)

xlabel('x (m)','fontsize',fsize)
ylabel('z (m)','fontsize',fsize)
set(gca,'xlim',xlimits, 'ylim',zlimits,'units','centimeters')
set(gca,'position',axes_posy)
set(gca,'xtick', xtickspots, 'ytick',ztickspots)
set(gca,'fontsize',fsize)
print([dirname,filesep,'example_68eq_calaveras_loc_xz.eps'], '-depsc') 

% First plot the catalogue locations
fh(9) = figure;
plot(hypodd_not68(:,6), hypodd_not68(:,7),'bo','markerfacecolor','b','markersize',msize)
hold on
for i = 1:nE
    plot(yend(i,ind), zend(i,ind),'ro','markerfacecolor','r','markersize',msize)
end
%     plot(yend(1,ind), zend(1,ind),'go','markerfacecolor','g','markersize',msize)
%     plot(yend(2,ind), zend(2,ind),'yo','markerfacecolor','y','markersize',msize)
%     plot(yend(3,ind), zend(3,ind),'co','markerfacecolor','c','markersize',msize)

xlabel('y (m)','fontsize',fsize)
ylabel('z (m)','fontsize',fsize)
set(gca,'xlim',xlimits, 'ylim',zlimits,'units','centimeters')
set(gca,'position',axes_posy)
set(gca,'xtick', ytickspots, 'ytick',ztickspots)
set(gca,'fontsize',fsize)
print([dirname,filesep,'example_68eq_calaveras_loc_yz.eps'], '-depsc') 


% Do some simple statistics
hypoDD_origx = E(:,1);
best_CWIx = xend(:,ind);
hypoDD_origy = E(:,2);
best_CWIy = yend(:,ind);
hypoDD_origz = E(:,3);
best_CWIz = zend(:,ind);

deltax = hypoDD_origx-best_CWIx;
deltay = hypoDD_origy-best_CWIy;
deltaz = hypoDD_origz-best_CWIz;
delta_coord = [deltax;deltay;deltaz];

CWIstats.minval = min(abs(delta_coord));
CWIstats.maxval = max(abs(delta_coord));
CWIstats.meanval = mean(abs(delta_coord));
CWIstats.stdval = std(abs(delta_coord));
CWIstats.nE = nE; 

% Now let us find the best fitting plane
points_in = [best_CWIx,best_CWIy,best_CWIz];
[x0, normal, resid, norm_resid] = lsplane(points_in);
CWIstats.plane.normal = normal;
CWIstats.plane.resid = resid;
CWIstats.plane.norm_resid = norm_resid;
CWIstats.plane.rms = norm_resid/sqrt(length(resid));


Locs.relocx = best_CWIx; 
Locs.relocy = best_CWIy; 
Locs.relocz = best_CWIz;
Locs.not68x = hypodd_not68(:,5);
Locs.not68y = hypodd_not68(:,6);
Locs.not68z = hypodd_not68(:,7);
