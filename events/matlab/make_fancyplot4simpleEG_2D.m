function [fig_handle,a1,a2] = make_fancyplot4simpleEG_2D(grid_x2, grid_y2, PDF_2D,profile, caxis_lims,user_colormap,fancyplot_zlims)
% Does the fancy plot for a 2D PDF of earthquake location. 
% Example Useage: 
%    make_fancyplot4simpleEG_2D(grid_x2_1, grid_y2_1, PDF_2D_EG1,profile_EG1, caxis_lims,tmp2)

fig_handle = figure;
a1 = subplot(3,1,[1:2]);
ind = find(grid_x2>=0);
sh_1 = surface(grid_x2(ind),grid_y2,PDF_2D(:,ind))
hold on
zdata_tmp = get(sh_1,'zdata');
l1 = plot3(zeros(size(profile(:,1)))-0.005,profile(:,1),zdata_tmp(:,1)+0.001,'k','linewidth',1)
%meshc(grid_x2_1,grid_y2_1,PDF_2D_EG1)
%xlabel('x')
%ylabel('y')
shading interp
caxis(caxis_lims)
colormap(user_colormap)
%axis equal
set(a1,'view', [-46 18])
set(a1,'zlim', fancyplot_zlims)
set(gca,'xlim', [min(grid_x2), max(grid_x2)], 'ylim', [min(grid_y2), max(grid_y2)])
grid on

fsize = 16;
h1 = ylabel('$\widetilde{y}$','interpreter','latex','fontsize',fsize)
yposy = get(h1,'position')
set(h1,'position', [yposy(1)+0.1,yposy(2)+0.5,yposy(3)])
h2 = xlabel('$\widetilde{x}$','Interpreter','LaTex','fontsize',fsize)
xposy = get(h2,'position')
set(h2,'position', [xposy(1)+0.1+0.4,xposy(2)+0.1,xposy(3)])
%zlabel('$P(E_2|E_1,\widetilde{\delta}_{CWIN})$','Interpreter','LaTex','fontsize',14)
set(gca,'fontsize',fsize)

a2 = subplot(3,1,3);
pcolor(grid_x2,grid_y2,PDF_2D)
shading interp
caxis(caxis_lims)
axis equal
hold on
plot([-1.2 1.2], [1.2 -1.2],'k--','linewidth',2)

set(gca,'zlim',[-0.00001 0.000001])
%set(a2,'view',[30 30])
axis tight
set(gca,'ztick',[],'ztickLabel',' ')
%xlabel('x')
%ylabel('y')
set(a2,'view', [-46 18])
set(a2,'position', [-0.1    0.05    1.5*0.7750    1.5*0.2157])

% h1 = ylabel('$\widetilde{y}$','fontsize',fsize)
% yposy = get(h1,'position')
% set(h1,'position', [yposy(1)+0.3,yposy(2)+0.5,yposy(3)],'Interpreter','LaTex')
% h2 = xlabel('$\widetilde{x}$','fontsize',fsize)
% xposy = get(h2,'position')
% set(h2,'position', [xposy(1)+0.7,xposy(2)+0.5,xposy(3)],'Interpreter','LaTex')
set(gca,'fontsize',fsize)
% 
