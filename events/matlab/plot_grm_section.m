function [grm_lims,hc] = plot_grm_section(fname,xcut,zcut,deltas,caxis_lims)

% function to plot a grm cross section (works with 2D only). 
% INPUTS:
% fname         [string] full path filename to media of interest
% xcut          [vector 2x1]  i.e. [lower_xlim, upper_xlim] for cut
% zcut          [vector 2x1]  i.e. [lower_zlim, upper_zlim] for cut
% deltas        [vector 3x1] [deltax, deltay, deltaz] in meters
% caxis_lims    [vector 2x1] lower and upper limits for caxis re-scale of
%               images. Note that scale is automatic if caxis_lims = [];

%grm_full = load(fname);

%grm = grm_full(zcut(1):zcut(2),xcut(1):xcut(2));

grm = dlmread(fname, '', [zcut(1),xcut(1),zcut(2),xcut(2)]);
[n m] = size(grm);

x = 0:deltas(1):(m-1)*deltas(1);
z = 0:deltas(3):(n-1)*deltas(3);

figure
h = pcolor(x/1000,z/1000,grm/1000);
set(gca,'fontsize',26)
shading interp
if ~isempty(caxis_lims)
    caxis([caxis_lims(1), caxis_lims(2)]);
end
hc = colorbar
set(hc, 'fontsize',26)
disp('   ')

grm_lims(1) = min(grm(:));
grm_lims(2) = max(grm(:));







