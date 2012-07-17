function [aho] = double_xticks(ahi,fhi,newticklocs, newlabels)
% This function creates a second axes (aho) underneath the original one
% (ahi) and uses it to create a second set of xticklabels
%
% INPUTS:
% ahi           [scalar] axis handle in
% fhi           [scalar] fhi = figure handle in
% newlabels     [double] 
%
% OUTPUTS: 
% aho           [scalar] axis handle out

if nargin ==0
   x = 0:0.1:25;
   y = exp(x);
   fhi = figure;
   plot(x,y)
   ahi = gca;
   set(ahi,'xtick',[0 5 10 15 20 25],'fontsize',10)
   newlabels = {'(n0)', '(n5)', '(n10)','(n15)', '(n20)','(n25)'};
end

ahi_position = get(ahi,'position');
ahi_fontsize = get(ahi,'fontsize');
ahi_xlim = get(ahi,'xlim')
ahi_xtick = get(ahi,'xtick')

voffset = 20; % as a percentage of the total height
aho_position(1)  = ahi_position(1);
aho_position(2)  = ahi_position(2)-voffset/100*ahi_position(4);  % move the new axes down slightly 
aho_position(3)  = ahi_position(3);
aho_position(4)  = 0.001;%voffset/(2*100)*ahi_position(4);  % make the height of the new axes very small

aho = axes('units','centimeters','position', aho_position);
set(aho,'xlim',ahi_xlim)
set(aho,'xtick', newticklocs)
set(aho,'xticklabel', newlabels)
set(aho,'ycolor','w')
%set(aho,'TickLength',[0 0])
set(aho,'fontsize',ahi_fontsize)
1+1
