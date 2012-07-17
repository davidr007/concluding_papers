function [inversion_results,fh] = plot_outcomes(dirname)
currentdir = pwd;
cd(dirname)

opttype = 'fortran'; % 'fortran' or 'matlab'
load(['randsearch_',opttype,'.mat'])
[sorted_x2, ind]=sort(xend(2,:));

[nevnts, nrand] = size(xend);

inversion_results.otherstuff.E = E;

% Have a quick look at the solution space
%fh.allsolutions = figure;
fh = [];
colorlookup = {'b.','r.','g.','c.','m.','y.','k.','bo','ro','go'};
colorlookup2 = {'b^','r^','g^','c^','m^','y^','k^','bd','rd','gd'};

%% Plot the joiner lines between actual location and optimisation solns
%for i = 1:nevnts
%    for j = 1:nrand
%        plot([xend(i,j) E(i,1)], [yend(i,j) E(i,2)],'Color','k')
%        hold on
%    end
%end
% Plot the optimisation solutions
%for i = 1:nevnts
%    plot(xend(i,:), yend(i,:),'ro','markersize',7)
%    hold on
%end
% Plot the actual solutions
%for i = 1:nevnts
%    plot(E(i,1),E(i,2), 'bo', 'markersize',7)
%    hold on
%end
%plot([ 0 0], get(gca,'ylim'),'k')
%plot(get(gca,'xlim'),[0 0],'k')

[nevnts, nrand] = size(xend);


% Do analysis of all solutions
ExMat = repmat(E(2:end,1), 1,nrand);  % Note we are omiting the 0 values
EyMat = repmat(E(3:end,2), 1,nrand);
EzMat = repmat(E(3:end,3), 1,nrand);
deltaxMat = ExMat - xend(2:end,:); 
deltayMat = EyMat - yend(3:end,:);
deltazMat = EzMat - zend(3:end,:);
deltax = deltaxMat(:);
deltay = deltayMat(:);
deltaz = deltazMat(:);
delta_coord = [deltax; deltay; deltaz];
inversion_results.allsol.nrand = nrand;
inversion_results.allsol.fminval = f;
inversion_results.allsol.minval = min(abs(delta_coord));
inversion_results.allsol.maxval = max(abs(delta_coord));
inversion_results.allsol.meanval = mean(abs(delta_coord));
inversion_results.allsol.stdval = std(abs(delta_coord));

% Do analysis for the best solution. 
clear deltax deltay deltaz delta_coord
[fminval,ind_fmin] = min(f);
deltax = E(2:end,1) - xend(2:end,ind_fmin);
deltay = E(3:end,2) - yend(3:end,ind_fmin);
deltaz = E(3:end,3) - zend(3:end,ind_fmin);
delta_coord = [deltax; deltay; deltaz];
inversion_results.bestsol.fminval = fminval;
inversion_results.bestsol.minval = min(abs(delta_coord));
inversion_results.bestsol.maxval = max(abs(delta_coord));
inversion_results.bestsol.meanval = mean(abs(delta_coord));
inversion_results.bestsol.stdval = std(abs(delta_coord));

% Draw the plot for the best solution 
%fh.bestsolution = figure;
% Plot the joiner lines between actual location and optimisation solns
%for i = 1:nevnts
%        plot([xend(i,ind_fmin) E(i,1)], [yend(i,ind_fmin) E(i,2)],'Color','k')
%        hold on
%end
%plot([ 0 0], get(gca,'ylim'),'k')
%plot(get(gca,'xlim'),[0 0],'k')




cd(currentdir)
% 



