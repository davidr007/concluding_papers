%% This script sets up example 1 and runs it

% Load, subsample and save the CWI_data

tic
%Calaveras EXAMPLE - 10eq 
fdom = 2.5;
vs = 2728.3;  % vs = 3300;

path2eventsoi ='c:/datafiles/workstuff/sandpit/davidr/thesis_version2/diags/eq_location_optimisation/events_oi.txt';
copyfile(path2eventsoi)
load events_oi.txt
[nE2,mE2] = size(events_oi);


hypoDD_loc68= load('C:\datafiles\workstuff\sandpit\davidr\concluding_papers\paper2\diags\extra_hypoDD3/hypoDD.reloc');
E = [];
count = 1;
for i = 1:length(events_oi);
    ind = find(events_oi(i)==hypoDD_loc68(:,1));
    if ~isempty(ind)
        E(count,:) = hypoDD_loc68(ind,5:7);
        sigmaE(count,:) = hypoDD_loc68(ind,8:10);
        count = count+1;
    end 
end

% events_oi = [events_oi(1:2); events_oi(5); events_oi(3:4); events_oi(6:end)];
[nE,tmp] = size(E);
nrand = 1;
npar = nE2*3;%-3-2-1;

%USING REAL EVENTS CALAVERAS
path2CWIstat_data ='c:/datafiles/workstuff/sandpit/davidr/thesis_version2/diags/revised_2dacoustics_1to5Hz/Calaveras_CWI_stat.mat';
%events_oi = [101361,102139,102776,105327];
[CWI_stat_small,CWI_stat_large] = subsample_CWI_stats(path2CWIstat_data,events_oi(:,1));
CWI_stat = CWI_stat_large;
CWI_stat = [CWI_stat,vs*ones(size(CWI_stat_large,1),1)];
% Now let's pull out all the contsraints with "large" mun
ind = find(CWI_stat(:,6)>0.1);
CWI_stat(ind,3:9) = [-99999*ones(length(ind),7)];
save CWI_stat.txt CWI_stat -ascii


xbounds = [-200*ones(nE2,1), 200*ones(nE2,1)]; 
%xbounds(2,:) = [0,1000];
ybounds = [-200*ones(nE2,1), 200*ones(nE2,1)]; 
%ybounds(3, :) = [0, 100];     % i.e. we force y3 to be positive
zbounds = [-200*ones(nE2,1), 200*ones(nE2,1)];


% Travel time constraints
ttconst = zeros(nE2,7);
E2 = zeros(nE2,4);
%ttconst_pointer = ones(nE2,1);
%ttconst_pointer(1:8) = [1, 1,1,1,1,1,1,1];
%ttconst_pointer(1:nE) = ones(nE,1);
%ttconst_pointer = ones(nE,1);
Hypodd_cart = E;%events_oi_tmp2(i,5:7);
Hypodd_loc = Hypodd_cart;
%Hypodd_loc = do_flips4inversion(Hypodd_loc_tmp);
for i = 1:nE2
    ind = find(hypoDD_loc68(:,1)== events_oi(i));
    if ~isempty(ind)  % add the event location contsraint from hypoDD
        % Now build the travel time constraint matrix
        ttconst(i,:) = [ events_oi(i),Hypodd_loc(ind,:),sigmaE(ind,:)];%19.5,19.5,15];
        E2(i,:) = [Hypodd_loc(ind,:),1];
    elseif isempty(ind)
        ttconst(i,:) = [ events_oi(i),-99999*ones(1,6)];
        E2(i,:) = [xbounds(i,1) + (xbounds(i,2)- xbounds(i,1)).*rand(1,1), ...
            ybounds(i,1) + (ybounds(i,2)- ybounds(i,1)).*rand(1,1), ...
            zbounds(i,1) + (zbounds(i,2)- zbounds(i,1)).*rand(1,1),0];
    else
        error('ERROR 121: You should not be here')
    end
end
ttconst(ttconst==0) = -99999;
save ttconst.txt ttconst -ascii



%zbounds(4, :) = [0, 100];     % i.e. we force y4 to be positive
randjumpx = zeros(nE2,nrand);
randjumpy = zeros(nE2,nrand);
randjumpz = zeros(nE2,nrand);
for i = 1:nrand
   randjumpx(:,i) = xbounds(i,1) + (xbounds(i,2)- xbounds(i,1)).*rand(1,nE2);
   randjumpy(:,i) = ybounds(i,1) + (ybounds(i,2)- ybounds(i,1)).*rand(1,nE2);
   randjumpz(:,i) = zbounds(i,1) + (zbounds(i,2)- zbounds(i,1)).*rand(1,nE2);
end

xstart = repmat(E2(:,1),1,nrand)+randjumpx;
ystart = repmat(E2(:,2),1,nrand)+randjumpy;
zstart = repmat(E2(:,3),1,nrand)+randjumpz;


%xstart = E(:,1) + 10;
%ystart = E(:,2) -10;
%zstart = E(:,3) -5;

% create a matricy of the xvec start locations
xvecstart = zeros(npar,nrand);
for i = 1:nrand
    xvecstart_tmp = [];
    %xvecstart_tmp = [xstart(2,i), xstart(3,i), ystart(3,i)];
    for j= 1    :nE2
        xvecstart_tmp = [xvecstart_tmp, xstart(j,i), ystart(j,i),zstart(j,i)];
    end
    xvecstart(:,i) = xvecstart_tmp';
end
save xvecstart.txt xvecstart -ascii


if strcmp(filesep,'/') % only do this if on LINUX)
    currentdir = pwd;  %store the current directory
    cd ~/storage/sandpit/davidr/events/fortran/
    %eval('!make clean')    % clean the slate in teh fortran events area
    eval('!make optimise_lstar_3D_cartesian') % re-make the program to ensure you get the latest version
    cd(currentdir) % cd back to the working area
    eval('!cp ~/storage/sandpit/davidr/events/fortran/optimise_lstar_3D_cartesian ./') % get the binary
    eval(['! echo "',num2str(nrand),' ',num2str(nE2),'" | ./optimise_lstar_3D_cartesian']) % run the binary
    

else
    error(['MESSAGE 122: you must manually move the files to LINUX and run the optimiser ...', ...
        'THEN run the remainder of this script manually'])
end


xendvec_tmp = load('xend.txt');
xendvec = zeros(nrand, npar);
count = 3;
for i = 1:nrand
    xendvec(i,:) = xendvec_tmp(count:count+npar-1);
    count = count +npar;
end
ftmp = load('fend.txt');
f = ftmp(:,1);


% Unwrap the solution into separate matrices for x and y coordinates
xend_tmp = zeros(size(xstart));
yend_tmp = zeros(size(ystart));
zend_tmp = zeros(size(zstart));
%xend_tmp(2,:) = xendvec(:,1)';
%xend_tmp(3,:) = xendvec(:,2)';
%yend_tmp(3,:) = xendvec(:,3)';
count = 1;
for j = 1:nE2
    xend_tmp(j,:) =  xendvec(:,count)';
    count = count+1;
    yend_tmp(j,:) =  xendvec(:,count)';
    count = count+1;
    zend_tmp(j,:) = xendvec(:,count)';
    count = count+1;
end

% Note that we are working with the cartesian coordinates so there is no need
% to do any converting
xend = xend_tmp;
yend = yend_tmp;
zend = zend_tmp;
%xend = zeros(size(xstart));
%yend = zeros(size(ystart));
%zend = zeros(size(zstart));
%for i = 1:nrand
%    [cart_out_tmp] = convert_eventcoordinates(E(1,:),E(2,:),E(3,:), 'loc2cart', [xend_tmp(:,i),yend_tmp(:,i), zend_tmp(:,i)])
%    xend(:,i) = cart_out_tmp(:,1);
%    yend(:,i) = cart_out_tmp(:,2);
%    zend(:,i) = cart_out_tmp(:,3);
%end


toc

disp(['runtime for original - i.e. not staggered'])


opttype = 'fortran';

E = E2;

save(['randsearch_',opttype,'.mat'],'xend','yend','zend','f','E','xstart','ystart','ttconst','events_oi')

















% save CWI_stat.mat CWI_stat events_oi
%
% n = length(events_oi(:,1));
% xstart = 100*rand(n*3-5,1);
% lb = [0 -100 0 0];  % i.e. x2>0, x3=free, y3>0, z3>0
% ub = 100*ones(size(xstart'));
% count = 5
% for i = 4:n
%     lb(count) = -100;  % i.e. x4,x5,x6,.. = free
%     lb(count+1) = -100; % i.e. y4,y5,y6,.. = free
%     lb(count+2) = 0; % i.e. z4,z5,z6,.. >0
%     count = count+3;
% end
%
% OPTIONS = optimset('MaxFunEvals',2000,'gradobj','on')%,'HessUpdate','dfp','TolX',10^(-40));
%
%
% [x,fval,exitflag,output] = fmincon(@Lstar,xstart,[],[],[],[],lb,ub);
% %[f grad] = Lstar(xstart);
%
%
% xvec = [0; x(1);x(2:3:end)];
% yvec = [0; 0; x(3:3:end)];
% zvec = [0; 0;x(4:3:end)];


