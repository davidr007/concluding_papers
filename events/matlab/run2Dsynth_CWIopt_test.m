function [] = run2Dsynth_CWIopt_test(nE,opttype,npar,nrand,E)
% Runs the optimisation for each of the 2D synthetic data cases in 
% */thesis_version2/diags/eq_location_optimisation/synth2Dmulti
% (see run_script for the wrapper which calls it)



% The purpose of this script is to analyse the results given a large
% number of random starts. 



% loads: tw, xout_mean_twp75,  xout_sigma_twp75

xbounds = [-100*ones(nE,1), 100*ones(nE,1)]; 
xbounds(2,:) = [0,100];
ybounds = [-100*ones(nE,1), 100*ones(nE,1)]; 
ybounds(3, :) = [0, 100];     % i.e. we force y3 to be positive
xstart = zeros(nE,nrand);
ystart = zeros(nE,nrand);
for i = 1:nE
   xstart(i,:) = xbounds(i,1) + (xbounds(i,2)- xbounds(i,1)).*rand(1,nrand);
   ystart(i,:) = xbounds(i,1) + (xbounds(i,2)- xbounds(i,1)).*rand(1,nrand);
end

switch opttype
    case 'fortran'
        % create a matricy of the xvec start locations 
        xvecstart = zeros(npar,nrand);
        for i = 1:nrand
            xvecstart_tmp = [xstart(2,i)];
            for j= 3    :nE
                xvecstart_tmp = [xvecstart_tmp, xstart(j,i), ystart(j,i)];
            end
            xvecstart(:,i) = xvecstart_tmp';
        end
        save xvecstart.txt xvecstart -ascii
        
        currentdir = pwd;  %store the current directory
        cd ~/storage/sandpit/davidr/events/fortran/
        %eval('!make clean')    % clean the slate in teh fortran events area
        eval('!make optimise_lstar_2D') % re-make the program to ensure you get the latest version
        cd(currentdir) % cd back to the working area
        eval('!cp ~/storage/sandpit/davidr/events/fortran/optimise_lstar_2D ./') % get the binary  
        eval(['! echo "',num2str(nrand),' ',num2str(nE),'" | ./optimise_lstar_2D']) % run the binary
        
        xendvec_tmp = load('xend.txt');
        xendvec = zeros(nrand, npar);
        count = 3;
        for i = 1:nrand
            xendvec(i,:) = xendvec_tmp(count:count+npar-1);
            count = count +npar;
        end
        ftmp = load('fend.txt');
        f = ftmp(:,1);
    case 'matlab'
        load('/export/storage/davidr/sandpit/davidr/thesis_version2/diags/revised_2dacoustics/foldederrorbar_fit3_normtype1.mat')
        nvec = nE*2-2-1;
        lb = -Inf*ones(1,nvec);
        lb(1) = 0;  % event x2 constrained onto positive x-acis
        lb(3) = 0; %event 3 y-value constrained to be positive.
        xendvec = zeros(nrand, npar);


        for i = 1:nrand
            disp(['Doing i = ',num2str(i),' of ', num2str(nrand), ' - setup']);
            xvecstart = [xstart(2,i)];
            for j= 3    :nE
                xvecstart = [xvecstart, xstart(j,i), ystart(j,i)];
            end
            disp(['Doing i = ',num2str(i),' of ', num2str(nrand), ' - optimistation']);
            OPTIONS = optimset('MaxFunEvals',50,'gradobj','on');%,'TolFun', 0.01);
            %[xendvec(i,:),f(i),exitflag,output] = fmincon(@Lstar_2D,xvecstart,[],[],[],[],lb,[], [], OPTIONS);
            [xendvec(i,:),f(i),exitflag,output] = fminunc(@Lstar_2D,xvecstart,OPTIONS);
        end


end


% Unwrap the solution into separate matrices for x and y coordinates
xend = zeros(size(xstart));
yend = zeros(size(ystart));
xend(2,:) = xendvec(:,1)';
count = 2;
for j = 3:nE
    xend(j,:) =  xendvec(:,count)';
    count = count+1;
    yend(j,:) =  xendvec(:,count)';
    count = count+1;
end


% Now as a final step we should flip yend to make sure yend(3,j) has the
% same sign as E(3,2)
for j = 1:nrand
    if (E(3,2)<0 & yend(3,j)>0) | (E(3,2)>0 & yend(3,j)<0)
        yend(:,j) = -yend(:,j);
    end
end
% we also need to flip xend to ensure that xend has the correct sign
for j = 1:nrand
    if (E(2,1)<0 & xend(2,j)>0) | (E(2,1)>0 & xend(2,j)<0)
        xend(:,j) = -xend(:,j);
    end
end


toc

disp(['runtime for original - i.e. not staggered'])




save(['randsearch_',opttype,'.mat'],'xend','yend','f','E','xstart','ystart')



