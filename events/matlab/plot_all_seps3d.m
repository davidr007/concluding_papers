function [] = plot_all_seps3d(astruc,expected,col,plotstruct,startval)

% plot_all_seps3d plots all of the separation values for a 3D simulation
% INPUTS
% astruc
% expected  [scalar] the expected separation. 
% col       [scalar] column to use in the matrices (i.e. in sepStruc)
% plotstruct [structure] plotting preferences
%                   plotstruct.xlimits
%                   plotstruct.xticks
%                   plotstruct.ylimits
%                   plotstruct.yticks
% startval  [scalar] the value to start at when computing the mean

figure
colors = [0 0 0;1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 0.5 0.5 0.5; 0.992 0.918 0.796]; 
colors = 0.5*ones(9,3); 

expossy =  plotstruct.xlimits(2) - 5/100*( plotstruct.xlimits(2)- plotstruct.xlimits(1));

col100ch1 = [];
col100ch2 = [];
col100ch3 = [];
col200ch1 = [];
col200ch2 = [];
col200ch3 = [];
col300ch1 = [];
col300ch2 = [];
col300ch3 = [];
col400ch1 = [];
col400ch2 = [];
col400ch3 = [];

for j = 1:9
    % Draw the first layer
    subplot(4,3,1)
    sep_tmp2 = astruc.(['station10',num2str(j)]).chan1;
    sep_tmp2(isnan(sep_tmp2)) = 0;
    plot(sep_tmp2(:,1),sep_tmp2(:,2),'color',colors(j,:))
    hold on
    if j ==1
        plot(plotstruct.xlimits, expected,'k--','linewidth',2)
        title('x channel')
        ylabel('<r> (surface)')
        set(gca,'xlim',plotstruct.xlimits,'ylim',plotstruct.ylimits, ...
            'xtick',plotstruct.xticks,'ytick',plotstruct.yticks)
    end
    col100ch1 = [col100ch1,sep_tmp2(startval:end,2)'];
    if j==9
        [mean100ch1, std100ch1] = compute_stats(col100ch1);
        plot(expossy,mean100ch1,'ko')
        errorbar(expossy,mean100ch1,std100ch1,'k')
    end
    subplot(4,3,2)
    sep_tmp2 = astruc.(['station10',num2str(j)]).chan2;
    sep_tmp2(isnan(sep_tmp2)) = 0;
    plot(sep_tmp2(:,1),sep_tmp2(:,2),'color',colors(j,:))
    hold on
    if j ==1
        plot(plotstruct.xlimits, expected,'k--','linewidth',2)
        title('y channel')
        set(gca,'xlim',plotstruct.xlimits,'ylim',plotstruct.ylimits, ...
            'xtick',plotstruct.xticks,'ytick',plotstruct.yticks)
    end
    col100ch2 = [col100ch2,sep_tmp2(startval:end,2)'];
    if j==9
        [mean100ch2, std100ch2] = compute_stats(col100ch2);
        plot(expossy,mean100ch2,'ko')
        errorbar(expossy,mean100ch2,std100ch2,'k')
    end
    subplot(4,3,3)
    sep_tmp2 = astruc.(['station10',num2str(j)]).chan3;
    sep_tmp2(isnan(sep_tmp2)) = 0;
    plot(sep_tmp2(:,1),sep_tmp2(:,2),'color',colors(j,:))
    hold on
    if j ==1
        plot(plotstruct.xlimits, expected,'k--','linewidth',2)
        title('z channel')
        set(gca,'xlim',plotstruct.xlimits,'ylim',plotstruct.ylimits, ...
            'xtick',plotstruct.xticks,'ytick',plotstruct.yticks)
    end
    col100ch3 = [col100ch3,sep_tmp2(startval:end,2)'];
    if j==9
        [mean100ch3, std100ch3] = compute_stats(col100ch3);
        plot(expossy,mean100ch3,'ko')
        errorbar(expossy,mean100ch3,std100ch3,'k')
    end

    % Draw the second layer
    subplot(4,3,4)
    sep_tmp2 = astruc.(['station20',num2str(j)]).chan1;
    sep_tmp2(isnan(sep_tmp2)) = 0;
    plot(sep_tmp2(:,1),sep_tmp2(:,2),'color',colors(j,:))
    hold on
    if j ==1
        plot(plotstruct.xlimits, expected,'k--','linewidth',2)
        set(gca,'xlim',plotstruct.xlimits,'ylim',plotstruct.ylimits, ...
            'xtick',plotstruct.xticks,'ytick',plotstruct.yticks)
        ylabel('<r> (layer 2)')
    end
    col200ch1 = [col200ch1,sep_tmp2(startval:end,2)'];
    if j==9
        [mean200ch1, std200ch1] = compute_stats(col200ch1);
        plot(expossy,mean200ch1,'ko')
        errorbar(expossy,mean200ch1,std200ch1,'k')
    end
    subplot(4,3,5)
    sep_tmp2 = astruc.(['station20',num2str(j)]).chan2;
    sep_tmp2(isnan(sep_tmp2)) = 0;
    plot(sep_tmp2(:,1),sep_tmp2(:,2),'color',colors(j,:))
    hold on
    if j ==1
        plot(plotstruct.xlimits, expected,'k--','linewidth',2)
        set(gca,'xlim',plotstruct.xlimits,'ylim',plotstruct.ylimits, ...
            'xtick',plotstruct.xticks,'ytick',plotstruct.yticks)
    end
    col200ch2 = [col200ch2,sep_tmp2(startval:end,2)'];
    if j==9
        [mean200ch2, std200ch2] = compute_stats(col200ch2);
        plot(expossy,mean200ch2,'ko')
        errorbar(expossy,mean200ch2,std200ch2,'k')
    end
    subplot(4,3,6)
    sep_tmp2 = astruc.(['station20',num2str(j)]).chan3;
    sep_tmp2(isnan(sep_tmp2)) = 0;
    plot(sep_tmp2(:,1),sep_tmp2(:,2),'color',colors(j,:))
    hold on
    if j ==1
        plot(plotstruct.xlimits, expected,'k--','linewidth',2)
        set(gca,'xlim',plotstruct.xlimits,'ylim',plotstruct.ylimits, ...
            'xtick',plotstruct.xticks,'ytick',plotstruct.yticks)
    end
    col200ch3 = [col200ch3,sep_tmp2(startval:end,2)'];
    if j==9
        [mean200ch3, std200ch3] = compute_stats(col200ch3);
        plot(expossy,mean200ch3,'ko')
        errorbar(expossy,mean200ch3,std200ch3,'k')
    end
    % Draw the third layer
    if j~=5
        subplot(4,3,7)
        sep_tmp2 = astruc.(['station30',num2str(j)]).chan1;
        sep_tmp2(isnan(sep_tmp2)) = 0;
        plot(sep_tmp2(:,1),sep_tmp2(:,2),'color',colors(j,:))
        hold on
        if j ==1
            plot(plotstruct.xlimits, expected,'k--','linewidth',2)
            set(gca,'xlim',plotstruct.xlimits,'ylim',plotstruct.ylimits, ...
            'xtick',plotstruct.xticks,'ytick',plotstruct.yticks)
            ylabel('<r> (layer 3)')
        end
        col300ch1 = [col300ch1,sep_tmp2(startval:end,2)'];
        if j==9
            [mean300ch1, std300ch1] = compute_stats(col300ch1);
            plot(expossy,mean300ch1,'ko')
            errorbar(expossy,mean300ch1,std300ch1,'k')
        end
        subplot(4,3,8)
        sep_tmp2 = astruc.(['station30',num2str(j)]).chan2;
        sep_tmp2(isnan(sep_tmp2)) = 0;
        plot(sep_tmp2(:,1),sep_tmp2(:,2),'color',colors(j,:))
        hold on
        if j ==1
            plot(plotstruct.xlimits, expected,'k--','linewidth',2)
            set(gca,'xlim',plotstruct.xlimits,'ylim',plotstruct.ylimits, ...
            'xtick',plotstruct.xticks,'ytick',plotstruct.yticks)
        end
        col300ch2 = [col300ch2,sep_tmp2(startval:end,2)'];
        if j==9
            [mean300ch2, std300ch2] = compute_stats(col300ch2);
            plot(expossy,mean300ch2,'ko')
            errorbar(expossy,mean300ch2,std300ch2,'k')
        end
        subplot(4,3,9)
        sep_tmp2 = astruc.(['station30',num2str(j)]).chan3;
        sep_tmp2(isnan(sep_tmp2)) = 0;
        plot(sep_tmp2(:,1),sep_tmp2(:,2),'color',colors(j,:))
        hold on
        if j ==1
            plot(plotstruct.xlimits, expected,'k--','linewidth',2)
            set(gca,'xlim',plotstruct.xlimits,'ylim',plotstruct.ylimits, ...
            'xtick',plotstruct.xticks,'ytick',plotstruct.yticks)
        end
        col300ch3 = [col300ch3,sep_tmp2(startval:end,2)'];
        if j==9
            [mean300ch3, std300ch3] = compute_stats(col300ch3);
            plot(expossy,mean300ch3,'ko')
            errorbar(expossy,mean300ch3,std300ch3,'k')
        end
    end
   
    
    % Draw the fourth layer
    subplot(4,3,10)
    sep_tmp2 = astruc.(['station20',num2str(j)]).chan1;
    sep_tmp2(isnan(sep_tmp2)) = 0;
    plot(sep_tmp2(:,1),sep_tmp2(:,2),'color',colors(j,:))
    hold on
    if j ==1
        plot(plotstruct.xlimits, expected,'k--','linewidth',2)
        set(gca,'xlim',plotstruct.xlimits,'ylim',plotstruct.ylimits, ...
            'xtick',plotstruct.xticks,'ytick',plotstruct.yticks)
        ylabel('<r> (layer 4)')
    end
    col400ch1 = [col400ch1,sep_tmp2(startval:end,2)'];
    if j==9
        [mean400ch1, std400ch1] = compute_stats(col400ch1);
        plot(expossy,mean400ch1,'ko')
        errorbar(expossy,mean400ch1,std400ch1,'k')
    end
    subplot(4,3,11)
    sep_tmp2 = astruc.(['station20',num2str(j)]).chan2;
    sep_tmp2(isnan(sep_tmp2)) = 0;
    plot(sep_tmp2(:,1),sep_tmp2(:,2),'color',colors(j,:))
    hold on
    if j ==1
        plot(plotstruct.xlimits, expected,'k--','linewidth',2)
        set(gca,'xlim',plotstruct.xlimits,'ylim',plotstruct.ylimits, ...
            'xtick',plotstruct.xticks,'ytick',plotstruct.yticks)
    end
    col400ch2 = [col400ch2,sep_tmp2(startval:end,2)'];
    if j==9
        [mean400ch2, std400ch2] = compute_stats(col400ch2);
        plot(expossy,mean400ch2,'ko')
        errorbar(expossy,mean400ch2,std400ch2,'k')
    end
    subplot(4,3,12)
    sep_tmp2 = astruc.(['station20',num2str(j)]).chan3;
    sep_tmp2(isnan(sep_tmp2)) = 0;
    plot(sep_tmp2(:,1),sep_tmp2(:,2),'color',colors(j,:))
    hold on
    if j ==1
        plot(plotstruct.xlimits, expected,'k--','linewidth',2)
        set(gca,'xlim',plotstruct.xlimits,'ylim',plotstruct.ylimits, ...
            'xtick',plotstruct.xticks,'ytick',plotstruct.yticks)
    end
    col400ch3 = [col400ch3,sep_tmp2(startval:end,2)'];
    if j==9
        [mean400ch3, std400ch3] = compute_stats(col400ch3);
        plot(expossy,mean400ch3,'ko')
        errorbar(expossy,mean400ch3,std400ch3,'k')
    end
end

function [data_mean, data_std] = compute_stats(datain);
    ind = ~isnan(datain);
    data_mean = mean(datain(ind));
    data_std = std(datain(ind));
    

