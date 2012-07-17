function [fitted_par,fh] = analyse_CWIest_pdftype(data,col, normtype,Hcut,Lcut,outputdir,outputfileprefix,nbins,qa_switch)

% INPUTS: 
% data    [structure or string]
%              if data = structure => use data as it is
%              if data = string => load data from file given by string
% col     [integer = 2,3,4] column to use from CWI separation output
%               2 => Taylo series approximation
%               3 => auto-correlation method (unperturbved wave)
%               4 => auto-correlation method (unperturbved wave)
% normtype [integer = 0,1]
%               0 => do not normalise separation
%               1 => normalise separation by wavelength
% Hcut      [vector] lower bound for Butterworth bandpass used on data.
%           Loops through all values.
% Lcut      [vector] upper bound for Butterworth bandpass used on data. 
% outputdir         [string] full path to director in which to save plots 
%                   or './' 
% outputfileprefix  [string] prefix to be used with filename of saved plots
% nbins     [scalar] number of bins to be used when drawing the histograms
% qa_switch [scalar]
%               1=> print figure
%               0=> do not print figure (note figure created and then
%               closed. 
%
% OUTPUTS
% fitted_par    [structure] stores the fitted parameters in the following
%               form:
%                   fitted_par.for1to5Hz(i,1) => mean for known_sep(i)
%                   fitted_par.for1to5Hz(i,2) => std for known_sep(i)
%                   fitted_par.for1to5Hz(i,2) => amplitude scaling for plotting for known_sep(i)
%   ** printed figures saved to file
% 
% USEAGE
% analyse_CWIest_pdftype('errorbar_data_col3_normtype1.mat',3,1,[1 2 3 4],[2 3 4 5],'./','pdf_hist_')


if isstruct(data)
    errorbar_data = data;
elseif ischar(data)
    load(data)
end

if length(Hcut) ~=length(Lcut)
    error('Length of Hcut must equal length of Lcut in ANALYSE_CWIEST_PDFTYPE')
end

%% Original inputs when this was a script
% col = 3; % 2=> use taylor, 3=> use autocorrelation (unpert.), 4=> use autocorrelation (pert.)
% normtype = 1; % 0 = > use actual separation estimates or 1=> normalise by wavelength (domfreq/velocity)
% fname = ['errorbar_data_col',num2str(col),'_normtype', num2str(normtype),'.mat'];
% load(fname)
% Hcut = [1 2 3 4];  % 1
% Lcut = [2 3 4 5];  % 5
%% Setup the Bounds for the fitting process
if normtype ==1
    lbounds = [0 0 0];
    upperbounds = [1 1 10000];
elseif normtype ==0
    lbounds = [0 0 0];
    upperbounds = [2000 2000 10000];
end



colors = {'r','b','g','k'};
markers = {'o','^','s','d'};

for j = 1:length(Hcut)
    folded_NPDF.known_sep_perfreq.(['for_',num2str(Hcut(j)),'to',num2str(Lcut(j)),'Hz']) = ...
        errorbar_data.known_sep_perfreq.(['for_',num2str(Hcut(j)),'to',num2str(Lcut(j)),'Hz']);
    for i = 1:15
        fh(i) = figure
        data = errorbar_data.agg_data.(['for_',num2str(Hcut(j)),'to',num2str(Lcut(j)),'Hz']).(['knownsep',num2str(i)]);
        %hist(data,100);
        [N X] = hist(data,nbins);
        Narea = trapz(X,N);
        %N = N/Narea;
        bar(X(2:end),N(2:end),1)
        %Nnorm = N/Narea;
        hold on

        % figure
        % plot(X,Narea)

        % Now let's try fitting
        %[out1,resnorm,residual,exitflag,output] = lsqcurvefit(@wrap_folded_normal_pdf,[mean(data),std(data),1],X,N);
        [out1,resnorm,residual,exitflag,output] = lsqcurvefit(@wrap_folded_normal_pdf,[mean(data),std(data),1],X(2:end),N(2:end),lbounds, upperbounds);
        fitted_par.(['for_',num2str(Hcut(j)),'to',num2str(Lcut(j)),'Hz'])(i,:) =  out1;
            
         if median(data) == out1(1)
             disp('Warning median(data) == out1(1)')
         end
         if std(data) == out1(2)
             disp('Warning std(data) == out1(2)')
         end
         
         
        Y = folded_normal_pdf(X,out1(1),out1(2));
        Y = Y/max(Y)*max(N(2:end));
        l1 = plot(X,Y,'r','linewidth',2)

        folded_NPDF.mean.(['for_',num2str(Hcut(j)),'to',num2str(Lcut(j)),'Hz'])(i) = out1(1); 
        folded_NPDF.std.(['for_',num2str(Hcut(j)),'to',num2str(Lcut(j)),'Hz'])(i) = out1(2); 

        
        norm_mu = mean(data);
        norm_sigma = std(data);
        Y2 = normpdf(X,norm_mu,norm_sigma);
        Y2 = Y2/max(Y2)*max(N(2:end));
        l2 = plot(X,Y2,'g','linewidth',2)
        set(gca,'fontsize',35)
        
        if qa_switch ==1
            set(gca,'ylim',[0 100],'xlim',[0 1])
            known_seps = [57, 113,170,226,255,283,311,339,368,396,424,453,566,679,792];
            text(0.55,85,['\delta_t = ',num2str(known_seps(i)),' m'],'fontsize',35)
            print('-depsc',[outputdir, filesep, outputfileprefix,'norm', num2str(normtype),'_col', num2str(col), '_', num2str(Hcut(j)),'to',num2str(Lcut(j)),'Hz-','knownsep',num2str(i),'.eps'])
        else
            close(fh(i))
        end
    end
end
%legend([l1,l2],{'folded normal PDF','normal PDF'})

fname = ['folded_errorbar_data_col',num2str(col),'_normtype', num2str(normtype),'.mat'];
save(fname, 'folded_NPDF')

