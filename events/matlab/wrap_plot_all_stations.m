function [xcorrStruc,sepStruc,domfreq,omega2Struc] =wrap_plot_all_stations(path, fname_u,...
            fname_multi_p,known_sep,Lcut, Hcut,diag_vec,plot_par)

% This is a wrapper for plot_all_stations. It creates a suite of plots for
% the CWI separation estimates for 2D finite difference. 
%
% INPUTS:
% path              [string] full path for location of input files.
% fname_u           [string] filename for unperturbed waves. Note that the
%                   file must satisfy the format of PLOT_ALL_STATIONS. That
%                   is - the file must have 12 columns. The first column is
%                   time, the second to twelth columns are waveforms measured
%                   at different receivers (i.e. the file usually comes from
%                   ac2d_for.f90)
% fname_multi_p     [cell array 1Xn] each element is a filename for a 
%                   perturbed event to be compared with FNAME_U
%
% OUTPUTS:
% xcorrStruc
% sepStruc
% domfreq           [structure] contains the dominant frequency for the
%                   reference and all the perturbed waveforms


for i = 1:length(fname_multi_p)

    disp(['Now doing the ', fname_u(end-17:end-4), '-',fname_multi_p{i}(end-17:end-4), ...
        ' event pair'])
    [xcorrStruc.([fname_u(end-17:end-4), '_',fname_multi_p{i}(end-17:end-4)]) ...
        , h1, ...
        sepStruc.([fname_u(end-17:end-4), '_',fname_multi_p{i}(end-17:end-4)]) ...
        , h2, ...
        extras, omega2Struc] ...
               = plot_all_stations([path,'/',fname_u],...
                 [path,'/',fname_multi_p{i}],known_sep(i), ...
                 Lcut, Hcut,diag_vec,plot_par);
    for j =1:11
        if i ==1
            
            %disp('doing ref')
            
            %domfreq.ref.(['station', num2str(j)]) = extras.(['station', num2str(j)]).end.u.fdom;
            domfreq.(['ref_',fname_u(end-17:end-4)]).(['station', num2str(j)]) = extras.(['station', num2str(j)]).end.u.fdom;
            %domfreq.(['ref_',fname_u(end-17:end-4)])
            %domfreq
            
        end
        tmp1001 = extras.(['station', num2str(j)]).end.p.fdom;
        domfreq.(fname_multi_p{i}(end-17:end-4)).(['station', num2str(j)]) = tmp1001;%extras.(['station', num2str(j)]).end.p.fdom;
    end
    
    moveval=0;
    try, cd diags, moveval = 1, catch, end
    if ~isempty(h1)
        figure(h1)
        print(h1,'-depsc', ['xcorr_All', fname_u(end-17:end-4), '_',fname_multi_p{i}(end-17:end-4), '_',num2str(Hcut), 'to', num2str(Lcut),'Hz.eps']);
        figure(h2)
        print(h2,'-depsc', ['sep_All', fname_u(end-17:end-4), '_',fname_multi_p{i}(end-17:end-4),'_',num2str(Hcut), 'to', num2str(Lcut),'Hz.eps']);
    end
    if moveval ==1, cd .., end
    
end



