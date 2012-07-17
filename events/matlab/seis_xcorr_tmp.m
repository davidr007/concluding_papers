function [norm_maxcorr,sep, tmp_sep, tmp_xcorr,omega2] = seis_xcorr(wave1, wave2, time1, time2,plot_par,windowtype,diag_switch)
% This is the worhouse of the CWI method. It computes the normalised cross
% correlation between two waveforms and uses it to estimate the separation
% (or RMS change in source orientation) between the two events.
%
% Equation reference are from Snieder and Vrijlandt (2005), JGR.
%
% INPUTS:
%
% wave1 & wave2     The two seismic waves to be analysed. Note that atleast
% one of
%                   them must be horizontally (time) translated so that they align well.
%                   Events must have the same sample rate.
% time1 & time2     The horizontally translated time vectors for events 1 and 2, respectively.
% plot_par          Structure controlling application of CWI
%                       plot_par.twindow =  width of sliding window (sec)
%                       plot_par.tpu =      time of P arrival on unperturbed wave
%                       plot_par.tpp =      time of P arrival on perturbed wave
%                       plot_par.ztapper =  length of ztapper (index) used
%                                           to filter waves
%                       plot_par.trans1 =   translation to apply to
%                                           unperturbed wave (typically
%                                           this is zero if using
%                                           ALLIGN_WAVEFORMS
%                       plot_par.trans2 =   translation to apply to
%                                           perturbed wave
%                       plot_par.start_arrivwin =   Optional Diagnostics - start of window
%                                                   for plotting the direct arrival
%                       plot_par.end_arrivwin =     Optional Diagnostics - start of window
%                                                   for plotting the direct
%                                                   arrival
%                       plot_par.start_codawin =    Optional Diagnostics - start of window
%                                                   for plotting the coda
%                       plot_par.end_codawin =      Optional Diagnostics - end of window for
%                                                   plotting the coda
%                       plot_par.parrival =         ??? NOT OPERATIONAL
%                       plot_par.event1name =       Optional Diagnostics - label used for
%                                                   un-perturbed event
%                       plot_par.event2name =       Optional Diagnostics - label used for
%                                                   perturbed event
%                       plot_par.rem_noise =        Control the removal of noise
%                                                   0=> do not remove noise (for Synthetics)
%                                                   1=> QUICK AND DIRTY -  noise computed from
%                                                       waveform cut-out which is of length
%                                                       indwinlen and which is located in the
%                                                       centre between tp and the start of the
%                                                       waveform (see code for potential flaws)
%                                                   2=> PREFFERED - uses size of waveform to
%                                                       compute the noise and avoids overlap
%                                                       with ztapper and tp.%                                                   m
%                                                   3=> QUICK AND DIRTY - removes noise starting
%                                                       at the beggining of the record using one
%                                                       window length (indwinlen).
%                       plot_par.lbnd4nenergy =     ??? NOT OPERATIONAL
%                       plot_par.ubnd4nenergy =     ??? NOT OPERATIONAL
%                       plot_par.src_type =         Type of excitation and CWI
%                                                   1 =>    DISPLACEMENT -  point source in a 2D
%                                                                           acoustic medium
%                                                   2 =>    DISPLACEMENT -  double couple source in
%                                                                           a 3D elastic medium
%                                                                           displaced in the source
%                                                                           plane
%                                                   3 =>    SOURCE VARIATION - double couple - co-located
%                                                                              and rotated
%                       plot_par.mean_vs =         Average S-wave velocity near source
%                       plot_par.mean_vp =         Average P-wave velocity near source
% windowtype        [case switch]
%                   1) 'independent' or 'i'
%                       =>  independent window (non-overlapping)
%                   2) 'non independent' or 'ni'
%                       =>  non independent windows (overlapping) sliding by
%                           one index at a time.
% diag_switch       [integer] - OPTIONAL set to 0 if not provided
%                           1 => draw plots
%                           0 => do not draw plots
% OUTPUTS:
% norm_maxcorr      [nx4] maximum normalised cross correlation (col2) as a
%                   function of time at the centroid of the n sliding windows.
%                   Note that maximum is taken over all lag times in the
%                   sliding window. Noise removal depends on value of
%                   plot_par.rem_noise.
% sep               [nx4] separation or RMS  change in source orientation
%                     ** separation in metres (src_type = 1 or 2)
%                       col1 => time at sliding window centroid (as with norm_maxcorr)
%                       col2 => separation using taylor series approximation (eq6)
%                       col3 => separation estimate using auto-correlation of the unperturbed
%                               wave (wave1)
%                       col4 => separation estimate using auto-correlation of the perturbed
%                               wave (wave2)
%                     ** RMS  change in source orientation (src_type =3)
%                       col1 => time at sliding window centroid (as with norm_maxcorr)
%                       col2 => RMS change in source orientation
%                       col3 => NOT APPLICABLE
%                       col4 => NOT APPLICABLE
% tmp_sep           [nx4] Extra information - useful for debugging
%                       col1 => omega2_numerator (eq8)
%                       col2 => omega2_denominator (eq8)
%                       col3 => 2*(1-maxcorr) (re-arrangement of eq6)
%                       col4 => omega2 (eq8)
% tmp_xcorr         [nx4] Extra information - useful for debugging
%                       col1 => Value of cross correlation at location of norm_maxcorr(i).
%                               That is, it is the maximum of the numerator of eq3 without
%                               normalisation or noise removal.
%                       col2 => Denominator of eq3
% omega2            [nx1] omega^2 values as a function of sliding window (i.e. as used to
%                   compute the separation)
%
% David Robinson
% April 2005
% Major Overhaul - August 2007


if nargin <7 | ~exist('diag_switch')
    diag_switch = 0
end

timewindow = plot_par.twindow;
sample_rate = abs(time1(2)-time1(1));  % Must be = abs(time2(2)-time2(1))
if sample_rate-abs(time2(2)-time2(1))>10^(-10);
    error('Sample rate of event1 and event2 must be the same in SEIS_XCORR');
end

%length of time window in vector index space
indwinlen = find_closest(time1-time1(1),timewindow,'euclidean');


%number of independent (non-overlapping windows)
%ERRORRRRRRRR - should I have an extra set of brackets inside the floor????
num_indwindows = min(   [floor((time1(end)-time1(1))/timewindow)-1, ...
    floor((time2(end)-time2(1))/timewindow)-1]);
% Now we just check that it works (some waveforms cause difficulties)
tester = 0;
while tester == 0
    if indwinlen + num_indwindows*indwinlen > max([length(time1),length(time2)])
        % test fails reduce num_indwindows by 1
        num_indwindows = num_indwindows-1;
    else % test passes
        tester =1;
    end
end

[noise_energy1,noise_energy2,ind_noisestart,ind_noiseend] = LOC_compute_noise(wave1,wave2,time1,time2,indwinlen,plot_par);
%maximum number of windows (sliding by one index at a time)
%num_slidwindows = min([length(time1)-indwinlen,length(time2)-indwinlen]);
num_slidwindows = min([length(wave1),length(wave2)])-2*indwinlen-1;


indsecondlen = find_closest(time1-time1(1),1,'euclidean');
num_secondwindows = floor(length(time1)/indsecondlen)-timewindow;

if ~isempty(noise_energy1) & ~isempty(noise_energy2)  %make sure that no problems were identified in LOC_compute_noise

    if num_slidwindows > 2  % make sure there are enough sliding windows
        indstart = 0;
        switch windowtype
            case {'i','independent'}
                maxcorr = zeros(num_indwindows,2);
                sep = zeros(num_indwindows,4);
                numwindows = num_indwindows;
                for i=1:numwindows
                    indstart(i) = (i-1)*indwinlen+1 + floor(indwinlen/2);
                    indend(i) = i*indwinlen + ceil(indwinlen/2);   %not this is a ceil (above is floor) to ensure window length = indwinlen
                    maxcorr(i,1) = (i-1)*timewindow + timewindow/2;
                    sep(i,1) = maxcorr(i,1);
                    [maxcorr(i,2), tmp_xcorr(i,:)] = LOC_compute_correlation(wave1,wave2,indstart(i),indend(i),time1,time2,plot_par,noise_energy1,noise_energy2,ind_noisestart,ind_noiseend);
                end
            case {'non independent','ni'}
                maxcorr = zeros(num_slidwindows,2);
                sep = zeros(num_slidwindows,4);
                numwindows = num_slidwindows;
                for i=1:numwindows
                    indstart(i) = (i-1)+1 + floor(indwinlen/2);
                    indend(i) = indstart(i) +indwinlen;
                    maxcorr(i,1) = (i-1)*sample_rate + timewindow/2;       %(i-1)*timewindow + timewindow/2;
                    sep(i,1) = maxcorr(i,1);
                    [maxcorr(i,2), tmp_xcorr(i,:)] = LOC_compute_correlation(wave1,wave2,indstart(i),indend(i),time1,time2,plot_par,noise_energy1,noise_energy2,ind_noisestart,ind_noiseend);
                end
            case {'seconds', 's'}
                maxcorr = zeros(num_secondwindows,2);
                sep = zeros(num_secondwindows,4);
                numwindows = num_secondwindows;
                for i=1:numwindows
                    indstart(i) = (i-1)*indsecondlen+1 + 10;  % 10 = slight offset to account time lagging of windows
                    indend(i) = indstart(i) +indwinlen;
                    maxcorr(i,1) = mean([time1(indstart(i)),time1(indend(i)) ]);%(i-1)*sample_rate + timewindow/2;       %(i-1)*timewindow + timewindow/2;
                    sep(i,1) = maxcorr(i,1);
                    [maxcorr(i,2), tmp_xcorr(i,:)] = LOC_compute_correlation(wave1,wave2,indstart(i),indend(i),time1,time2,plot_par,noise_energy1,noise_energy2,ind_noisestart,ind_noiseend);
                end
                
                
      
            otherwise
                error('Inavlid value for input WINDOWTYPE in SEIS_XCORR');
        end
        norm_maxcorr = maxcorr;

        for i=1:numwindows
            [sep(i,2),sep(i,3), sep(i,4), tmp_sep(i,:), omega2(i)] = LOC_compute_distance(wave1,wave2,time1,sample_rate,indstart(i),indend(i),norm_maxcorr(i,2), plot_par, diag_switch);
        end

    else % the waveforms are too short for the size of the sliding window
        norm_maxcorr=[];
        sep=[];
        tmp_sep=[];
        tmp_xcorr =[];
        omega2=[];
    end

else    % i.e. problems were identified in LOC_compute_noise
    norm_maxcorr=[];
    sep=[];
    tmp_sep=[];
    tmp_xcorr =[];
    omega2=[];
end     

% =====================================================
%% Managers the guts of the work for each time window
function [maxcorr,sep] = LOC_do_work(wave1,wave2,indstart,indend,time1,time2,sample_rate)
[maxcorr] = LOC_compute_correlation(wave1,wave2,indstart,indend,time1,time2);
[sep] = LOC_compute_distance(wave1,time1,sample_rate,indstart,indend,maxcorr);

% =====================================================
%% Compute the the maximum normalised cross correlation for each time window
function [maxcorr, tmp_xcorr] = LOC_compute_correlation(wave1,wave2,indstart,indend,time1,time2,plot_par,noise_energy1,noise_energy2,ind_noisestart,ind_noiseend)
%count = 0;
num_ts = indend-indstart;
count = -floor(num_ts/2);

% find index equivalent of 0.05sec
ind = find(time1-time1(1)>0.05);


[maxcorr] = corr_2traces(wave1,wave2,indstart,indend,ind_noisestart,ind_noiseend,ind(1));
tmp_xcorr = [NaN NaN];


% C = zeros(1,num_ts);
% for i=1:num_ts
%     C(i) = trapz(time1(indstart:indend),wave1(indstart:indend).*wave2(indstart+count:indend+count));
%     [energy1(i)] = compute_energy(wave1(indstart:indend),time1(indstart:indend));
%     [energy2(i)] = compute_energy(wave2(indstart+count:indend+count),time2(indstart+count:indend+count));
%     count = count+1;
% end

% ind = find(energy1.*energy2>0);
% if ~isempty(ind)
%     if plot_par.rem_noise == 0
%         normC = C./sqrt(energy1.*energy2);
%     elseif plot_par.rem_noise ==1 | plot_par.rem_noise ==2 | plot_par.rem_noise ==3
%         %normC = (C./sqrt(energy1.*energy2))./(sqrt(1-noise_energy1^2./energy1.^2).*sqrt(1-noise_energy2^2./energy2.^2));
% 
%         tmp1 = noise_energy1^2./energy1.^2;
%         noiseP1 = zeros(size(tmp1));
%         noiseP1(tmp1<=1) = sqrt(1-noise_energy1^2./energy1(tmp1<=1).^2);
%         noiseP1(tmp1>1) = NaN;
%         tmp2 = noise_energy2^2./energy2.^2;
%         noiseP2 = zeros(size(tmp2));
%         noiseP2(tmp2<=1) = sqrt(1-noise_energy2^2./energy2(tmp2<=1).^2);
%         noiseP2(tmp2>1) = NaN;
%         normC = (C./sqrt(energy1.*energy2))./(noiseP1.*noiseP2);
%         if ~isempty(find(tmp1>1)) | ~isempty(find(tmp2>1))
%             warning('WARNING: noise energy exceeds signal energy');
%         end
% 
%     else
%         error('Invalid value for plot_par.rem_nois in SEIS_XCORR')
%     end
% 
%     ind = find(normC<=1);
%     if ~isempty(ind)
%         [maxcorr,indmax] = max(normC(ind)); %??? note have tried max(abs(normC)) here as well but this is not (?) correct because we are not interested in anti-correlation
%         tmp_xcorr = zeros(1,2);
%         tmp_xcorr = [C(indmax), sqrt(energy1(indmax).*energy2(indmax))];
% 
%     else
%         maxcorr = NaN;
%         tmp_xcorr = [NaN NaN];
%     end
% else
%     maxcorr =NaN;
%     tmp_xcorr = [NaN NaN];
% end


% =====================================================
%% Compute the distance between events for each time window
function [sep_taylor,sep_auto_u, sep_auto_p, tmp_sep, omega2] = LOC_compute_distance(wave_u,wave_p,time_u,sample_rate_td,indstart,indend,maxcorr,plot_par, diag_switch);
% note maxcorr is a scalar here
% energy_numerator = gradient(wave1,sample_rate);
grad_wave_u =gradient(wave_u,sample_rate_td);
omega2_numerator = mean(grad_wave_u(indstart:indend).*grad_wave_u(indstart:indend)); %compute_energy(grad_wave_u(indstart:indend),time_u(indstart:indend));
omega2_denominator= mean(wave_u(indstart:indend).*wave_u(indstart:indend));%compute_energy(wave_u(indstart:indend),time_u(indstart:indend));

nu = length(wave_u);

if omega2_denominator~=0 & ~isnan(maxcorr) & maxcorr > 0  % last condition is required to deal with some emergent signals?
    omega2 = omega2_numerator/omega2_denominator;
    if maxcorr >1 & maxcorr <1.0001  % deal with rounding errors when maxcorr is roughly 1
        maxcorr = 1;
    end
    sigma2_taylor = 2*(1-maxcorr)/omega2; % from Eq. 6

    [autocorr_u,lags_u] = xcorr(wave_u);
    lags_u = lags_u*sample_rate_td;
    autocorr_u = autocorr_u/max(autocorr_u);
    ind = find(autocorr_u<maxcorr);  % find everywhere where the a-corr<0
    ind2 = find(ind>nu);
    xinter_top = ind(ind2(1));  % lag index just to the right of first intercept with maxcorr
    %sigma2_auto = (lags_u(xinter_top) + lags_u(xinter_top-1))/2;
    % now do the linear interpolation
    slope_u = (autocorr_u(xinter_top)- autocorr_u(xinter_top-1))/(lags_u(xinter_top)- lags_u(xinter_top-1));
    intercept_u = autocorr_u(xinter_top) - slope_u*lags_u(xinter_top);
    sigma_auto_u = (maxcorr-intercept_u)/slope_u;
    sigma2_auto_u = sigma_auto_u^2;

    [autocorr_p,lags_p] = xcorr(wave_p);
    lags_p = lags_p*sample_rate_td;
    autocorr_p = autocorr_p/max(autocorr_p);
    ind = find(autocorr_p<maxcorr);  % find everywhere where the a-corr<0
    ind2 = find(ind>nu);
    xinter_top = ind(ind2(1));  % lag index just to the right of first intercept with maxcorr
    %sigma2_auto = (lags_u(xinter_top) + lags_u(xinter_top-1))/2;
    % now do the linear interpolation
    slope_p = (autocorr_p(xinter_top)- autocorr_p(xinter_top-1))/(lags_p(xinter_top)- lags_p(xinter_top-1));
    intercept_p = autocorr_p(xinter_top) - slope_p*lags_p(xinter_top);
    sigma_auto_p = (maxcorr-intercept_p)/slope_p;
    sigma2_auto_p = sigma_auto_p^2;

    if diag_switch ==1
        f1 = figure;
        ind = find(autocorr_u<0);  % find everywhere where the a-corr<0
        ind2 = find(ind>nu);
        xinter_top = ind(ind2(1));  % just to the right of the first (high) x-intercept
        ind3 = find(ind<nu);
        xinter_bot = ind(ind3(end));  % just to the left of the first (low) x-intercept
        lags4plot_u= linspace(lags_u(xinter_top), lags_u(xinter_bot),1000);
        R_est_u = 1-omega2*lags4plot_u.^2/2;
        h1 = plot(lags_u,autocorr_u/max(autocorr_u),'b');
        hold on
        h2 = plot(lags_p,autocorr_p/max(autocorr_p),'r');
        plot([lags_u(xinter_bot), lags_u(xinter_top)], [maxcorr maxcorr], '--')
        plot(sigma_auto_u,maxcorr,'bo')
        plot(sigma_auto_p,maxcorr,'ro')
        h3 = plot(lags4plot_u,R_est_u,'k');
        plot(sqrt(sigma2_taylor) , maxcorr, 'ko');
        set(gca,'xlim',[lags_u(xinter_bot), lags_u(xinter_top)])
        set(gca,'ylim',[0,1])
        xlabel('omega')
        ylabel('Correlation')
        legend([h1,h2,h3],{'EU AutoCorr', 'EP AutoCorr', 'Taylor'})
        title(['omega: T = ', num2str(sqrt(sigma2_taylor)),' EU = ', num2str(sigma_auto_u), ' EP = ', num2str(sigma_auto_p)])
        pause
        close(f1)
    end



    switch plot_par.src_type
        case 1
            sep_taylor = sqrt(2*6000^2*sigma2_taylor);  % from Eq. 54
            sep_auto_u = sqrt(2*6000^2*sigma2_auto_u);
            sep_auto_p = sqrt(2*6000^2*sigma2_auto_p);
        case 2

            vs = plot_par.mean_vs;
            vp = plot_par.mean_vp;
            part1 = 7.0*(2.0/(vp^6)+3.0/(vs^6));
            part2 = 6.0/(vp^8)+7.0/(vs^8);
            part3_taylor  = part1*sigma2_taylor/part2;
            sep_taylor = sqrt(part3_taylor);
            part3_auto_u  = part1*sigma2_auto_u/part2;
            sep_auto_u = sqrt(part3_auto_u);
            part3_auto_p  = part1*sigma2_auto_p/part2;
            sep_auto_p = sqrt(part3_auto_p);

        case 3  % Note that becuase we do not have the auto-correlation here
            % we assign all values to the same appr
            sep_taylor = maxcorr -1;
            sep_auto_p = sep_taylor;
            sep_auto_u = sep_taylor;
            tmp_sep = [sep_taylor sep_taylor sep_taylor sep_taylor];

    end

    if sigma2_taylor<0
        disp('sigma2<0')
        try, disp(['calculated sep is: ', num2str(sep)]), catch, end
        disp('sep re-assigned to NaN')
        sep_taylor = NaN;
    end
else
    sep_taylor = NaN;
    omega2 = NaN;
    sep_auto_u = NaN;
    sep_auto_p = NaN;
end
tmp_sep = zeros(1,4);
tmp_sep = [omega2_numerator, omega2_denominator, 2*(1-maxcorr),omega2];

function     [noise_energy1,noise_energy2,ind_noisestart, ind_noiseend] = LOC_compute_noise(wave1,wave2,time1,time2,indwinlen,plot_par)
% The purpose of this sub-function is to compute the noise
%Inputs:
% wave1         [vector] - unperturbed wave
% wave2         [vector] - perturbed wave
% time1         [vector] - time for unperturbed wave
% time2         [vector] = time for perturbed wave
% indwinlen     [scalar] window length in index space
% plotpar       [structure] - standard structure which includes
%                   plot_par.twindow = width of timewindow sec)
%                   plot_par.tpu = time of P arrival on unperturbed wave
%                   plot_par.tpp = time of P arrival on perturbed wave
%                   plot_par.ztapper = length of ztapper (index) used
%                   to filter waves

sample_rate = time1(2)-time1(1);
errorhandle = 0;  %initialise an error handle - if this is set to 1 them all outputs of LOC_compute_noise are empty

ind_noisestart = [];
ind_noiseend = [];

if plot_par.rem_noise == 0  % we do not remove any noise
    noise_energy1 = 0;
    noise_energy2 = 0;
else % we remove the noise
    % First we find some important locations in the signals
    indtpu = find_closest(time1,plot_par.tpu,'euclidean');
    indtpp = find_closest(time2,plot_par.tpp,'euclidean');
    indeztaper = find_closest(time1,plot_par.ztapper*sample_rate + time1(1),'euclidean');
    if plot_par.rem_noise ==1  % quick and dirty option
        % noise computed from waveform cut-out which is of length
        % indwinlen and which is located in the centre between tp and
        % the start of the waveform
        % KNOWN PROBLEM - if indwinlen exceeds 1/2(tp+indeztaper) then
        % window used to compute noise will extend past the parrival.
        indstartu = ceil((indtpu+indeztaper)/2)-ceil(indwinlen/2);
        indendu = indstartu+indwinlen;
        ind_lennoise_u = indwinlen;
        indstartp = ceil((indtpp+indeztaper)/2)-ceil(indwinlen/2);
        indendp = indstartp+indwinlen;
        ind_lennoise_p = indwinlen;
    elseif plot_par.rem_noise == 2  % PREFFERED technique
        %   uses the maximum ammount of data available for computing the noise
        %   noise is later re-normalised so that it is representative
        %   of noise in a single window of length indwinlen
        %   Note that this option will also scale up i.e. if indwinlen>
        %   tp-startwave
        ind5percentu = ceil(0.05*(indtpu-indeztaper));  %    used to ensure noise is
        %                                                   completely clear of the
        %                                                       end of the ztapper and the
        %                                                       start of the p-arrival
        ind5percentp = ceil(0.05*(indtpp-indeztaper));
        if ind5percentp<0 | ind5percentu<0  % check to make sure that ztapper does not cross p-arrival
            %disp('Ztapper length cannnot be larger than tp-startwave')
            %error('Problem in noise calculation')
            errorhandle = 1;
        end
        indstartu = indeztaper + ind5percentu;
        indendu = indtpu - ind5percentu;
        ind_lennoise_u = indendu-indstartu;   % used for re-normalisation of noise1 to window length
        indstartp = indeztaper + ind5percentp;
        indendp = indtpp - ind5percentp;
        ind_lennoise_p = indendp-indstartp;
        if indstartu~=indstartp
            indstartu
            indstartp
            error('problem in noise calc - 1')
        else
            ind_noisestart = indstartu;
        end
        if indendu~=indendp
            indendu
            indendp
            error('problem in noise calc - 2')
        else
            ind_noiseend = indendu;
        end
    elseif plot_par.rem_noise == 3 % quick and dirty -
        % take frst time window in waveform
        % KNOWN PROBLEM - this approach computes noise from a window that
        % overlaps the ztapper used in filtering
        indstartu = 1;
        indendu = indstartu+indwinlen;
        ind_lennoise_u = indwinlen;
        indstartp = 1;
        indendp = indstartp+indwinlen;
        ind_lennoise_p = indwinlen;
    else
        error('ERROR: Invalid value of plot_par.rem_noise')
    end
    % Now we can compute the noise
    if indendu > indstartu
        noise_energy1 = compute_energy(wave1(indstartu:indendu),time1(indstartu:indendu));
        noise_energy1 = noise_energy1*indwinlen/ind_lennoise_u; %ensure that computed noise is consistenet with window length
        noise_energy2 = compute_energy(wave2(indstartp:indendp),time2(indstartu:indendp));
        noise_energy2 = noise_energy2*indwinlen/ind_lennoise_p; %ensure that computed noise is consistenet with window length
    else
        noise_energy1 = [];
        noise_energy2 = [];
    end

    if errorhandle ==1 % there is a problem so we must re-set the outputs to empty
        noise_energy1 = [];
        noise_energy2 = [];
    end
end