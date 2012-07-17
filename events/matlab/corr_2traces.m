function [r] = corr_2traces(d1,d2,i1,i2,n1,n2, ilb)

% function to compute the correlation coefficients of data traces d1 and d2
% between samples i1 and i2. 
%
% usage: [r] = corr_2traces(d1,d2,i1,i2)
%
% input:
%         d1 = first time series
%         d2 = second time serie
%         i1 = first sample to use for cross-correlation
%         i2 = last sample to use for cross-correlation
%         n1 = first sample to use for noise (added by DR)
%         n2 = last sample to use for noise (added by DR)
%         ilb = index for lag bounds - typically set to index equivalent of 0.5sec  
%
% output:
%         cor = matrix with correlation coefficients
%         r(i) = average of cor(j,j+i) over j 
%         sigr(i) = standard deviation of this average
%
% Roel Snieder, 12/30/2005


% find index equivalent of 0.05 sec


% get correlation matrix  
for i = -ilb:ilb
%     whos d1
%     i1 
%     i2
    dum =  corrcoef( d1(i1:i2) , d2(i1+i:i2+i) );
    rtmp(i+ilb+1) = dum(1,2);
end
r = max(rtmp);

% carry out noise correction, first determine average signal and noise
% levels. THE NOISE IS DETERMINED FROM THE FIRST 'nnoise' SAMPLES!
nnoise = 190;

noise1 = mean( d1(n1:n2) .* d1(n1:n2) ); % noise1 = mean( d1(1:nnoise) .* d1(1:nnoise) );
signal1 = mean( d1(i1:i2) .* d1(i1:i2) );
noise2 = mean( d2(n1:n2) .* d2(n1:n2) );  % noise2 = mean( d2(1:nnoise) .* d2(1:nnoise) );
signal2 = mean( d2(i1:i2) .* d2(i1:i2) );
if noise1 < signal1 ;
   corfac1 = sqrt( 1 - noise1/signal1 );
else
   corfac1 = 1;
   warning('WARNING: Noise larger than signal, correction factor set to 1!')
end
if noise2 < signal2 ;
   corfac2 = sqrt( 1 - noise2/signal2 );
else
   corfac2 = 1;
   warning('WARNING: Noise larger than signal, correction factor set to 1!')
end

r = r / (corfac1*corfac2);

