function signal_out = apply_butter(signal_in, butter_order,type, domain,cut_freq, sample_rate)
% signal_in
% butter_order      [integer] order of the filter e.g. 5 or 9
% type              [string]
%                       'low' = low pass
%                       'high' = high pass
%                       'band' = band pass (NOTE - does not seem stable???)
% domain            [string]
%                       'time' = do filtering in time domain
%                       'freq' = do filtering in frequency domain
% cut_freq          [double] cut-off frequency
% sample_rate       sampling rate in Hz



if nargin ==1 % only the input signal provided
    butter_order = 5;
    type = 'low';
    domain = 'time';
    cut_freq = 5;
    sample_rate = 200;
end

%n = 100000;
%signal_in_long =[zeros(1000,1); signal_in(1001:end)];
signal_in_long = signal_in;


TotalSamps = length(signal_in_long);
Nyq = sample_rate/2;
Wn = cut_freq/Nyq;
    
    
if exist('butter')==0
    error('Signal Processing Toolbox is required to use APPLY_BUTTER')
end

switch type
    case 'low'
    [b,a]=butter(5,Wn);
  
    case 'high'
    [b,a]=butter(5,Wn,'high');
    
    case 'band'
        if length(Wn)~=2
            error('Length of Wn must be 2 for a bandpass filter')
        end
    display('WARNING - I have had problems with bandpass (try low then high pass)')
    [b,a]=butter(5,Wn);
    %[b,a]=cheby1(5,0.5,Wn);
end

    
 switch domain
    case 'freq'
        [h,w]=freqz(b,a, 'whole', TotalSamps, 200);
        FFTinputmAcc1 = hLowPass.*fft(inputmAccOld(:,2));   %required for filter
        inputmAcc1 = real(ifft(FFTinputmAcc1));  
    case 'time'     % apply filter in time domain (using filter.m)
        signal_out = filter(b,a,signal_in_long);  
end   

% if exist('butter')~=0   % Apply butterworth filter if signal processing toolbox is installed
%     % setup
%     inputmAccOld=inputmAcc;
%     Nyq = 1/(2*SampleRate);
%     FFTinputmAccOld = fft(inputmAccOld(:,2));
%     f=[0:(TotalSamps/2)]./(TotalSamps*SampleRate);
% 

%signal_out = signal_in_long(n+1:end);
  