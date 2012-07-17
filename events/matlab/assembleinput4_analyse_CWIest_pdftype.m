function [assembled_data] = assembleinput4_analyse_CWIest_pdftype(datain, knownsep, Hcut, Lcut)
% The purpose of this function is to assemble data for use with analyse_CWIest_pdftype
%
% INPUTS: 
% datain        [matrix] 1 column for each value in knownsep. Rows contain the
%               CWI separation estimates of interest
% knownsep      [vector]  contains all of the known separations
% Hcut          [vector] lower bound for Butterworth bandpass used on data.
%               Loops through all values.
% Lcut          [vector] upper bound for Butterworth bandpass used on data. 
%
%
% OUTPUTS: 
% assembled_data [structure] data input for analyse_CWIest_pdftype
%
% Note that this can only be used for 1 frequency band at a time. Any
% organisation of 

% David Robinson

[n m] = size(datain.((['for_',num2str(Hcut(1)),'to', num2str(Lcut(1)),'Hz']))); 
for j = 1:length(Hcut)
    assembled_data.known_sep_perfreq.(['for_',num2str(Hcut(j)),'to', num2str(Lcut(j)),'Hz']) = ...
        knownsep.(['for_',num2str(Hcut(j)),'to', num2str(Lcut(j)),'Hz']); 
    for i = 1:m
        ind = ~isnan(datain.(['for_',num2str(Hcut(j)),'to', num2str(Lcut(j)),'Hz'])(:,i));
         assembled_data.agg_data.(['for_',num2str(Hcut(j)),'to', num2str(Lcut(j)),'Hz']).(['knownsep',num2str(i)]) = ...
             datain.(['for_',num2str(Hcut(j)),'to', num2str(Lcut(j)),'Hz'])(ind,i);
    end
end

    
