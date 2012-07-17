function [VariablesBinned] = bin(Variables,BinCentroids);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bin.m takes the the elements in Variables and places them in bins based on the 
% elements in BinCentroid. Note that the end points of each bin are asumed to 
% be half way between the bin centroids. The first and last bins extend to negative
% and positive infinity respectively.
%
% Inputs:
% Variables         A vector of variables to be binned
% BinCentroids      A vector of bin centroids
%
% Output:
% VariablesBinned   A 2 column matrix The second column contains the binned variables. 
%                   Note that each of the values in VariablesBinned(:,2) take on the value
%                   of a bin centroid (BinCentroids). The first column contains a logical index
%                   whose number is a pointer to the relevant index in BinCentroid.
%
%                                                       David Robinson
%                                                       Last Update: 19/03/02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(BinCentroids)==1
    VariablesBinned = [ones(length(Variables),1),BinCentroids.*ones(length(Variables),1)];
else
    
    Binends = [BinCentroids(1:length(BinCentroids)-1)+diff(BinCentroids)./2];
    index = zeros(length(Variables),1);   % initialising index vector
    VariablesBinned = zeros(length(Variables),1);  % initialising the VariablesBinned vector
    
    
    % the first bin
    I = find(Variables<Binends(1));
    VariablesBinned(I) = BinCentroids(1);
    index(I) = 1;
    clear I
    %[index(1:20) Variables(1:20) VariablesBinned(1:20)]   % diagnostics
    
    
    % all middle bins
    for ind1 = 2: length(Binends)
        I = find(Variables>=Binends(ind1-1) & Variables<Binends(ind1));
        index(I) = ind1;
        VariablesBinned(I) = BinCentroids(ind1);
        %[index(1:20) Variables(1:20) VariablesBinned(1:20)]   % diagnostics
        clear I
    end
    
    %Binends    % diagnostics
    %BinCentroids  % diagnostics
    % the last bin
    ind1 = length(BinCentroids);
    I = find(Variables>=Binends(length(Binends)));
    VariablesBinned(I) = BinCentroids(ind1);
    index(I) = length(Binends)+1;
    %[index(1:20) Variables(1:20) VariablesBinned(1:20)]   % diagnostics
    %wow = Binends(length(Binends))   % diagnostics
    clear I
    
    VariablesBinned = [index VariablesBinned];
end