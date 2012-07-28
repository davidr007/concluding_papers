function EQdata4gmt(filein,fileout)
%Loads a csv file created from the online GA earthquake catalogue and 
%creates a file for convenient plotting in GMT.
%
%NOTE - you must remove the column "Approximate Location" 
% from the CSV file because it has extra commas in it which cause
%problems with textread.
%
%INPUTS
%filein     [string] path to input file created from the online GA
%           earthquake database
% fileout   [string] full path name of output file
%
% David Robinson
% 25 July 2012


if nargin ==0
    filein ='GAcat_1Jan1992to24July2012_M2_10.csv';
    fileout = 'temp.dat';
end

% Read the CSV file
fid = fopen(filein);
tmp = textscan(fid,'%s%s%s%s%s%s%s%s%s%s%s%s%s%s','delimiter',',','EmptyValue', -Inf);
%tmp{1} = magnitude
%tmp{2} = UTC date
%tmp{3} = UTC Time
%tmp{4} = Sydney Date
%tmp{5} = Sydney Time
%tmp{6} = Latitude
%tmp{7} = Longitude
%tmp{8} = Significant
%tmp{9} = Approximate location
%tmp{10} = Depth (km)
%tmp{11} = Solution finalised
%tmp{12} = Source
%tmp{13} = EVENT ID
%tmp{14} = ORIGIN ID



magnitudes = str2num(char(tmp{1}(2:end)));
longitudes = str2num(char(tmp{7}(2:end)));
latitudes = str2num(char(tmp{6}(2:end)));

magnitudeMarker = zeros(size(magnitudes));
magnitudeMarker(magnitudes>=5) = 0.12;
magnitudeMarker(magnitudes<5) = 0.1;
magnitudeMarker(magnitudes<4) = 0.08;
magnitudeMarker(magnitudes<3) = 0.04;

out = [longitudes,latitudes,magnitudeMarker];
save eq_1992to2012_GA.dat out -ascii


