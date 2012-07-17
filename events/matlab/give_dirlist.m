function [filelist, dirlist] = give_dirlist(directory)
%
% Gives a listing of all files and directories in a directory.
% Note that '.' and '..' are removed/
%
% INPUTS:
% directory     [string] directory to be searched
%
% OUTPUTS:
% filelist      [structure array] list of files in the directory
% dirlist       [structure array] list of directories in directory
%
% David Robinson
% 15 November 2006

filelist = [];
counterfl = 1;

dirlist = [];
counterdl = 1;
tmp = dir(directory);
for i = 1:length(tmp)
   if ~strcmp(tmp(i).name,'.') & ~strcmp(tmp(i).name,'..')
      if tmp(i).isdir ==1  & ~strcmp(tmp(i).name,'.svn')% add directory name to dirlist
          dirlist(counterdl).name = tmp(i).name;
          counterdl = counterdl +1;
      elseif  tmp(i).isdir ~=1  & ~strcmp(tmp(i).name,'.svn')% add file name to filelist
          filelist(counterfl).name = tmp(i).name;
          counterfl = counterfl +1;
      end
   end    
end