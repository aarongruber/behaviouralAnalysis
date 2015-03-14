function dealBEHIntoDirectory(varargin)
%
% takes .mat files in the specified directory (or present directory if no inputs) and 
% places each one in a directory of the
% same name as the file.
%

if nargin==1
    path = varargin{1};
else
   path = pwd; 
end

 d = dir([path,'\*.mat']);
 for i=1:length(d)
    dirname = d(i).name(1:end-4);
    mkdir(dirname);
    copyfile(d(i).name,[path,'\',dirname,'\',d(i).name(1:end-4),'.BEH.mat']);
 end
 