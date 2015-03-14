function [ses, varargout] = loadBEH(paths, varargin)
%
%  ses = loadBEH(paths, mode, extention)
%
%  Load BEH behavioral files.
%
%  INPUTS:
%       
%       paths:  filepath(s) of parent direcories. String or cell array of
%               strings if more than one path specified
%       
%       options:
%            -r: recursive mode looks in all folders
%            -ext: specifies extension to look for e.g. 'BEHv1.mat'. Default
%               is '.mat' if not otherwise specified. Format is 'switch',
%               'value'.
%            -fpattern: specifies any part of file to match
%            -dpattern: specifies any part of the directory to match
%            -registry: an AnimalRegistry object - for accounting (optional)
%
%       If no inputs; attempts to load mat files in pwd. Check that they
%       are BEH files. Registry is created.
%
%   OUTPUTS:
%       ses:  array of session objects
%
%   EXAMPLE:
%       ses = loadBEH({'C:\Data\Folder1'; 'C:\Data\Folder2'}, '-ext', 'BEHv1.mat'}
%           will load all files with extension 'BEHv1.mat' in Folder1 and
%           Folder2.
%
%       ses = loadBEH('C:\Data\Folder1', '-r', '.mat'}
%           will load all *.mat files in all subfolders of Folder1.
%
%       ses = loadBEH(pwd,'-r','-dpattern','Gr','fpattern','_14');
%           will load all *mat files recursively under pwd with directories
%           containing 'Gr' and filenames containint '_14'
%
%  Aaron Gruber; 2011_9_10
%
% TODO:  
%   - set option for filtering (i.e. loading) based on some tag in the BEH
%   file (e.g. 'Date', 'Rat_id', 'type',...)
%   - maybe break into find files & load files?
%

% default vals
recurse = false;
extens = '.mat';
fpat = '';
dpat = '';
registry = [];

% parse options 
if nargin < 1;
   paths = pwd;     % if no inputs, load all BEH files in pwd
end
if nargin < 2;
   registry = AnimalRegistry; % make a registry - needed for session constructor
end

if(isempty(paths))
   paths = pwd; 
end

% because animal registry is not char, string functions crap out. So make a
% copy with placeholders for non-char types
varargin_txt = varargin;        
varargin_txt(~cellfun(@ischar, varargin)) = {'placeholder'};

if nargin > 1; 
   optIndx = find(cellfun(@isempty, strfind(varargin_txt,'-'))==0);
   opts_cell = varargin_txt(optIndx);
   for opt=opts_cell
       switch opt{:}
           case '-r'
               recurse=true;
           case '-ext'
               indx = find(strcmp('-ext',varargin_txt)==1);
               extens = varargin{indx+1};
           case '-fpattern'
               indx = find(strcmp('-fpattern',varargin_txt)==1);
               fpat = ['*',varargin{indx+1}];
           case '-dpattern'
               indx = find(strcmp('-dpattern',varargin_txt)==1);
               dpat = ['*',varargin{indx+1},'*\'];
           case '-registry'
               indx = find(strcmp('-registry',varargin_txt)==1);
               registry = varargin{indx+1};
           otherwise
               if(strcmp(opt{:}(1),'-'))
                   error([opt{:}, ' is not a valid option']);
               end
       end 
   end
end

if(~iscell(paths))
    paths = {paths};
end
if(isempty(strmatch(class(registry),'AnimalRegistry'))) % test if registy is proper class
    registry = AnimalRegistry;
    varargout(1) = registry;                    % if we made it, send back out
end

% find files - recursive looks in
reshape(paths,numel(paths),1); % ensure nx1 for easy visualization and so transpose (next line) works
for pn=paths'
    pn = char(pn);
    if(pn(end)~='\')
        pn(end+1) = '\';       % make sure last char is '\'
    end
    if(recurse)
        d = rdir([pn,'\**\',dpat,fpat,'*',extens]); 
    else
        d = rdir([pn,fpat,'*',extens]);
    end
end

if(isempty(d))
   disp('no files found'); 
   ses = [];
else
    numLoaded = 0;
    for i=1:length(d)                    % load sessions
       msg = ['loading BEH file', d(i).name, ' : ',num2str(i,'%03d'),'/',num2str(length(d),'%03d')];
       disp(msg)
       f=whos('-FILE',d(i).name); % get file info (rather than load whole file)
       %if(strmatch(f.name,'session'))                    % if struct name is session, load it
       if(strncmp('session', f.name,7)) 
           numLoaded = numLoaded+1;
           ses(numLoaded) = Session(d(i).name,registry); 
       end
       fprintf(1,repmat('\b',1,numel(msg)+1));
    end
end