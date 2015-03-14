function ds = makeDS(propName, TrCell, gpFields, avgLevel, trialGpNames, varargin)
%
%  ds = makeDS(propName, TrCell, gpFields, avgLevel, trialGpNames, '-option', 'optParm')
%
% Create a dataset array with data reqired for
% plotting and stats, and to pre-process data by taking average at
% specified level (Ses_id, Rat_id, ...).
%
% propName = property on which to do stats; if it is a session-level
%              property, then averaging is done by session automatically
%
% TrCell = cell array of one or more groups of class {trials}
%
% gpFields = properties over which to do grouping {'Rat_id', 'Date',...}
%
% avgLevel = properties for averaging before stats (e.g. average over trials or animals before doing
%           states){'Ses_id', 'Rat_id',...}; usually this is {'Ses_id'}, or
%           {'Ses_id','Rat_id} - the subsequent plotting will compute the variance at
%           the next highest level. The averaging will be done in the order
%           given e.g. {'Ses_id','Rat_id} will be session for each rat, then rat
%
% trialGpNames = names if more than one group of trials passed in via 'TrCell'
%
%
% Optional inputs
%
%   '-StatType' <{'mean'} (default),{'median'}, {'sem'}, {'meanci'}, .. or other supported by 
%               grpstats.m, or custom function as {@(x) f(x)} (in cell array) > 
%               specify a stat (mean is default).  
%               example custom function: {@(x) quantile(x,.2)} is 20th percentile
%
%
% example: 
%   ds = makeDS2('numTrials', Data_Trials, {'Group_id','Mode' }, {'Ses_id','Rat_id'});
%
%   where ds is a dataset of 'numTrials' parsed by 'Group_id' and
%   'Mode', after being averaged over sesssion and then for each rat.
%
%
%  Aaron Gruber,   
%                  2015_02_5   added -StatType as an option

%
%
%
% TO DO: add switch -stat (e.g. 'mean', 'medial') stat specifies which
% method to compute central tendency if pre-averaging is specified
%

%% defaults
statType = {'mean'}; % by default, compute median; unless optionally changed

%% process optional imputs
if nargin > 5; 
    varargin_txt = varargin;        
    varargin_txt(~cellfun(@ischar, varargin)) = {'placeholder'};
   optIndx = find(cellfun(@isempty, strfind(varargin_txt,'-'))==0);
   opts_cell = varargin_txt(optIndx);
   for opt=opts_cell
       switch opt{:}
           case '-StatType'
               indx = find(strcmp('-StatType',varargin_txt)==1);
               statType = varargin{indx+1};
            otherwise
               if(strcmp(opt{:}(1),'-'))
                   error([opt{:}, ' is not a valid option']);
               end
       end 
   end
end
%%

% some formatting stuff
if(~iscell(TrCell)); 
    TrCell = {TrCell};
end

if(~exist('avgLevel'))
    avgLevel = {'none'};
elseif(isempty(avgLevel))
    avgLevel = {'none'};
end

if(~iscell(avgLevel)); 
    avgLevel = {avgLevel};
end

% check if propName is a session level property
if(isprop(TrCell{1}(1).Session,propName)); 
    avgLevel = unique([{'Ses_id'},avgLevel], 'stable'); %- Ses_id needs to go first in the list - use 'stable' to prevent reordering
    disp(['Property: ', propName,' is a session-level property so only one per session is used' ])
end

% add property for the level of averaging
for i=avgLevel
    if(ispropInAnyLevel(TrCell{1}(1), i{1}))
        gpFields=union(gpFields,i);
    elseif(~strcmpi(i{1},'none'))
    %else
        error([i{1}, ' is not a valid property of Trials']) 
    end
end

% handle multiple Trial groups passed in - set this as frist grouping parameter which
% will be coded in color in plots 
if(numel(TrCell)>1)
    if(~exist('trialGpNames')); trialGpNames = ''; end
    if(isempty(trialGpNames))
        disp(['Warning: Multiple trials groups passed in, but no group names passed into makeDS: generic ones are used'])
        trGp =  strcat({'gp'},num2str((1:numel(TrCell))','%-d'));  % generic names
    else
        trGp = trialGpNames;
    end
    for k=1:length(TrCell)
        TrCell{k}.setTempGroupingTag(trGp{k});      % set grouping tag for different groups 
    end
    gpFields = {'TempGroupingTag',gpFields{:}};      % first criteria is grouping tag
end

allTr = [TrCell{:}]; % concatonate trials

%build array
for k=1:numel(gpFields)
    d{:,k} = nominal(allTr.getPropVals(gpFields{k}))';
end

% get the data for the property
dat = [allTr.getPropVals(propName)]'; % this is cell array - convert to numeric in cell array
%dat(cellfun('isempty',dat)) = {NaN}; % set empty to NaN
%dat = [dat{:}]';

ds = dataset(d{:},dat);    % if dataset exsists:   % ds.(propName) = cat(1,dat{:});
ds=set(ds,'VarNames',{gpFields{:}, propName});   
ds = set(ds, 'UserData', avgLevel);

% loop over gpFields, taking them out of the grouping param for each
% sucsessive iteration - will have an average for the last one
for i=1:numel(setdiff(avgLevel,'none'))+1
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% changed here %%%%%%%%%%%%%%%%%%%%%%%%   
    %ds=grpstats(ds, gpFields,{'mean'},'datavars',propName); % make average
    ds=grpstats(ds, gpFields, statType,'datavars',propName); % make stats, indicated by cell statType
    
    if(~strcmpi(statType,'mean')) % give a shout if not mean, just to prevent errant use.
        disp(['in MakeDS, statType is ', statType])
    end
    %disp('now doing median in makeDS');    statType
    %ds=grpstats(ds, gpFields,{'median'},'datavars',propName); % make average
    %ds=grpstats(ds, gpFields,{@(x) quantile(x,.2)},'datavars',propName); %do 20% percentile of distribution
    
   % dsB=grpstats(ds, gpFields,{'mean','@nanmean'},'datavars',propName); % make average
    % rename fields - probably an eaiser way to do this
    fn=get(ds, 'VarNames'); % reset name of variable that was averaged
    for k=1:numel(fn)
       if(~isempty(strfind(fn{k},propName)))
           fn{k} = propName;
           ds=set(ds,'VarNames',fn);
       end
    end
    % restrict the names in gp Fields to make new average (drop one name
    % each iteration, cuasing averaging over that level on next interation
    if i<numel(avgLevel)+1 % prevent err on last iteration
        gpFields = setdiff(gpFields,avgLevel(1:i));
    end
    % generate warning if no averaging was done on a particular iteration
    % (e.g. if 'Rat_id' is given prior to 'Ses_id' in the avgLevel list)
    if(numel(avgLevel)>0 && sum(ds.GroupCount) < numel(ds.GroupCount))
        warning(['No averaging done in dataset for property: ',propName])
    end
end

allTr.setTempGroupingTag();  % clear temp tags
