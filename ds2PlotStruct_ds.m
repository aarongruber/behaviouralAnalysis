function s = ds2PlotStruct_ds(ds, groups, datavar, varargin)

% if only one group tag - make it the second; make fisrt one empty

% if groups{1} is empty - add column to ds with same dummy tag (all)
if(isempty(groups{1}))
    [len,~]=size(ds);
    TmpName = 'all'
    ds = [ds, dataset({nominal(ones(len,1)),TmpName})];  % add common data to all observations
    groups(1) = {TmpName};
end

if(numel(groups)==3)
    if(isempty(groups{3}))   % ignore 3rd grouping var if empty
        groups = groups([1,2]);
    end
end

ctMeth = 'mean';  % method for central tendancy 'mean', 'median', ...
varMeth = 'sem';  % method for variation: 'std', 'stderr',  

% Options:
%  -plot


% Group var is parsed as: x color, x var, mult trace var

% if any group var is empty [] or '', make a new colum that has same tag

if(numel(groups)>3)
   error('this function can not parse more than 3 dimensions') 
end

%% find if we need to average by testing if all group tags are unique
numGp = numel(groups);
id = '';
for k=1:numGp
    id = strcat(id,char(ds.(groups{k})));
    
    %ds.(groups{k}) = nominal(ds.(groups{k}));  % ensure that labels are nominal - mostly for correct labels
end

% lazy way - this will average if non-unique elements are present, and will add column 
% for variance in either case. If code is slow, use logic to determine if
% average must be made, and if not add the column of varianc=0 manually.

ds=grpstats(ds,groups,{varMeth,ctMeth},'datavars',datavar); 
vardatavar = [varMeth,'_',datavar]; % grpstats will make this column name
datavar = [ctMeth,'_',datavar];

% if(numel(id) ~= numel(unique(id)));
%     % need to average
%     ds=grpstats(ds,groups,{varMeth,ctMeth},'datavars',datavar);
%     % % dataset no longer has original datavar name -> has prefix of stat method
%     % % Change it back for consistency with code below;
%     %VarNames = get(ds,'VarNames');
%     %VarNidx = ~cellfun(@isempty, strfind(VarNames,[ctMeth,'_',datavar]));
%     %VarNames(VarNidx) = {datavar};  % change var name back to original
%     %set(ds,'VarNames',VarNames); %
%     % % or just change datavar name
%     datavar = [ctMeth,'_',datavar];
% else
%     
% end

%% find index for each element
ds=sortrows(ds,groups(1)); % sort by color group (1st grouping tag)
gpi = nan(length(ds.(groups{k})),numGp);  % preallocate to NaN (unassigned data will not be plotted)

for k=1:min(2,numGp)
    gpi(:,k) =  grp2idx(ds.(groups{k}));
end
% For the 3rd grouping var (usually data for every animal in group), we
% don't want grouping index over all dataset - just within color group (1st
% group). Therefore assign 3rd grouping index within each group (1st) 
if(numGp==3)
    for k=1:max(gpi(:,1))           % loop over all levels of 1st grouping var
        idx=gpi(:,1)==k;
        gpi(idx,3) = grp2idx(char(ds.(groups{3})(idx))); % have to cast to char, otherwise it knows about the other members and reproduces the behavior we want to avoid (rat ID from one group are not in the other and result in columns of nan)
    end
end
C = num2cell(gpi);   % cast to cell array b/c can extract comma list

[numDat,~] = size(gpi);
if(numGp>1)
    m = nan(max(gpi)); % preallocate for speed only if more than one grouping variable
    v = nan(max(gpi));
end
% loop through all observations in ds and assign to position in matrix
for k=1:numDat
    m(C{k,:}) = ds.(datavar)(k);  % could probably assign using vector index, but RHS will have to be reshaped to match (error prone?)
    v(C{k,:}) = ds.(vardatavar)(k);
end

%% make a structure for plotting
numGp1Var = max(gpi(:,1));
for k=1:numGp1Var
    if(numGp==3)
        s.gp(k).dat = squeeze(m(k,:,:));
        s.gp(k).var = squeeze(v(k,:,:));
    elseif(numGp==2)
       s.gp(k).dat = m(k,:)';
       s.gp(k).var = v(k,:)';
    else
       s.gp(k).dat = m(k);
       s.gp(k).var = v(k);
   end
   
end
s.varMeth = varMeth;
s.ctMeth = ctMeth;
s.datavars = datavar;
s.yLabel = datavar;
s.cLabel = groups{1};
s.cLabelVals = getlabels(nominal(ds.(groups{1})));

if(numGp>1)
    s.xLabel = groups{2};
    s.xdat = [1:max(gpi(:,2))]';   
    s.xLabels = getlabels(nominal(ds.(groups{2})));  
    if(numGp==3)
        s.yLabel = [s.yLabel, ' [for each ',groups{3},']'];
    end
else
    s.xLabel = groups{1};
    s.xLabels = {''};
    s.xdat = [1];
end

if(0)
    set(0, 'DefaulttextInterpreter', 'none')
    figure
    for k=1:numGp1Var
       [numXpts,numtrace] = size(s.gp(k).dat);
       xincr = 0.15*(exp(-(numGp1Var)/6));   % just a func to control spacing between points betwen groups; ensures that all points fit between +-0.2 %0.15; %(numGp1Var-1)/20 - (numGp1Var-1)/2
       xdelt = -(numGp1Var-1)/2*xincr + (k-1)*xincr;
       errorbar(s.xdat*ones(1,numtrace)+xdelt, s.gp(k).dat, s.gp(k).var,'o', 'color', getColor(k,'num'),'MarkerFaceColor',getColor(k,'num'));
       hold on
    end
    ylabel([s.yLabel, ';  ',s.ctMeth,' +/- ',s.varMeth])
    xlabel(s.xLabel)
    legend(s.cLabelVals{:})
    set(gca,'XTick',s.xdat,'XTickLabel',s.xLabels)
end



%%  notes %%
% could instead 1: separate into datasets by group1 tag;
%               2: make 2d matricies of other tags
%               3: 
% plot - offset on x axis

