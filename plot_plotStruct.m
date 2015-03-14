function plot_plotStruct(hax, s, varargin);
%
%   plots plotStruct
%
% Opts: 
%      -GroupOrder:  specifies that following input is a cell array of group_id in order of plot
%      -PlotLine:   adds line between plots
%
%   Aaron Gruber;   2011_09_23

set(0, 'DefaulttextInterpreter', 'none')
if(isempty(hax))
    figure
    hax=axes;
end
axes(hax)


pltSymStr = '.';  % default to no line
pltMarkerSeries = {'o','+','<','d','s'};
numGp1Var = length(s.gp);
pltOrder = 1:numGp1Var;     % vector of indicies for plotting order
gpOrder = [];
Xorder = 1:numel(s.xLabels);
Legend_flg = 1;             % flag for indicating legend

varargin_txt = varargin;        
varargin_txt(~cellfun(@ischar, varargin)) = {'placeholder'};

if nargin > 1; 
   optIndx = find(cellfun(@isempty, strfind(varargin_txt,'-'))==0);
   opts_cell = varargin_txt(optIndx);
   for opt=opts_cell
       switch opt{:}
           case '-GroupOrder'
               indx = find(strcmp('-GroupOrder',varargin_txt)==1);
               gpOrder = varargin{indx+1};
           case '-PlotLine'
               pltSymStr = ['-'];    % adds line between points
           case '-NoLegend'
               Legend_flg = 0;
           case '-Xorder'  % specify the order for x points
               indx = find(strcmp('-Xorder',varargin_txt)==1);
               XorderIn = varargin{indx+1};
               if(numel(XorderIn)==numel(Xorder))
                   Xorder = XorderIn;
               else
                   disp('wrong number of elements for -Xorder; ignoring this command');
               end
           otherwise
               if(strcmp(opt{:}(1),'-'))
                   error([opt{:}, ' is not a valid option']);
               end
       end 
   end
end
 
%% this section is for re-ordering the plot order via option 'GroupOrder'
labelError = 0;
if(length(gpOrder)==numGp1Var)
    for k=1:numGp1Var
       indx = strmatch(gpOrder{k},s.cLabelVals); % find order 
       if(indx > 0)
           pltOrder(k) = indx;
       else
           disp(['GroupOrder item: ''',gpOrder{k},''' does not match grouping labels in plot structure: ignoring ordering command '])
           labelError = 1;
       end
    end
else
   disp('list of groups specified in -GroupOrder do not match group names in plot structure, ignoring the ordering command'); 
end
if(labelError)
    pltOrder = 1:numGp1Var; % if any errors during re-ordering, then ignor and recreate default ordering vector
end


    for i=1:numel(pltOrder)
       k=pltOrder(i);
       [numXpts,numtrace] = size(s.gp(k).dat);
       xincr = 0.15*(exp(-(numGp1Var)/6));   % just a func to control spacing between points betwen groups; ensures that all points fit between +-0.2 %0.15; %(numGp1Var-1)/20 - (numGp1Var-1)/2
       xdelt = -(numGp1Var-1)/2*xincr + (k-1)*xincr;
       %errorbar(hax, s.xdat*ones(1,numtrace)+xdelt, s.gp(k).dat, s.gp(k).var, pltSymStr,'Marker',pltMarkerSeries{mod(i-1,numel(pltMarkerSeries))+1}, 'color', getColor(i,'num'),'MarkerFaceColor',getColorFace(i,'num'));
       errorbar(hax, s.xdat*ones(1,numtrace)+xdelt, s.gp(k).dat(Xorder), s.gp(k).var(Xorder), pltSymStr,'Marker',pltMarkerSeries{mod(i-1,numel(pltMarkerSeries))+1}, 'color', getColor(i,'num'),'MarkerFaceColor',getColorFace(i,'num'));
       hold on
    end
    ylabel([s.yLabel, ';  ',s.ctMeth,' +/- ',s.varMeth])
    xlabel(s.xLabel)
    
    set(gca,'XTick',s.xdat,'XTickLabel',s.xLabels(Xorder))
    %set(gca,'XTick',s.xdat,'XTickLabel',s.xLabels)
    if(Legend_flg); legend(s.cLabelVals{pltOrder}); end
    