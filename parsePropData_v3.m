function d = parsePropData_v3(TrCell, gpCrit, propName, AvgLevel, TrLabels)
% groups data into 2d cell array for plotting and stats
%
%
% input             CellArray of trials, 
%                   grouping tags{'Rat_id','Date','Group',...}, 
%                   property name, collapseLevel {animal, session, trial} 
%                   level of averaging {n=none (all data returned),ses=session, all = across all trials}
%                   TrLabels: labels for TrCell (only if dim >1)
%

% if(lenght(TrCell) + length(gpTagCell) <3)
%     error('function parsePropData2d only works for 3 dimensions or less')
% end
% 
% if(length(TrCell)>1 && length(gpTagCell) ==1)
%     for k=1:length(TrCell)
%         TrCell{k}.setTempGroupingTag(k);    % set grouping tag for different groups 
%     end
%     criteria{1} = 'TempGroupingTag';        % first criteria is grouping tag
%     criteria{:} = [criteria, gpTagCell{:}]; % set remainder of criteria
% elseif( length(TrCell)==1 && lentgh(gpTagCell)==2 )
%     
% end

for k=1:length(gpTagCell)
    crit{k} = unique(([TrCell{:}]).getPropVals(gpCrit{k})); % cell array of criteria for sorting
end

for k=1:length(crit)
    
end

for k=1:length(TrCell)
    
    
end

% check that length of all grouping dimensions is <3: length(TrCell)+length(gpTagCell)
% for all TrCell objects
% set TempGroupingTag of TrCell
% if length TrCell >1, set TempGrouping as first grouping criteria
    % set TrialLabels as next criteria
% else if trCell ==1 && length(gpTagCell) ==2;
    % set gpTagCell as the grouping criteria

% make a table of grouping criteria
    % find unique values (in all TrCell) for each criteria :: convert cell
    % array of objects to concatonated list:
    % 'criteria{k} = unique([TrCell{:}].getPropVals(criteria{k}))'
    
    % for loop to make every permutation
    
% loop over criteria list
    % loop over criteria list 2
        % find intersection: trials = [TrCell{:}].getTrials(criteria{k}) & [TrCell{:}].getTrials(criteria{m})
        % if(sesAvgFlag)
            % sesCellArray = trials.parseIntoSessions();
            % for m=1:length(sesCellArray)
                  %tmp2(m)= (sesCellArray{m}.getPropVal(propName));
            % end
            % d.dat_m = mean(tmp2);
            % d.dat_median = meedian(tmp2);
            % d.dat_std = stdev(tmp2);
        %else
            % d.dat_m = mean(trials.getPropVal(propName));
            % d.dat_std = stdev(trials.getPropVal(propName));
            % d.dat_median = median(trials.getPropVal(propName));
%       end
  
for k=1:length(TrCell)
    
    
end

% build groupBy({group, animal, session, trial})
% groupingTag & method for setting & clearing