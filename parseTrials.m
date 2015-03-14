function [trials txt] = parseTrials(TrCell, levels, TrialLabels)
%
% INPUTS:
%   TrCell: Cell array of Trial objects (n=1 to inf)
%   levels: Cell array of code names for analysis {Group, Animal,
%   Date, Session, Trial} or custom property

% output: plotDataObject, data

% to do: add input for object tags, 
%        1)add field to return data instead of trials objects
%        2)output field for category ids and labels
%        3) output a field of tags for each one (for anova)
%       ** figuure out how to return trials - maybe just transpose rather
%       than trying to find...
%    ** fix searching by group 
%       
% for each field, concatonate values for each object member
% find unique values as categories across all input objects
% find matches (intersect) to full list for each input object


    
for k=1:length(levels)
    category{k} = [];    % Unique categories for sorting data (concatonated across input Trial objects so output is alligned)
    for i=1:length(TrCell)
        switch levels{k}

            case 'Group'
                txt.label{k} = 'Group';   
                vals{k,i} = TrCell{i}.getSessionInfo('Group_id');
            case 'Animal'
                txt.label{k} = 'Animal';
                vals{k,i} = TrCell{i}.getSessionInfo('Rat_id');  
            case 'Session'
                txt.label{k} = 'Session';   
                vals{k,i} = TrCell{i}.getSessionInfo('Ses_id');
            case 'Date'
                txt.label{k} = 'Date'; 
                vals{k,i} = TrCell{i}.getSessionInfo('Date');
            case 'Trial'
                txt.label{k} = 'Trial';        
            otherwise
                txt.label{k} = levels{k};   
                vals{k,i} = TrCell{i}.getSessionInfo(char(levels{k}));  % this could be the only line needed -> pass in 'Rat_id'... instead of Animal
        end
    end
    category{k} = unique([category{k}, vals{k,:}]);
end

if(length(TrCell)>1 && length(levels)==1)
    txt.x = category{:};
    txt.y = TrialLabels;
    txt.label{k+1} = 'obj';
    for i=1:length(TrCell)
       for k=1:length(levels)
           for m = 1:length(category{k}) % for each category
                indx = strfind(vals{k,i},char(category{k}(m))); 
                trials{i,m} = TrCell{i}([cellfun(@(x) ~isempty(x), indx)]);
           end
       end
    end
elseif(length(TrCell)==1 && length(levels)==2) % nope this has to include combinatorial search either in two steps or as one: best approach is probalby to move search to Trials
%     txt.x = category{1};
%     txt.y = category{2};
%     for i=1:length(TrCell)
%        for k=1:length(category{2})
%            for m = 1:length(category{1}) % for each category
%                 indx = strfind(vals{k,i},char(category{k}(m))); 
%                 trials{k,m} = TrCell{i}([cellfun(@(x) ~isempty(x), indx)]);
%            end
%        end
%     end
else
    error('This combination of Trials objects and levels not yet implemented')

end
'hi'
