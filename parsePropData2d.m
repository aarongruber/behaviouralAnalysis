function d = parsePropData2d(TrCell, gpCrit, propName, AvgLevel, TrLabels)
% groups data into 2d cell array for plotting and stats
%    AvgLevel is one 'rat', 'session', 'all', 'none'
%
% input CellArray of trials, grouping tags, property name, collapseLevel {animal, session, trial} 
%
%
% add average by animal and group
%
% add check if session data - don't average for trials


if(isprop(TrCell{1}(1).Session,propName) && ~strcmpi(AvgLevel,'session'))
    AvgLevel = 'session';
    disp(['Property: ', propName,' is a session-level property so only one per session is used' ])
end

if(length(TrCell)>1 && length(gpCrit) ==1)
    for k=1:length(TrCell)
        TrCell{k}.setTempGroupingTag(k);    % set grouping tag for different groups 
    end
    critTag = {'TempGroupingTag',gpCrit{:}};        % first criteria is grouping tag
elseif( length(TrCell)==1 && length(gpCrit)==2 )
    critTag = gpCrit;
else
    error('function parsePropData2d only works for 2 dimensions or less')
end

allTr = [TrCell{:}];
for k=1:length(critTag)
    allVals = allTr.getPropVals(critTag{k});
    if(~isempty(allVals))
        if( isnumeric(allVals{1}))
            allVals = cellfun(@(x) {num2str(x)},allTr.getPropVals(critTag{k}));
        end
        crit{k} = unique(allVals);
    else
        crit{k} = '';
    end
end


for k=1:length(crit{1})
   for m=1:length(crit{2})
      disp(['parsing ',num2str(length(crit{2})*(k-1)+m),' / ', num2str(length(crit{1})*length(crit{2}))])
      d.label{k,m} = {crit{1}{k}, crit{2}{m}};
      trials = allTr.getTrials(critTag{1}, crit{1}{k}) & allTr.getTrials(critTag{2}, crit{2}{m});
      switch AvgLevel
          case 'session'                                % make session averages
              sesCellArray = trials.parse1d('Ses_id'); % cell array parsed into sessions
              for n=1:length(sesCellArray)
                  sesVals_m(n)= mean(cell2mat(sesCellArray{n}.getPropVals(propName)));
                  sesVals_median(n)= median(cell2mat(sesCellArray{n}.getPropVals(propName)));
              end
              d.dat{k,m} = sesVals_m;
              d.dat_mean{k,m} = mean(sesVals_m);
              d.dat_median_m{k,m} = mean(sesVals_median);
              d.dat_std{k,m} = std(sesVals_m);
         case 'rat'                                % make session averages
              ratCellArray = trials.parse1d('Rat_id'); % cell array parsed into sessions
              for n=1:length(ratCellArray)
                  ratVals_m(n)= mean(cell2mat(ratCellArray{n}.getPropVals(propName)));
                  ratVals_median(n)= median(cell2mat(ratCellArray{n}.getPropVals(propName)));
              end
              d.dat{k,m} = ratVals_m;
              d.dat_mean{k,m} = mean(ratVals_m);
              d.dat_median_m{k,m} = mean(ratVals_median);
              d.dat_std{k,m} = std(ratVals_m);
          case 'all'
              sesVals= trials.getPropVals(propName);
              d.dat{k,m} = cell2mat(sesVals);
              d.dat_mean{k,m} = mean(cell2mat(sesVals));
              d.dat_median{k,m} = median(cell2mat(sesVals));
              d.dat_std{k,m} = std(cell2mat(sesVals));
          case 'none'
              d.dat{k,m} = cell2mat(trials.getPropVals(propName));
          otherwise
      end
   end
end
allTr.setTempGroupingTag();  % clear temp tags
d.dimLabel = critTag;
d.dimTickLabel = crit;

