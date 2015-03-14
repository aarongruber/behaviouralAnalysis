classdef Trial < dynamicprops % handle 
    properties (SetAccess = private)
        
        Session = Session.empty;
       
        % Trial Data;  first letter is caplitalized whether or not in data
        
        Odors = nan;
        Trial_start_ts = nan;
        L_light_on_ts = nan;
        L_light_off_ts = nan;
        R_light_on_ts = nan;
        R_light_off_ts = nan;
        Odor_port_on_ts = nan;
        Odor_port_off_ts = nan;
        Odor1_on_ts = nan;
        Odor1_off_ts = nan;
        Odor2_on_ts = nan;
        Odor2_off_ts = nan;
        Tone_on_ts = nan;
        Tone_off_ts = nan;
        Well_on_ts = nan;
        Well_off_ts = nan;
        Well_id = nan;
        L_reward_ts = nan;
        R_reward_ts = nan;
        L_lick_ts = nan;
        R_lick_ts = nan;
        Error_id = nan;
        Error_string = nan;
        Error_flg = nan;
        Odor_port_timeout_flg = nan;
        Odor_port_early_go_flg = nan;
        Wrong_choice_flg = nan;
        Well_early_go_flg = nan;
        Stimulation_on_ts = nan;
        Stimulation_off_ts = nan;
        House_light_on_ts = nan;
        House_light_off_ts = nan;
        Odor1_id = nan;
        Unique_trial_id = nan;  % unique identifyer
        Trial_indx = nan;        % trial index
        %Params = struct();
        Mode = nan;                   
        Odor2_id = nan;
        Tone_period = nan;
        Rewarded_well = nan;
        L_reward_id = nan;
        L_reward_num = nan;
        L_reward_delay = nan;  % this should be 'L_reward_delay'l; it was 'L_bolus_delay'
        R_reward_id = nan;
        R_reward_delay = nan;
        R_reward_num = nan;
        Debounce_up = nan;
        Debounce_down_odor_port = nan;
        Debounce_down_well = nan;
        Date = nan;   % copy here for speed
        Rat_id = nan; % copy here for speed
        Ses_id = nan;
        Group_id = nan; % for speed; due to frequency of re-grouping, this one makes a big difference
        Drug = nan;
        
        Well_vacuum_on_ts = nan;  % for reversal learning
        Well_vacuum_off_ts = nan;
        Odor_port_all_on_ts = nan;
        Odor_port_all_off_ts = nan;
        R_well_all_on_ts = nan;
        R_well_all_off_ts = nan;
        L_well_all_on_ts = nan;
        L_well_all_off_ts = nan;
        
        GoNoGoType = nan;      % indicator for go/no-go discrimination
        OkGo_flg = nan;            % correct go
        ErrGo_flg = nan;           % err go
        TempGroupingTag =nan;  % for analysis
        Block_id = nan;
    end
    properties (Dependent)
        Group
    end
    properties
       Data = nan;              % for holding data 
       NumPreTrialPokes = nan;  % number odor pokes prior to trial start
       NumPreTrialRewardWellPokes = nan' % entries into reward reward wells prior to trial start
       InitTime = nan;          % time to initiate trial
       ReactTime = nan;
       ResponseTime = nan;
       LicksBeforeRew = nan;
       LicksTotal = nan;
       NumRewards = nan; % number of reward on a trial
       %RespEntropy = nan;
       %ProbWSLS = nan;  % probability of win stay lose switch strategy
       %ProbSame = nan;  % probability of repeating choice
       %LogitRegB = nan; % logistical regression coefficients
       %LogitRegFitAcc = nan; % logistical regression accuracy
       HighValWell_id = nan;   % well with higher value
       ChoseHighVal_flg = nan;  % in
       Block_trial_indx = nan;  % index within block
       Block_number = nan;      % order of block in session
       Grouped_trial_indx = nan; % trial index in groups of some integer (e.g. 1 1 1 2 2 2 3 3 3....) for stats & plotting
       TestDay_idx = nan;   % sequential ordering of test day - for comparison across groups tested at different times *** must take care b/c missing day will offset all subsequent data
       ConsecTrialAcrossSess_indx = nan; % consecutive trials across sessions
       TrialInSession_indx = nan;  % trial number
       ElapsedTime = nan;     % total elapsed time across all sessions at a given trial (seconds)(for finding total time to some criterion)
       SessionTime = nan; % time since beginning of session (seconds)
       SessionCountdownTimeUntillSessionEnd = nan; % time until session end (seconds)
       UncountedGO = nan; % number of uncounted trials: poke -> well before lights come on; ** new **
       PokeBeforeLight = nan; % nose pokes before light came on (the show up in the previous trial)
       TimeInWell = nan; % elaped time from first entry to reward well until last exit  ** new **
       TotalPokes = nan; % total number of pokes on a trial
       InterTrialInterval = nan; % time from reinforcement (or omission of expected reinforcement) to next poke time - time of short memroy for reinforcement & choice
       SampledBothWellsBeforeTrial_flg = nan; % logical if rat checked both wells
       LastWellCheckedBeforeTrial_id = nan;   % ID of last well checked - e.g. if rat sampeld multiple wells
       IsSequentialFromGoodTrial_flg = nan; % logical for trial following good trial (e.g. no double sampling of wells)
       LoseSwitch_BasedOnOutOfTrialResp_flg = nan; % lose-switch computed based on last sampled well - e.g. if rat sampled both wells;   need to put it here to compute effect of delay
       LoseSwitch_NoMultWellResp_flg = nan; % lose switch only on trials without prior multiple well entries (e.g. to both sides)
       LoseSwitch_all_flg = nan;        % lose switch response using all trials
       WinStay_all_flg = nan;           % win stay over all trials
       WinStay_BasedOnOutOfTrialResp_flg = nan; % win stay using last well sampled (even if 2nd)
       WinStay_NoMultWellResp_flg = nan; % win stay flag for trials with only one well choice
       LastTrialRewarded_flg = nan;      % flag indicating if last trial was rewarded 
       ChoiceSwichFromLastTrial = nan;  % glag indicating that choice switched.

    end
    methods
        function group = get.Group(obj)
            if ~isempty(obj.Session)
                group = obj.Session.Animal.Group;
            else
                group = Group.empty;
            end
        end
    end
    methods
        function obj = Trial(session,tr)
            % Two ways to have input arguments:
            %       SESSION, TR
            %           SESSION is a Session object
            %           TR is a trial structure
            
            
            % If there are no input arguments, create the default Trial
            % object.  (This is needed to implement syntaxes such as
            % Trial.empty and initializing an array with "myTrials(10) =
            % Trial").
            
            % If there are input arguments...
            if nargin
                
                % Initialize a Trial array with the same length as the TR structure
                obj(length(tr)) = Trial;
                
                % If desired properties are in a structure, copy them up
                % a level to tr.* for algroithmic simplicity
%                 for fl = fieldnames(tr)'
%                     if(strcmp(fl,'params'))
%                        'hi' 
%                     end
%                    valididx = find(~cellfun(@isempty,{tr.(char(fl))})==1); % find non-empty ones
%                    if(~isempty(valididx))
%                        if(isstruct(tr(valididx(1)).(char(fl))))                % if it is a structure
%                                 %%%% the following is if the fieldnames of sub
%                                 %%%% structure is not the same from object to
%                                 %%%% object
%                                 %valididx = find(~cellfun(@isempty,{tr.(char(fl))})==1); % returns indicies of non empty params
%                                 %validtr = tr(valididx);
%                                 %allfields = cellfun(@(x) fieldnames(x)',{validtr.params},'UniformOutput',0);
%                                 %paramfields = unique([allfields{:}]);
%                                 %paramfields = intersect(lower(paramfields),lower(properties(obj)));
%                                 paramfields = intersect(lower(fieldnames(tr(valididx(1)).(char(fl)))),lower(properties(obj)));
%                            for fl_sub = paramfields'
%                                [tr.(char(fl))] = deal(tr.(char(fl)).(char(fl_sub))); % copy data up a level
%                            end
%                        end
%                    end
%                 end
                
              % Loop through the fieldnames on TR.
                % this code could be replaced with 'for fl = fieldnames(tr)' if you don't care about removing fields
                trfldnm =fieldnames(tr);  % get fieldnames of tr
                fl_list = trfldnm(ismember(lower(fieldnames(tr))',lower(fieldnames(obj)))); % list of fieldnames both in tr (data file) and in session.TR object.
                fl_list = [fl_list; trfldnm(cellfun(@(x) isstruct(tr(1).(char(x))), fieldnames(tr)))]; % add on structs (e.g. params) that are in tr
                
                missing_fl = setxor(lower(fl_list),lower(trfldnm)); % get any fieldnames that are in tr that are not in obj (e.g. those excluded in the ismember operation above);
                for fn = missing_fl'
                    %warning(['For file: ',[char(session.Rat_id),'-',char(session.Date)] ,';  property ''',fn{:},''' not read in b/c it is not listed as a Trial property'])
                end
                % conclusion of finding fieldnames in both tr data and Trial
                
                %for fl = fieldnames(tr)'   % this is superceeded by the above lines that pick out any fields not in 'Trial'
                for fl = fl_list'  
                    % the following code is overly complex b/c of having to
                    % copy data out of a substructure (params). If the
                    % files were first processed for a proper format (cap
                    % field names) and no substructures, this code could be
                    % reduced to only a couple of lines.
                    
                    % For each field name, deal the elements of TR to
                    % the properties of the object. Convert non-numeric to
                    % nominal to allow consistent use of equlity operators
                    % ( ==, <, ...)
                    objFName = char(fl);
                    objFName = [upper(objFName(1)), objFName(2:end)]; % capitalize first letter of property name
                    nanIndx = cellfun(@isempty, {tr.(char(fl))});   % find indicies of empty entries
                    %nanIndx = cellfun(@(x) numel(x)<1|isnan(x), {tr.(char(fl))},'UniformOutput',1);  % added isnan part 2013_2_21 because when all fields are empty, not working properly
                    if(~isstruct(tr(1).(char(fl)))) 
                        if(isnumeric([tr.(char(fl))]) || islogical([tr.(char(fl))]))
                            [tr(nanIndx).(char(fl))] = deal(nan); % convert empty entries to numeric nan
                            [obj.(objFName)] = deal(tr.(char(fl)));
                        else  
                            [tr(nanIndx).(char(fl))] = deal('nan'); % convert empty entries to numeric nan
                            fvals = num2cell(nominal({tr.(char(fl))}));
                            [obj.(objFName)] = fvals{:};
                            %[obj.(objFName)] = deal(nominal({tr.(char(fl))})); % cast as nominal
                        end
                    end
                    
                end
                
                % make a unique trial identifier for each trial
                IDnames = num2cell(nominal(strcat({[char(session.Rat_id),'-',char(session.Date),'-',char(session.Time),'-']},num2str((1:length(obj))','%04d'))));
                [obj.Unique_trial_id]=IDnames{:};
                TrIdx = num2cell([1:length(obj)]');
                [obj.Trial_indx] = TrIdx{:};
               % for k=1:length(obj)
               %     obj(k).Unique_trial_id = nominal([char(session.Rat_id),'-',char(session.Time),'-',num2str(k,'%04d')]); % if I knew how do replicate the string, could deal this rather than a loop
               %     obj(k).Trial_indx = k;
               % end
             
           if(isfield(tr,'params'))  % see if there is a structure of params
                % Handle the PARAMS field as a special case
                % Get the TRs with a non-empty PARAMS field
                %valididx = ~cellfun(@isempty,{tr.params});
                valididx = find(~cellfun(@isempty,{tr.params})==1); % returns indicies of non empty params
                validtr = tr(valididx);
                
                % Get the fieldnames of PARAMS that are also
                % properties of the class
                
                % The following works only if PARAMS always has the same
                % fieldnames:
                % paramfields = intersect(fieldnames([validtr.params]),properties(obj));
                
                % Since PARAMS may be different for each element of
                % validtr, use this:
                allfields = cellfun(@(x) fieldnames(x)',{validtr.params},'UniformOutput',0);
                paramfields_all = unique([allfields{:}]);
                paramfields = intersect(lower(paramfields_all),lower(properties(obj)));
                
                % For each of the properties of the class that
                % are in the PARAMS field, assign them to the
                % new object array
                
                %for pr = paramfields
                for ii = 1:numel(paramfields)  % the long loop for expression is broken in R2013b, so go back to a dumb way.
                    pr = paramfields{ii};
                    
                    pr = paramfields_all(strcmpi(pr,paramfields_all)); % take care of capitalization missmatch by finding original name (in tr.param)
                    fieldvalue = cellfun(@(x) x.(char(pr)),{validtr.params},'UniformOutput',0);

                     % convert fieldvalue from a cell array to a comma-separated list 
                    objFName = char(pr);
                    objFName = [upper(objFName(1)), objFName(2:end)]; 
                    
                    nanIndx = cellfun(@isempty, fieldvalue);
                    if(isnumeric([fieldvalue{:}])||islogical([fieldvalue{:}]))
                        [fieldvalue{nanIndx}] = deal(nan);
                        [obj(valididx).(objFName)] = fieldvalue{:};
                    else
                        [fieldvalue{nanIndx}] = deal('nan');
                        fieldvalue = num2cell(nominal(fieldvalue));
                        [obj(valididx).(objFName)] = fieldvalue{:};
                    end
                       
                    %[fieldvalue{cellfun(@isempty, fieldvalue)}] = deal;
                    %[obj(valididx).(objFName)] = deal(fieldvalue{:});                   
                    
                    %%%% or do it the the loop way
                    %for vidx = valididx
                    %    obj(vidx).(char(pr)) = fieldvalue{vidx};
                    %end
                end
                
           end % loop for params   
                
                % Create a link between the newly created object(s) and the Session
                [obj.Session] = deal(session);
                % copy this info from session level to each trial for speed
                [obj.Date] = deal(nominal(session.Date));
                [obj.Rat_id] = deal(nominal(session.Rat_id));
                [obj.Ses_id] = deal(session.Ses_id);
                [obj.Drug] = deal(session.Drug);
                
                obj = obj( cellfun(@(x) ~isnan(x), {obj.Trial_start_ts}) ); % get rid of any trials with empty start times
                
            end
        end
%%        
        function trials = and(obj1, obj2) % logical 'and' (intersection of sets): return trils that are in both objects
            if(isempty(obj1) || isempty(obj2))
                trials = [];
            else
                [~, obj1_indx, obj2_indx] = intersect([obj1.Unique_trial_id], [obj2.Unique_trial_id]);
                %[~, obj1_indx, obj2_indx] = intersect(cellfun(@(x) char(x), {obj1.Unique_trial_id}, 'UniformOutput',0),  ...
                %                                      cellfun(@(x) char(x), {obj2.Unique_trial_id}, 'UniformOutput',0));
            trials = obj1(obj1_indx);
            end
        end
%%        
        function trials = or(obj1, obj2) % logical 'or' (union of sets)
            [~, obj1_indx, obj2_indx] = union([obj1.Unique_trial_id], [obj2.Unique_trial_id]);
            trials = [obj1(obj1_indx),obj2(obj2_indx)];
        end
%%        
        function trials = ne(obj1, obj2) % logical ~= 'not equal to' (exclude from obj1 any that are in obj2)
            if(isempty(obj1))
                trials = [];
            elseif(isempty(obj2))
                trials = obj1;
            else
                [~, obj1_indx] = setdiff([obj1.Unique_trial_id], [obj2.Unique_trial_id]);
                trials = obj1(obj1_indx);
            end
        end
%%        
        function trials = plus(obj, val) % get trials in the future !!!! needs to be verified
            trials(length(obj)) = Trial;
            %trials = Trial;
            for k = 1:length(obj)
                trmax = length(obj(k).Session.Trials);  % total length of session
                ind = find([obj(k).Session.Trials.Trial_start_ts]==obj(k).Trial_start_ts); % find index in full trial structure: have to do it this way rather than trial id because trials may be filtered (remove no responses) so indicies may not be consecutive
                if(~isempty(ind))
                    if( ind+val >0 && ind+val <= trmax ) % prevent index out of bounds
                        trials(k) = obj(k).Session.Trials(ind+val);
                    end
                else
                    error('Trial does not exsist in parent session: were trials filtered (e.g. removed early go trials)?')
                end
            end
   %         trials = trials([cellfun(@(x) ~isempty(x), {trials.Trial_start_ts})]); % get rid of empty ones (e.g. out of bounds index that are not assigned)
            trials = trials([cellfun(@(x)  ~(isempty(x)|isnumeric(x)), {trials.Unique_trial_id})]);   % the isnumeric is to get rid of nan trials     
        end
%%        
        function trials = minus(obj, val) % get trials in the past
            trials = plus(obj, -val);
        end
%%        
        function setTempGroupingTag(obj, k) % set temporary tags for grouping purposes
            if(nargin==1)
                [obj.TempGroupingTag] = deal([]);
            else
                if(isnumeric(k))
                    [obj.TempGroupingTag] = deal(k);
                else
                    [obj.TempGroupingTag] = deal(nominal(k)); % if k is non-numeric make it nominal to be consistent
                end
            end
        end
%%        
        function setGroupID(obj, label)  % permenant grouping tag - this is a copy of the ses level info, but include at Trial level for speed
           [obj.Group_id] = deal(nominal(label));
        end
%%        
        function sessArray = parseIntoSessions(obj) % return cell array of trials grouped by session
                                                    % now replaced by parse1d
            sessArray = parse1d(obj,'Ses_id');
        end
%%        
        function sessArray = parse1d(obj, param) % return cell array of trials parsed by property {Ses_id, Rat_id, Date, ....}
             sessIDs = getPropVals(obj, param);         % get IDs
             uniqueIDs = getlabels(sessIDs);      % this is for nominals
             %uniqueIDs = unique(sessIDs);              % this is for strings
             for k=1:length(uniqueIDs)
                 indx = sessIDs==uniqueIDs{k};
                 sessArray(k) = {obj(indx)};
                 %indx = strfind(sessIDs,uniqueIDs{k});  % find all indicies with same session ID
                 %sessArray(k) = {obj(~[cellfun(@isempty, indx)])};
                 %sessArray(k) = {obj([cellfun(@(x) ~isempty(x), indx)])}; % turn into vector for indexing
             end
        end
%%
        
%     function sessArray = parse_cell(obj, param_cell)
%         
%     end
%     
%     function sessArray = parseNd(obj, param_cell) % return cell array of trials parsed by property {Ses_id, Rat_id, Date, ....}
%             %
%             % sessArray = parseNd({obj_cell}, {param_cell})
%             %
%             % parese input cell array of Trial objects by the tags in cell array
%             % 
%             % figure out how to imput obj as cell array!! - functional
%             % form?
%             %
%             % TODO: if not a property, try to feval the cell contents so
%             % as to allow expressions
%             %
%             % TODO **** allow more than one object in the input ****
%              if (numel(param_cell)>1)
%                  sessArray = parseNd(obj, param_cell(2:end)); % recursive call with one less parameter !!!
%              else
%                  for i=1:numel(obj) % loop over cell array 
%                     sessArray{i} = parse1d(obj, param_cell{1}); % parse each input **** will end up with ixn cell array
%                  end
%                  
%              end
%              sessArray = reshape(sessArray); %% have to reshape it b/c won't know how many levels there are
%              
%              
% %             k=0;
% %             indx_list = 1:numel(param_cell);
% %             for i=1:numel(param_cell)-1
% %                 TRarray = parse1d(obj, param_cell{i})
% % % %                      sessIDs = getPropVals(obj, param);         % get IDs
% % % %                      uniqueIDs = getlabels(sessIDs);      % this is for nominals
% % % %                      %uniqueIDs = unique(sessIDs);              % this is for strings
% % % %                      for k=1:length(uniqueIDs)
% % % %                          indx = sessIDs==uniqueIDs{k};
% % % %                          sessArray(k) = {obj(indx)};
% % % %                          %indx = strfind(sessIDs,uniqueIDs{k});  % find all indicies with same session ID
% % % %                          %sessArray(k) = {obj(~[cellfun(@isempty, indx)])};
% % % %                          %sessArray(k) = {obj([cellfun(@(x) ~isempty(x), indx)])}; % turn into vector for indexing
% % % %                      end
% %             end
%         end
        %%
        function trials = getTrials(obj,param,val)
           % returns trials in obj (Trials obj) for which param (in Trials, Params, or Session) is val
           % (val can be number, char or empty '')
           if nargin==1 || isempty(param)
               trials = obj;
           else
               allVals = getPropVals(obj,param);
               if(~isempty(val))  % if not empty, find matching values, otherwise return the empty ones (b/c strcmpi does not match '')
                   
                   % added the following on 2013_2_13 because trials with all nan for some param were causing 
                   %%%if(sum(isnan(allVals)) == numel(allVals)) % check if all of the values are nan - then a problem, as we don't know if to make it numeric or char
                   if(strcmpi(val,'nan')&sum(isnumeric([obj.(param)]))>0) % if the object is numeric (rather than nominal), then use isnan funct, otherwise can use ==
                     trials = obj(isnan([obj.(param)])); 
                   else
                     trials = obj([obj.(param)]==val);   
                   end
                   
                % trials = obj([obj.(param)]==val);  % this was the only line previous to 2013_2_21
                   %%% allVals = cellfun(@num2str,allVals,'UniformOutput',0); % convert to char
                   %%% trials = [obj([strcmpi(allVals, num2str(val))])]; % string comp is more robust since it will allow multiple digit vals or strings
               else
                   if(isnumeric([obj.(param)]) || islogical([obj.(param)]))
                        trials = obj([obj.(param)]==nan);
                   else
                       trials = obj([obj.(param)]=='nan');
                   end
                  %  allVals = [cellfun(@(x) isempty(x), allVals)];
                  %  trials = [obj(allVals)];
               end        
           end
        end
%%
        function allVals = getPropVals(obj,param)
           % returns values in obj (Trials obj) for property 'param' in Trials, Params, or Session levels 

           if nargin==1 || isempty(param)
               error('error in getPropVals:  must specify a property')
           else
               %allVals = cell(1,length(obj)); % preallocate to increase speed
               %%%% preallocate on nominals? or can make cell array and
               %%%% transform into nominal?
               if(isprop(obj(1), param))  % first look in Trials; don't need it if not looking only in trials (this and next 3 lines are all that are needed if no search in Params)
                   allVals = [obj.(param)]; 
%                elseif(isfield(obj(1).Params, param)) % then look in params struct  
%                     for k=1:length(obj)
%                         if(isfield(obj(k).Params, param)) % have to re-check because sometimes params may be empty
%                             allVals{k} = obj(k).Params.(param);
%                         end
%                     end
               elseif(isprop(obj(1).Session, param)) % then look in parent Session  
                   %allVals = nominal(nan(1,length(obj)));
                    for k=1:length(obj)
                        allVals(k) = obj(k).Session.(param);
                    end
               else
                   error(['error: field named "',param{:},'" not found' ]);
               end
           end
        end
%%
        function trials = findTrials(obj, expr, range) % get trials based on some criteria; Range is index limits to avoid index out of range [0,-1] *** test
            sessArray = parse1d(obj,'Ses_id');         % parse into sessions b/c any comp that looks across trials (k-1) will otherwise have cross-session artifacts
            trials = [];
            for k=1:length(sessArray)
                trials = [trials, findTrialsHelper(sessArray{k}, expr, range)];
            end
        end
%%
        function trials = findTrialsHelper(obj, expr, range) % helper function that works on session-by-session data; set up this way to simplify calling syntax as 'obj.x == y'
           if nargin<3; range = [0,0]; end;
           indx= zeros(1, length(obj));   % initialize for speed
           for k=1+range(1):length(obj)+range(2)
              if(eval(expr))
                  indx(k) = 1;
              end
           end
           trials = obj(logical(indx));
        end
%%
        function vec = ispropInAnyLevel(obj, pname) % test if it is a property in any level
            try
                [~] = getPropVals(obj, pname);
                vec = ones(1,length(obj));
            catch
                vec = zeros(1,length(obj));
            end
        end
%%
        function trials = findDS(obj, expr, range);  % use dataset array - might not need since findTrials does most of this but a bit slower
        
        end
%%
        function trials = removeTrialsWithEmptyOrNaNField(obj,propName); % removes trials that have no data or are nan for values in property 'propName'
            out = cellfun(@double, {obj.(propName)}, 'Uniformoutput', 0);
            indx = ~( cellfun(@isnan, out) | cellfun(@isempty, out) );
            trials = obj(logical(indx));
            
            %%% rewrote with above code b/c if mix of nominal & double,
            %%% can't concatonate & throws error; also should be faster
%             indx = ones(1,numel(obj));
%             switch class([obj.(propName)])
%                 case 'nominal'
%                     for m=1:numel(indx); 
%                         indx(m) = ~isempty(obj(m).(propName));  % should be a smarter way to do this - can't figure out how to deal output of 'list' to cell array
%                     end;
%                     trials = obj(logical(indx));
%                 otherwise 
%                     % have to do it this way b/c some props may have mutliple numbers (liking or rewards) and concatonation messes up indexing
%                     for m=1:numel(indx);
%                         if(isempty(obj(m).(propName))); indx(m)=0;
%                         elseif(isnan(obj(m).(propName))); indx(m) = 0;
%                         end
%                     end;
%                     trials = obj(logical(indx));
%             end
            
        end
%%
        function check_flg = checkTrialOrderInSession(obj) % checks for ascending trial index order in each session
            sessArray = parse1d(obj,'Ses_id');    
            check_flg = ones(1,length(sessArray));
            for k=1:length(sessArray)
                %check_flg(k) = 1;
                if(~isempty(find(diff([sessArray{k}.Trial_indx])<1,1)))
                    disp(['Trials out of order in session: ',char(sessArray{k}(1).Session.Ses_id)])
                    check_flg(k) = 0; % could make this an error or warning
                end
            end
            
        end
%%
        function setConsecutiveTrialsAndTestDays(obj) % sets the testing day (first date is 1, and so on) and trials as if all sessions were concatonated; also sets cumulitive time  
           % set cumulitive time within each session (seconds)
           sessArray1 = parse1d(obj,'Ses_id');  
           for k=1:numel(sessArray1)
               TrialStartTimes = num2cell(([sessArray1{k}.Trial_start_ts]-sessArray1{k}(1).Trial_start_ts)./(1000));
               [sessArray1{k}.SessionTime] = TrialStartTimes{:};  %
               % compute countown time until session ends (negative seconds)
               TrialStartTimes_fromLast = num2cell(([sessArray1{k}.Trial_start_ts]-sessArray1{k}(end).Trial_start_ts)./(1000));
               [sessArray1{k}.SessionCountdownTimeUntillSessionEnd] = TrialStartTimes_fromLast{:};
           end
           for k=1:numel(sessArray1)  % set trial number
            TrialIndx = num2cell([1:numel(sessArray1{k})]);
            [sessArray1{k}.TrialInSession_indx] = TrialIndx{:};
           end
           % compute cumulative time of session i since the first ever trial of that rat (across sessions)
           sessArray = parse1d(obj,'Rat_id');      % parse into rats
           for k=1:length(sessArray)
               sesPerRat = parse1d(sessArray{k},'Date'); % parse data from each rat into dates
               dateCell = [];
               for m=1:length(sesPerRat)
                    dateCell = [dateCell, {sesPerRat{m}(1).Date}];  % get dates
               end
               [~,Indx] = sort([dateCell{:}]);
               sortSesPerRat = sesPerRat(Indx);
               allTR = []; all_st_delt = [];
               for m=1:length(sortSesPerRat)
                   [~,I] = sort([sortSesPerRat{m}.Trial_indx]);  % make sure trials are in order
                   sortSesPerRat{m} = sortSesPerRat{m}(I);           
                   [sortSesPerRat{m}.TestDay_idx] = deal(m);
                   st = [sortSesPerRat{m}.Trial_start_ts];   % get start times of trials
                   st_delt = [0, diff(st)]./1000;            % time between trials in seconds
                   all_st_delt = [all_st_delt, st_delt];     % concatonate across sessions here
                   allTR = [allTR, sortSesPerRat{m}.getTrials()];
               end
               % compute time elapsed across sessions
               ConsTrialIndx = num2cell([1:numel(allTR)]);   % set index of consecutive trials
               [allTR.ConsecTrialAcrossSess_indx] = ConsTrialIndx{:};
               %%% indx = cumsum([1, cellfun(@numel, sortSesPerRat)]); these are the indicies of the start of each session (should be 0 in elapsed time per session
               %%% indx = indx(1:end-1);
               st_cum = num2cell(cumsum(all_st_delt));   % cumulitive sum of elapsed time
               [allTR.ElapsedTime] = st_cum{:};
           end   
        end
%%      %task specific functions
        function computeNumRewards(obj, selinoidID) % cout the number of rewards given for a particular selinoid on all wells
            for k=1:length(obj)
                obj(k).NumRewards = sum(~isnan([[obj(k).L_reward_ts],[obj(k).R_reward_ts]])); % # have to cahnge if adding more reward wells
            end
        end
%%        
        function computeRespTimeAndLicks(obj)   % compute reaction times, response times, and number of licks
            ReactTime=num2cell([obj.Odor_port_off_ts]-[obj.Tone_on_ts]);
            %ReactTime = cellfun(@(x,y) x-y, obj.getPropVals('Odor_port_off_ts'), obj.getPropVals('Tone_on_ts'),'UniformOutput',0); % cellfun faster than for loop
            [obj.ReactTime] = ReactTime{:};   % alternate to 'deal'
            %RespTime = cellfun(@(x,y) x-y, obj.getPropVals('Well_on_ts'), obj.getPropVals('Odor_port_off_ts'),'UniformOutput',0); % cellfun faster than for loop
            RespTime=num2cell([obj.Well_on_ts]-[obj.Odor_port_off_ts]);
            [obj.ResponseTime] = RespTime{:};   % alternate to 'deal'
            InitTime = num2cell([obj.Odor_port_on_ts]-[obj.Trial_start_ts]);
            [obj.InitTime] = InitTime{:};
            for k=1:length(obj)
                obj(k).NumPreTrialPokes = sum(obj(k).Odor_port_all_on_ts<obj(k).Odor_port_on_ts); % count odor port entires prior to the one after start of trial
                obj(k).TotalPokes = numel(obj(k).Odor_port_all_on_ts);
                obj(k).NumPreTrialRewardWellPokes = sum([obj(k).R_well_all_on_ts < obj(k).Odor_port_on_ts, obj(k).L_well_all_on_ts < obj(k).Odor_port_on_ts]); % cound well entries prior to the good odor poke
            end
    
           % compute total time in all reward wells
           for k=1:length(obj)
                    obj(k).TimeInWell =  max([obj(k).R_well_all_off_ts, obj(k).L_well_all_off_ts]) - obj(k).Well_on_ts; % total time in/around well is last offset - onset
           end
              
           
           % compute inter-trial interval (from reward time to poke)
           for k=2:length(obj)
               if(obj(k).Ses_id == obj(k-1).Ses_id) % only within same session
                    obj(k).InterTrialInterval = obj(k).Odor_port_on_ts - max([obj(k-1).House_light_on_ts, obj(k-1).L_reward_ts, obj(k-1).R_reward_ts]); % this is the reward memory time, from reinforcement to next poke
               else
                   obj(k).InterTrialInterval = nan;
               end
           end
           % Odor_port_on_ts
           
           % count any poke->go behaviours outside of valid trials (uncounted trials)
           [obj.UncountedGO] = deal(0); % assume none, then detect occurances: otherwise, some will be all NaN and plotting/stats functions will not work properly
           for k=1:length(obj)-1
              wellsOn = [obj(k).R_well_all_off_ts, obj(k).L_well_all_off_ts];
              obj(k).PokeBeforeLight = numel(find(wellsOn > min(obj(k).Odor_port_all_on_ts) & wellsOn <  max(obj(k).Odor_port_all_on_ts)));  % if the pattern .... poke -> well -> poke is found, this means rat poked before next trial 
              if(obj(k).PokeBeforeLight>0)   % if there is an extra poke at the end of the trial, look to see if rat goes to the well on this trial, or to the well before a nose poke on the next trial.
                  num_well_firstNextTrial = find(wellsOn > min(obj(k+1).Odor_port_all_on_ts));
                  num_well_afterExtraPokeThisTrial = find(wellsOn >  max(obj(k).Odor_port_all_on_ts));
                  if(~isempty([num_well_afterExtraPokeThisTrial, num_well_firstNextTrial]));
                     obj(k).UncountedGO = 1;
                  end
              end
           end       
            
           % compute lick stats 
            for k=1:length(obj)
                obj(k).LicksTotal = sum(~isnan([obj(k).L_lick_ts, obj(k).R_lick_ts])); % compute total licks
                if obj(k).LicksTotal == 0;
                    obj(k).LicksTotal = nan;    % this is so we don't average supurious '0's' e.g. not averaging 0 when there is no response
                else 
                    obj(k).LicksBeforeRew = sum([obj(k).L_lick_ts, obj(k).R_lick_ts]<min([obj(k).L_reward_ts,obj(k).R_reward_ts])); 
                end
            end
            % get rid out outlier lick counts
            indx = find([obj.LicksBeforeRew] > nanmean([obj.LicksBeforeRew]) + 4*nanstd([obj.LicksBeforeRew]));
            if(~isempty(indx))
                [obj(indx).LicksBeforeRew] = deal(nanmean([obj.LicksBeforeRew]));  % get rid of any outliers
            end
            indx = find([obj.LicksTotal] > nanmean([obj.LicksTotal])+4*nanstd([obj.LicksTotal]));
            if(~isempty(indx))    
                [obj(indx).LicksTotal] = deal(nanmean([obj.LicksTotal]));  % get rid of any outliers       
            end
            
        end
%%        
        function computeSessPerformance(obj) % compute fracton correct ....
            sessArray = parse1d(obj,'Ses_id');
            for k=1:length(sessArray)
                %sessArray{k}.Session.NumTrials = length(sessArray{k});
                sessArray{k}.Session.PercentCorrectChoice = length(sessArray{k}.getTrials('Wrong_choice_flg',0))/length(sessArray{k}); %*** this relys on aquisition code - need to check here?
                
                numBoli = [sessArray{k}.getPropVals('L_reward_ts'), sessArray{k}.getPropVals('R_reward_ts')];
                numBoli = sum(~isnan(numBoli)); % count the number of not nan
                sessArray{k}.Session.MeanRewardsPerTrial = numBoli/length(sessArray{k});   % number of boli expected per trial
            end
        end
        % todo - should change setBlockID to accept switch for property
        % defining blocks: size, odor mapping, delay, ...
%%        
        function setBlockID(obj, minNumberOfTrials)  % set block id based on some criteria; MinNumber (optional) is minimum consecutive number to consider (to filter out transient states)
            if nargin < 2; minNumberOfTrials=1; end
            [obj.Block_id] = deal(nominal('N'));  % initialize to Neither
            %RnumBol =  cell2mat(obj.getPropVals('R_reward_num'));
            %LnumBol =  cell2mat(obj.getPropVals('L_reward_num'));
            RnumBol = [obj.R_reward_num];
            LnumBol = [obj.L_reward_num];
            if (max([RnumBol>LnumBol])>0); [obj(logical([RnumBol>LnumBol])).Block_id] = deal(nominal('Big_R')); end % concatonate other criterial like long if both big & long - otherwise have other sequential criterial that follows these lines
            if (max([RnumBol<LnumBol])>0); [obj(logical([RnumBol<LnumBol])).Block_id] = deal(nominal('Big_L')); end
            sessArray = parse1d(obj,'Ses_id'); % parse into sessions
            for k=1:length(sessArray)          % loop over sessions
                blk = sessArray{k}.getPropVals('Block_id');
                %blk_trns = [1, find(strcmpi(blk(2:end),blk(1:end-1))==0)+1, length(blk)]; % indicies of block transitions; must append index 1 and end (length of Trials)
                blk_trns = [1, find(diff([blk(2:end)==blk(1:end-1)])>0), length(blk)]; %indicies of block transitions: appeend 1 and last index
                blk_num = 0;
                for m=2:length(blk_trns)
                    if blk_trns(m)-blk_trns(m-1) >= minNumberOfTrials % only blocks with at least specified number of trials
                        blk_num = blk_num +1;
                        [sessArray{k}(blk_trns(m-1):blk_trns(m)).Block_number] = deal(blk_num); % assign block number
                        blk_ind = num2cell([1:(blk_trns(m)-blk_trns(m-1))+1]);
                        [sessArray{k}(blk_trns(m-1):blk_trns(m)).Block_trial_indx] = blk_ind{:}; % assign block-relative indicies
                    else
                        [obj(blk_trns(m-1):blk_trns(m)).Block_id] = deal(nominal('N')); % if there were not at least min number of trials, then erase block id.
                    end
                end  
            end
        end
%%        
        function setBlockGoNoGo(obj, minNumberOfTrials)  % set block ID for odor reversal in go nogo; 
                        % todo: merge with function 'setBlockID'
                        %  This could be done with fewer lines by using
                        %  'unique' on concat of odor & reward_id; but this
                        %  misses error checking and min block size
                        %  
                        %
            if nargin < 2; minNumberOfTrials=1; end % minimum number of trials to define a block
            [obj.Block_id] = deal(nominal('N'));  % initialize to Neither
            sessArray = parse1d(obj,'Ses_id'); % parse into sessions
            for k=1:length(sessArray)          % loop over sessions
                % determine well
                numReinfR = sum([sessArray{k}.getPropVals('R_reward_num')]>0); % number reinforcements R
                numReinfL = sum([sessArray{k}.getPropVals('L_reward_num')]>0);
                if(numReinfR > numReinfL)
                    goWell_str = 'R';
                    NoGo_selID = max([sessArray{k}.getPropVals('R_reward_id')]);           % find selinoid for no-go (1=sucrose, so either 2 or 3)
                    NoGoInstruct_vec = [[sessArray{k}.getPropVals('R_reward_id')]]==NoGo_selID;  % find all trials where this selinoid would open if response
                    if(numReinfL > 5); warning(['G0/NoGO task: reinforcements can be obtained in multiple wells; ratID: ',char(sessArray{k}(1).Rat_id),'; Date: ', ...
                            char(sessArray{k}(1).Date),'; # for left: ',num2str(numReinfL),' !!!!']); end % do some error checking on behavioural setup
                else
                    goWell_str = 'L';
                    NoGo_selID = max([[sessArray{k}.getPropVals('L_reward_id')]]);           
                    NoGoInstruct_vec = [[sessArray{k}.getPropVals('L_reward_id')]]==NoGo_selID;  
                    if(numReinfR > 5); warning(['G0/NoGO task: reinforcements can be obtained in multiple wells; ratID: ',char(sessArray{k}(1).Rat_id),'; Date: ', ...
                            char(sessArray{k}(1).Date),'; # for right: ',num2str(numReinfR),' !!!!']); end % do some error checking on behavioural setup
                end
                indx = find(NoGoInstruct_vec>0);
                NoGoOdor = [sessArray{k}(indx).getPropVals('Odor1_id')];
                idx_i=find(diff(NoGoOdor)~=0);  % this is index of indicies; have to convert it back
                blk_trns = [1, indx(idx_i+1), length(NoGoInstruct_vec)]; % vector of indicies defining blocks (first and last index of array are obligatory); +1 is b/c diff operator shifts to left by 1
                blk_num = 0;
                for m=2:length(blk_trns)
                    if blk_trns(m)-blk_trns(m-1) >= minNumberOfTrials % only blocks with at least specified number of trials
                        blk_num = blk_num +1;
                        [sessArray{k}(blk_trns(m-1):blk_trns(m)).Block_number] = deal(blk_num); % assign block number
                        blk_ind = num2cell([1:(blk_trns(m)-blk_trns(m-1))+1]);
                        [sessArray{k}(blk_trns(m-1):blk_trns(m)).Block_trial_indx] = blk_ind{:}; % assign block-relative indicies
                    else
                        [obj(blk_trns(m-1):blk_trns(m)).Block_id] = deal(nominal('N')); % if there were not at least min number of trials, then erase block id.
                    end
                end
                % label trials as errant 'go' or 'nogo'; could do this
                % seperatey from block assignment, but would then have to
                % recompute or store the correct goWell
                nogoTR = sessArray{k}(NoGoInstruct_vec);    % trials instructed nogo (where quinine would be delivered)
                [nogoTR.ErrGo_flg] = deal(0);               % initialize
                if(sum([nogoTR.Well_id]==goWell_str) >1)    % have to prevent attempting to assign to empty set
                    [nogoTR( [nogoTR.Well_id]==goWell_str ).GoNoGoType] = deal('ErrGo');    % if they went to well -> err
                    [nogoTR( [nogoTR.Well_id]==goWell_str ).ErrGo_flg] = deal(1); 
                end
                if(sum([nogoTR.Well_id]~=goWell_str) >1)
                    [nogoTR( [nogoTR.Well_id]~=goWell_str ).GoNoGoType] = deal('OkNoGo');
                end
               
                goTR = sessArray{k}(~NoGoInstruct_vec);    % trials instr
                [goTR.OkGo_flg] = deal(0);                 % initialize
                if(sum([goTR.Well_id]==goWell_str) >1)
                    [goTR( [goTR.Well_id]==goWell_str ).GoNoGoType] = deal('OkGo');
                    [goTR( [goTR.Well_id]==goWell_str ).OkGo_flg] = deal(1);
                end
                if(sum([goTR.Well_id]~=goWell_str) >1)
                    [goTR( [goTR.Well_id]~=goWell_str ).GoNoGoType] = deal('ErrNoGo');
                end
                
            end
        end
        
       % function setNoGoLabels(obj)  % label trials in which rats performed errant 'go' or 'nogo' 
       %     
       %     
       % end
       %function computeNoResponse(obj) % compute prob of not going
       %    [obj.OkGo_flg] = deal(0);
       %    [obj.OkGo_flg] = deal(1);
       %end
%%       
        function setHighValWellBasedOnNumBolus(obj) % sets Property  
            [obj.HighValWell_id] = deal(nominal('N'));  % initialize to Neither
            RnumBol = [obj.R_reward_num];
            LnumBol = [obj.L_reward_num];
            if (max([RnumBol>LnumBol])>0); 
                [obj(logical([RnumBol>LnumBol])).HighValWell_id] = deal(nominal('R')); 
            end% if more boli to R
            if (max([RnumBol<LnumBol])>0); 
                [obj(logical([RnumBol<LnumBol])).HighValWell_id] = deal(nominal('L')); 
            end
              % set prop ChoseHighVal_flg =1 if rat chose high value action
            %out = cellfun(@(x,y) single(strcmpi(x,y)), obj.getPropVals('HighValWell_id'), obj.getPropVals('Well_id'),'UniformOutput',0); % cast to single to allow math later (not defined for logical)
            out = num2cell([obj.HighValWell_id]==[obj.Well_id]); % conver to cell in order to deal into structure
            [obj.ChoseHighVal_flg] = out{:};
        end
%%
        function computeResponseEntropy(obj,len) % len is how many responses to use for a 'word'
            sessArray = parse1d(obj,'Ses_id');
            for k=1:length(sessArray)
                resp_bin = [sessArray{k}.Well_id]; % binary response
                %[sessArray{k}.RespEntropy] = deal(CompEntropyOnBinSeq(resp_bin,len));
                sessArray{k}.Session.RespEntropy = CompEntropyOnBinSeq(resp_bin,len);
            end 
        end
%%        
        function computeWSLS(obj) % computes probability of win stay lose switch 
            sessArray = parse1d(obj,'Ses_id');
            for k=1:length(sessArray)
               % WSLS = [sessArray{k}(2:end).Well_id]==[sessArray{k}(1:end-1).Rewarded_well]; % Win stay lose switch is when rat choice on trial n follows the rewarded well on trian n-1
               % PWSLS = sum(WSLS)/numel(WSLS); 
                %%%[sessArray{k}.ProbWSLS] = deal(PWSLS);
                ch=[sessArray{k}.Well_id];
                del_ch = ch(2:end)~=ch(1:end-1);  % delta (diff) response
                err=[sessArray{k}.Error_flg];
                if(sum(err)<1)
                    warning(['Error in computeWSLS: no losses in session for subject', char(obj(1).Session.Rat_id),' on ', char(obj(1).Session.Date)]);
                    warning(['Mode is: [', num2str(unique([sessArray{k}.Mode])),']; and rewarded well is set to: [', char(unique([sessArray{k}.Rewarded_well])),']']);
                end
                PWSLS = sum(err(1:end-1) == del_ch)/numel(del_ch); % win stay is where err==0 and del_ch ==0; lose switch is where err==1 & del_ch ==1; so count common elements
                
                sessArray{k}.Session.ProbWSLS = PWSLS;
                sessArray{k}.Session.ProbLoseSwitch = sum(err(1:end-1) & del_ch )/sum(err); % probability of lose switch is where an error occured, and the next choice was different
                sessArray{k}.Session.ProbLoseStay = sum(err(1:end-1) & ~del_ch )/sum(err); % probability of lose stay is where an error occured, and the next choice was stay
                sessArray{k}.Session.ProbWinStay = sum(~err(1:end-1) & ~del_ch )/sum(~err); % probability of win stay is where no error occured (err=0), and choice did not switch (del_ch=0)
                
                
              %%%%% new computations for lose-stay on trial by trial basis
                tr = sessArray{k};
                                        
               % throw logical flag if rat sampled both wells prior to choice
               for k=2:length(tr)
                  R_entries_pre = [tr(k).R_well_all_off_ts]<tr(k).Odor_port_on_ts;
                  R_entries_postLastT = [tr(k-1).R_well_all_off_ts]>tr(k-1).Odor_port_on_ts;
                  L_entries_pre = [tr(k).L_well_all_off_ts]<tr(k).Odor_port_on_ts;
                  L_entries_postLastT = [tr(k-1).L_well_all_off_ts]>tr(k-1).Odor_port_on_ts;
                  if( nansum([[R_entries_pre],[R_entries_postLastT]])>0 && nansum([[L_entries_pre],[L_entries_postLastT]]) > 0); % if rats checked both R and L wells between last unpoke (last trial) and poke on current trial, then sampled both wells
                      tr(k).SampledBothWellsBeforeTrial_flg = 1;
                  else
                      tr(k).SampledBothWellsBeforeTrial_flg = 0;
                  end
               end

               % compute the last well sampled before the current trial
               for k=2:length(tr)
                  indx = max(find([tr(k).R_well_all_off_ts]<tr(k).Odor_port_on_ts));   if(isempty(indx)); R_entries_pre = []; else R_entries_pre = tr(k).R_well_all_off_ts(indx); end;
                  indx = max(find([tr(k-1).R_well_all_off_ts]>tr(k-1).Odor_port_on_ts)); if(isempty(indx)); R_entries_postLastT = []; else R_entries_postLastT = tr(k-1).R_well_all_off_ts(indx); end;
                  indx = max(find([tr(k).L_well_all_off_ts]<tr(k).Odor_port_on_ts));  if(isempty(indx)); L_entries_pre = []; else L_entries_pre = tr(k).L_well_all_off_ts(indx); end;
                  indx = max(find([tr(k-1).L_well_all_off_ts]>tr(k-1).Odor_port_on_ts)); if(isempty(indx)); L_entries_postLastT = []; else L_entries_postLastT = tr(k-1).L_well_all_off_ts(indx); end;

                  if(max([R_entries_pre,R_entries_postLastT,0]) > max([L_entries_pre,L_entries_postLastT,0])); % have to add 0 so as to avoid [] from max function if any elements are empty
                      tr(k).LastWellCheckedBeforeTrial_id = 'R';
                  elseif(max([R_entries_pre,R_entries_postLastT,0]) < max([L_entries_pre,L_entries_postLastT,0]))
                      tr(k).LastWellCheckedBeforeTrial_id = 'L';
                  end

               end     

               %IsSequentialFromGoodTrial_flg  
                for m=2:numel(tr)
                    % compute lose-switch for all trials
                    if( tr(m-1).Error_flg && tr(m).Well_id ~= tr(m-1).Well_id ) 
                        tr(m).LoseSwitch_all_flg = 1;
                    elseif( tr(m-1).Error_flg && tr(m).Well_id == tr(m-1).Well_id ) 
                        tr(m).LoseSwitch_all_flg = 0;
                    end
                    
                    % compute lose switch only for trials without sampling of both wells
                    if(tr(m).SampledBothWellsBeforeTrial_flg == 0) % only compute for trials in which the rat did not sample both wells
                        %if( tr(m-1).Error_flg && strcmp(char(tr(m).Well_id), char(tr(m-1).Well_id) )) 
                        if( tr(m-1).Error_flg && tr(m).Well_id ~= tr(m-1).Well_id ) 
                            tr(m).LoseSwitch_NoMultWellResp_flg = 1;
                        elseif( tr(m-1).Error_flg && tr(m).Well_id == tr(m-1).Well_id ) 
                            tr(m).LoseSwitch_NoMultWellResp_flg = 0;
                        end
                    end

                    % compute lose-switch based on last sampled well - even if out of the official trial
                    if(~isnan(tr(m).LastWellCheckedBeforeTrial_id))
                        %if( tr(m-1).Error_flg && ~strcmp(char(tr(m).Well_id), char(tr(m).LastWellCheckedBeforeTrial_id) )) % cast to strings to catch exception when one is NaN 
                        if( tr(m-1).Error_flg && tr(m).Well_id ~= tr(m).LastWellCheckedBeforeTrial_id )
                            tr(m).LoseSwitch_BasedOnOutOfTrialResp_flg = 1;
                        %elseif( tr(m-1).Error_flg && strcmp(char(tr(m).Well_id), char(tr(m).LastWellCheckedBeforeTrial_id) ) ) 
                        elseif( tr(m-1).Error_flg && tr(m).Well_id == tr(m).LastWellCheckedBeforeTrial_id )
                            tr(m).LoseSwitch_BasedOnOutOfTrialResp_flg = 0;
                        end
                    end
                    
                    % compute win stay for all trials
                    if( ~tr(m-1).Error_flg && tr(m).Well_id ~= tr(m-1).Well_id ) 
                        tr(m).WinStay_all_flg = 0;
                    elseif( ~tr(m-1).Error_flg && tr(m).Well_id == tr(m-1).Well_id ) 
                        tr(m).WinStay_all_flg = 1;
                    end
                    
                    % compute win-stay based on last sampled well - even if out of the official trial
                    if(~isnan(tr(m).LastWellCheckedBeforeTrial_id))
                        %if( tr(m-1).Error_flg && ~strcmp(char(tr(m).Well_id), char(tr(m).LastWellCheckedBeforeTrial_id) )) % cast to strings to catch exception when one is NaN 
                        if( ~tr(m-1).Error_flg && tr(m).Well_id == tr(m).LastWellCheckedBeforeTrial_id )
                            tr(m).WinStay_BasedOnOutOfTrialResp_flg = 1;
                        %elseif( tr(m-1).Error_flg && strcmp(char(tr(m).Well_id), char(tr(m).LastWellCheckedBeforeTrial_id) ) ) 
                        elseif( ~tr(m-1).Error_flg && tr(m).Well_id ~= tr(m).LastWellCheckedBeforeTrial_id )
                            tr(m).WinStay_BasedOnOutOfTrialResp_flg = 0;
                        end
                    end
                    
                    % compute win stay only for trials without sampling of both wells
                    if(tr(m).SampledBothWellsBeforeTrial_flg == 0) % only compute for trials in which the rat did not sample both wells
                        if( ~tr(m-1).Error_flg && tr(m).Well_id == tr(m-1).Well_id ) 
                            tr(m).WinStay_NoMultWellResp_flg = 1;
                        elseif( ~tr(m-1).Error_flg && tr(m).Well_id ~= tr(m-1).Well_id ) 
                            tr(m).WinStay_NoMultWellResp_flg = 0;
                        end
                    end

                    if(tr(m-1).Error_flg)
                        tr(m).LastTrialRewarded_flg == 0;
                    else
                        tr(m).LastTrialRewarded_flg == 1;
                    end
                    
                    if(tr(m).Well_id == tr(m-1).Well_id)
                        tr(m).ChoiceSwichFromLastTrial == 0;
                    else
                        tr(m).ChoiceSwichFromLastTrial == 1;
                    end
                    
                end
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                
            end 
        end
        
%%        
        function computeProbR(obj) % compute prob of response to the right
            sessArray = parse1d(obj,'Ses_id');
            for k=1:length(sessArray)
                R = [sessArray{k}(2:end).Well_id]=='R'; % Win stay lose switch is when rat choice on trial n follows the rewarded well on trian n-1
                PR = sum(R)/numel(R); 
                sessArray{k}.Session.ProbR = PR;
            end 
        end
%%        
        function computePropRepeat(obj) % compute probability of repeating same response in consecutive trials
            sessArray = parse1d(obj,'Ses_id');
            for k=1:length(sessArray)
                same = [sessArray{k}(2:end).Well_id]==[sessArray{k}(1:end-1).Well_id]; % true when rat makes same response in consecutive trials
                Psame = sum(same)/numel(same); 
                %[sessArray{k}.ProbSame] = deal(Psame);
                sessArray{k}.Session.ProbSame = Psame;
            end 
        end
%%        
        function computeLogitRegression_symmetricRew(obj, n) % n is the number of past trials to use in regression
            sessArray = parse1d(obj,'Ses_id');
%             for k=1:length(sessArray)
%                % make Y vector of choices
%                 %Y = [[sessArray{k}(n+1:end).Well_id]=='R']'; % make an n x c matric (number observations x number of responses: here is 2 for 'R' or 'L')
%                % make Y vector of switches
%                   Y = [[sessArray{k}(n+0:end).Well_id]=='R']'; % index from n+0 rather than n+1 (as before) because diff operation eliminates the first element
%                   Y = double(~diff(Y))+1;  % vector where switches are coded as '1' and repeated responses are '0'; has to be non-zero double for regression algorithm
%                % Y(:,2)=~Y;
%                %%% get total # of rewards - only for non-symmetric rewards
%                % numRew = zeros(1,length(sessArray{k}));
%                % for m=1:length(sessArray{k})
%                %     numRew(m) = length([sessArray{k}(m).R_reward_ts, sessArray{k}(m).L_reward_ts])-1;
%                % end
%                % make X vector of predictors 
%                for m=1:n
%                     X(2*m-1,:) = ([sessArray{k}(n-m+1:end-m).Well_id]=='R').*2-1;  % this codes response in terms of +1/-1
%                     X(2*m,:) = ([sessArray{k}(n-m+1:end-m).Rewarded_well]=='R').*2-1; % codes rewards as +1/-1;
%                     %X(2*m,:) = (([sessArray{k}(n-m+1:end-m).Rewarded_well]=='R').*2-1).*numRew(n-m+1:end-m); % codes reward by size (magnitude) and location +/-
%                     %X(2*m-1,:) = ([sessArray{k}(n-m+1:end-m).Well_id]=='R');  % this codes response in terms of +1/0
%                     %X(2*m,:) = ([sessArray{k}(n-m+1:end-m).Rewarded_well]=='R');
%  
%                end
%                 X = X';
%                 try
%                     [B, var, stats] = mnrfit(X,Y);
%                     fitY = mnrval(B,X);         % compute predicted (fit) output from coefecients and x data
%                     fitY = round(fitY(:,1));    % just look at first column (choice to 'R')
%                     accuracy = sum(fitY==Y(:,1))/numel(fitY); % fraction of correct predictions by model
% 
%                     sessArray{k}.Session.LogitRegStats = stats; %LogitRegFitAcc
%                     sessArray{k}.Session.LogitRegB = B; % B=[const, Rc(t-1), Cc(t-1), Rc(t-2),....
%                     sessArray{k}.Session.LogitRegFitAcc = accuracy;
%                 
%                     %[sessArray{k}.RespEntropy] = deal(CompEntropyOnBinSeq(resp_bin,len));                    
%                 catch
%                     warning(['can not process regression for ',char(sessArray{k}(1).Rat_id),'!']);
%                 end
%                 clear X Y fitY; 
%                 
%            end   
               for k=1:length(sessArray)
                   resp = [sessArray{k}.Well_id]=='R';
                   rewWell = [sessArray{k}.Rewarded_well]=='R';
                   rew = resp==rewWell;
                   
                 % for testing
                    %resp = 1:100;
                    %rew = 201:300;



               % make Y vector of choices
                %Y = [[sessArray{k}(n+1:end).Well_id]=='R']'; % make an n x c matric (number observations x number of responses: here is 2 for 'R' or 'L')
               % make Y vector of switches
                %  Y = resp(n:end)'; % index from n+0 rather than n+1 (as before) because diff operation eliminates the first element; have to start at n to have the previous n for prediction
                %  Yold = double(~diff(Y))+1;  % vector where switches are coded as '1' and repeated responses are '0'; has to be non-zero double for regression algorithm
                  
                  Ypred = [1, double(~diff(resp))*2-1]; % assign first response as not a switch
                  Y = (Ypred(n+1:end)'+1)/2+1; % start predicting on n+1 trial b/c need previous n for predictors
                  
                  % maxe a matrix of predictors - alternating terms for choice and reward: 1st value of Y is nth+1 trial - first row of X is [n, n-1, n-1, ... 1]
                   for m=1:n
                       % for testing
                        %%X(2*m-1,:) = resp(n-m+1:end-m);  
                        %%X(2*m,:) = rew(n-m+1:end-m);

                        %X(2*m-1,:) = resp(n-m+1:end-m).*2-1;  % this codes response in terms of +1/-1 for R/L
                        
                        X(2*m-1,:) = Ypred(n-m+1:end-m);  % this codes reponses as +1/-1 for switches from last trial
                        X(2*m,:) = rew(n-m+1:end-m).*2-1; % codes rewards as +1/-1;

                        %X(2*m,:) = (([sessArray{k}(n-m+1:end-m).Rewarded_well]=='R').*2-1).*numRew(n-m+1:end-m); % codes reward by size (magnitude) and location +/-
                        %X(2*m-1,:) = ([sessArray{k}(n-m+1:end-m).Well_id]=='R');  % this codes response in terms of +1/0
                        %X(2*m,:) = ([sessArray{k}(n-m+1:end-m).Rewarded_well]=='R');

                   end
                   % double check - another way to make the X matrix
%                    for mm=1:numel(Ypred)-n
%                        for jj = 1:n;
%                             XX(mm,jj*2-1) = Ypred(mm+n-jj);
%                             XX(mm,jj*2) = rew(mm+n-jj)*2-1;
%                        end
%                    end
                   
                    X = X';
                    try
                        [B, var, stats] = mnrfit(X,Y);
                        fitY = mnrval(B,X);         % compute predicted (fit) output from coefecients and x data
                        fitY = round(fitY(:,1))+1;    % just look at first column (choice to 'R')
                        accuracy = sum(fitY==Y(:,1))/numel(fitY) % fraction of correct predictions by model

                        sessArray{k}.Session.LogitRegStats = stats; %LogitRegFitAcc
                        sessArray{k}.Session.LogitRegB = B; % B=[const, Rc(t-n), Cc(t-n), Rc(t-n+1),....
                        sessArray{k}.Session.LogitRegFitAcc = accuracy;

                        %[sessArray{k}.RespEntropy] = deal(CompEntropyOnBinSeq(resp_bin,len));                    
                    catch
                        warning(['can not process regression for ',char(sessArray{k}(1).Rat_id),'!']);
                    end
                    clear X Y fitY; 
                
                end 
        end
%%        
        function setHighRiskByNumberBoli(obj) % when reward prob is not saved, have to look at # boli and occurance
            
            
        end
%%        
        function [NewPropName, out] = compute(obj, expr, level, range, NewPropName) % apply some expression to each
                                          % or set new prop Datax = output
                                          % DataxExpr = expr passed in
                                          % allow user to pass in 'session'
                                          % such that new data is in
                                          % session rather than each trial
                                          % (e.g. % correct)
                                          % level is 'session', or 'trial'
        end
        function RenameProp(obj, property, pattern, replacementName)
            C = cellstr(char([obj.(property)]')); % turn to cell array
            %indx = strcmpi(pattern,C);
            match_vec = strfind(C,pattern);
            indx = ~cellfun(@isempty, match_vec); % find the non-empty ones
            
            [obj(indx).(property)] = deal(nominal(replacementName));
            
        end
        
    end
    
   %%  
     methods(Static)
      function TRarray = parseTRCells(TRarray, ParamCell)
        if(~iscell(TRarray))
            %error('function parseTRCells requires TR input to be cell array');
            TRarray = {TRarray};
        end
        if(~iscell(ParamCell))
            error('function parseTRCells requires params input to be cell array');
        end

        for i=1:numel(ParamCell)
            for j = 1:numel(TRarray)
                TRarray1{j,:} = parse1d(TRarray{j},ParamCell{i});
                %TRarray(j,:) = parse1d(TRarray{j},ParamCell{i});
            end
           
            TRarray = [(TRarray1{:})];
           % TRarray = deal(TRarray1{:});
           % TRarray = reshape(TRarray{:}, numel(TRarray{:}), 1);
        end
        %TRarray = [(TRarray{:})];
        
      end
    end
    
end