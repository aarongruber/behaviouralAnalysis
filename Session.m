classdef Session < dynamicprops% handle
    properties (SetAccess = private)
        Animal = Animal.empty;
        Trials = Trial.empty;
        Directory
        
        % Session Data: must have cap first letter
        Rat_id
        Drug
        Dose
        Notes
        Program_version
        Date
        Time
        Engine
        Ses_id
        Group_id
        
        %Response % by commenting, it will not load this
    end
    properties
        NumTrials = nan;
        PercentCorrectChoice = nan;  % percentage of trials rewarded
        MeanRewardsPerTrial = nan;    % total rewards (integrates multiple boli of reward)
        RespEntropy = nan;
        ProbNoResp = nan;  % probability of poking but not responding in a any reward well
        ProbR = nan;     % prop of responding to 'R'
        ProbWSLS = nan;  % probability of win stay lose switch strategy
        ProbWinStay = nan; % probability of win stay
        ProbLoseSwitch = nan; % probability of lose switch
        %ProbLoseSwitch_InclOutOfTrialResp = nan; % compute prob lose-switch using last well visited - e.g. if more than one sampled
        %ProbLoseSwitch_NoMultWellResp = nan;  % prob of lose-switch including only trials with no multiple well checks prior to trial
        ProbLoseStay = nan; % probability of lose stay
        ProbSame = nan;  % probability of repeating choice
        LogitRegB = nan; % logistical regression coefficients
        LogitRegStats = nan; % logistical regression stats
        LogitRegFitAcc = nan; % logistical regression accuracy 
        TimeToCriterion = nan;
        TempGrouping = nan;
    end
    methods
       function obj = Session(sessiondir,animalsregistry)
           if nargin ~= 0 
               % Constructor, call with ether directory or a specific file.
               % If a directory, the *.BEH.mat file will be loaded. Multiple files is an error. 
               % The second argument is the animal registry
               % that links sessions with the subject indicated in rat_id
               % in the saved data structure.
               %
               %%% Get the name of the animal associated with the sessiondir
               %dircontents = dir(sessiondir);
               %filenames = {dircontents.name};
               %fileidx = cellfun(@(x) strncmp('20',x,2),filenames); % get indicies of fileneames that start with '20' as quick & dirty way to find BEH sessions
               %%% load variable SESSION
               %load(fullfile(sessiondir,filenames{fileidx})); % loads struct named 'session'
               %%% todo: throw message if no file is found
               if(isdir(sessiondir))        % if it is a directory, get all BEH files
                   dircontents = dir([sessiondir,'\*.BEH.mat']);
                   filename = {dircontents.name};
                   if(length(filename)>1); error(['Multiple BEH.mat files found in dir: ',sessiondir,'; use only one .BEH.mat file per directory or supply filename']); end
                   if(length(filename)<1); warning(['No BEH.mat files found in dir: ',sessiondir]); end
               else
                   nIdx = strfind(sessiondir,'\');   % could do in one step with regexp
                   filename = {sessiondir(nIdx(end)+1:end)};
                   sessiondir = sessiondir(1:nIdx(end));
               end

               load(fullfile(sessiondir,char(filename)));
               id = session.info.rat_id;
               animal = Animal.getAnimalFromID(animalsregistry,id); %#ok<PROP>
               if isempty(animal)
                   % We're creating a new animal, but not assigning a group.
                   % Call the constructor for Animal
                   animal = Animal(animalsregistry); %#ok<PROP>
                   animal.ID = session.info.rat_id;
               end
               obj.Animal = animal;
               obj.Animal.addSession(obj);
               obj.Ses_id = nominal([session.info.rat_id,':',session.info.date,':',session.info.time]);
               obj.addSessionData(session);
               obj.addTrials(session);
               obj.NumTrials = length(obj.Trials);
               obj.Directory = sessiondir;
               
           end
       end
       function addSessionData(obj,session)
           %for fl = fieldnames(session.info)'  % this uses all field names in the data
           for fl =  [intersect(fieldnames(session.info)', fieldnames(obj)'), intersect(fieldnames(session.info)', lower(fieldnames(obj))')]% this copies only file names in struct and in the file
               objFName = char(fl);             % cap
               objFName = [upper(objFName(1)), objFName(2:end)];    
               obj.(objFName) = nominal(session.info.(char(fl)));
           end
       end
       function addTrials(obj,session)
           newTrialsArray = Trial(obj,session.tr);  % for speed, pre-allocate memory
           obj.Trials = [obj.Trials newTrialsArray];
       end

       function trials = getTrials(obj,param,val) % try again using find in Trial methods
           trials = [];
           for k=1:length(obj)
               if nargin ==1
                   trials = [trials, obj(k).Trials.getTrials()];
               else
                   trials = [trials, obj(k).Trials.getTrials(param,val)];
               end
           end
       end
       function checkForDuplicateData(obj, param) % check if duplicate datasets are found (e.g. dataset loaded multiple times)
            IDs = [obj.Ses_id];
            IDs_unique = unique(IDs);
            if(numel(IDs)>numel(IDs_unique))
                indx = []; msg = [];
                for k=1:numel(IDs)-1
                    if(find(IDs(k) == IDs(k+1:end) >0));
                        indx = [indx, k];
                        msg = [msg, 'duplicate data sets found - session ID= ',char(IDs(k)), ';  ses indx= ', num2str(k'), '\n']
                    end
                end
                %string = strcat('duplicate data sets found - session ID= ',char(IDs(indx)), ';  ses indx= ', num2str(indx'));
                error(sprintf(msg))
            end   
       end
       
       function setDrugLabelsByDate(obj, DateID, label) % method for setting the 'drug' label on particular dates.
           for k=1:numel(obj)
               s_string = char(obj(k).Ses_id);   % have to pull the date out of the session ID because the date is screwed up in ses object for some reason.
               if(s_string(end-9:end) == DateID)
                   if(~isempty(obj(k).Drug))
                       warning(['Overwriting exsisting field "Drug" in session ID= ',char(obj(k).Ses_id), ';  ses indx= ', num2str(k'), '\n'])
                   end
                   obj(k).Drug = nominal(label); 
               end
           end
       end
       
       function setGroupLabels(obj, RatIdCell, label) % give a cell array of rat Ids and label to set Group_id field for each
           if nargin ==1
               [obj(:).Group_id] = deal([]); % if no inputs, clear all labels
               trs = [obj.Trials];
               trs.setGroupID([]);
           else
               for fl=RatIdCell
                   indx = [obj.Rat_id] == fl{:};
                   if(sum(indx)<1);  
                       error(['In setGroup, name ', fl{:},' is not found in the dataset'])
                   else
                      %indx = cellfun(@(x) ~isempty(x), strfind({char(obj.Rat_id)}, char(fl)));   
                      [obj(indx).Group_id] = deal(nominal(label));
                      trs = [obj(indx).Trials];
                      trs.setGroupID(label);
                   end
               end
           end
       end
       function checkGroupLabels(obj)
           for k=1:numel(obj)
               if(isempty(obj(k).Group_id))
                   disp(['Grouping name not set for session(',num2str(k),') = ',char(obj(k).Ses_id),'; Could be typo in Rat_id = ',char(obj(k).Rat_id)])
               end
           end
       end
       function ses = removeSessionsWithNoGroupID(obj)
           indx = zeros(1,numel(obj));
           for k=1:numel(obj)
               if(~isempty(obj(k).Group_id))
                    indx(k)=1;
               end
           end
           ses = [obj(logical(indx))];
       end
       function ses = removeTrialsWithEmptyField(obj,fieldn) % new session that has only responses (e.g. not: early go or no response)
           for k=1:length(obj)
              obj(k).Trials = obj(k).getTrials()~=obj(k).getTrials(fieldn,'nan');  % get all trials where Well_id is not empty
           end
       end
       function setProp(obj, propName, vals) % vals can be be vector or cell of numbers or strings
           if ~iscell(vals)
               vals = num2cell(vals);
           end
           [obj.(propName)] = vals{:};   % this syntax works in place of deal for cell arrays
       end
       function vals = compute(obj, expr)  % apply expression in session %% could also be function handle %% **** really should use trial method -> remake group >> allows filtering of trials first
            for k=1:length(obj)
               vals(k)= feval(expr); 
            end
       end
       
       %% some logical functions: find, &, |, ~=
       function getSessions(obj,param,val)% return sessions meeting criteria
           % returns sessions in obj (Trials obj) for which param (in Trials, Params, or Session) is val
           % (val can be number, char or empty '')
           if nargin==1 || isempty(param)
               ses = [];
           else
               if(~isempty(val))  % if not empty, find matching values, otherwise return the empty ones (b/c strcmpi does not match '')
                   ses = obj([obj.(param)]==val);
                   % allVals = cellfun(@num2str,allVals,'UniformOutput',0); % convert to char
                   % trials = [obj([strcmpi(allVals, num2str(val))])]; % string comp is more robust since it will allow multiple digit vals or strings
               else
                   if(isnumeric([obj.(param)]) || islogical([obj.(param)]))
                       ses = obj([obj.(param)]==nan);
                   else
                       ses = obj([obj.(param)]=='nan');
                   end
               end
           end
       end
       function ses = and(obj1, obj2) % logical 'and' (intersection of sets): return trils that are in both objects
           if(isempty(obj1) || isempty(obj2))
               ses = [];
           else
               [~, obj1_indx, obj2_indx] = intersect([obj1.Unique_trial_id], [obj2.Unique_trial_id]);
               %[~, obj1_indx, obj2_indx] = intersect(cellfun(@(x) char(x), {obj1.Unique_trial_id}, 'UniformOutput',0),  ...
               %                                      cellfun(@(x) char(x), {obj2.Unique_trial_id}, 'UniformOutput',0));
               ses = obj1(obj1_indx);
           end
       end
       function ses = or(obj1, obj2) % logical 'or' (union of sets)
           [~, obj1_indx, obj2_indx] = union([obj1.Ses_id], [obj2.Ses_id]);
           ses = [obj1(obj1_indx),obj2(obj2_indx)];
       end
       function ses = ne(obj1, obj2) % logical ~= 'not equal to' (exclude from obj1 any that are in obj2)
           if(isempty(obj1))
               ses = [];
           elseif(isempty(obj2))
               ses = obj1;
           else
               [~, obj1_indx] = setdiff([obj1.Ses_id], [obj2.Ses_id]);
               ses = obj1(obj1_indx);
           end
       end
       
       function ses = removeSessonsWithEmptyOrNaNField(obj,propName)
            out = cellfun(@double, {obj.(propName)}, 'Uniformoutput', 0);
            indx = ~( cellfun(@isempty, out) );
            for i=1:numel(indx);             % do it the dumb way to avoid indexing problems with empty cells
                if(indx(i) == 1)
                    if(isnan(out{i}))
                        indx(i) = 0;
                    end
                end
            end
            ses = obj(logical(indx));
       end
       
       function ses = unique(obj) % return unique session from an array of sessions (omitting ones with blank session ids)
           indx = ~( cellfun(@isempty, {obj.Ses_id}) ); % check if any are empty (as this will screw up indexing); this should never happen
           if(sum(indx)<numel(indx))    
               disp('in session.unique: eliminating sesions with no ses_id')
               obj = obj(indx);
           end
           unique_ses_id = getlevels([obj.Ses_id]);                 % find unique tags for session id
           [~,~,indx_ses] = intersect(unique_ses_id, [obj.Ses_id]); % find the indicies of examples
           ses = obj(indx_ses);
       end
       
       
       %% odor-box and behavior specific functions
       function computeRespTimeAndLicks(obj)
            for k=1:length(obj)
               obj(k).Trials.computeRespTimeAndLicks;   % kick it to the Trials method
            end
       end
       function computeSessPerformance(obj)
            for k=1:length(obj)
               obj(k).Trials.computeSessPerformance;   % kick it to the Trials method
            end
       end
       function setHighValWellBasedOnNumBolus(obj) % determine well giving more boli of reward
           for k=1:length(obj)
               obj(k).Trials.setHighValWellBasedOnNumBolus;
           end
       end
       function setBlockID(obj, minNumberOfTrials)  % set BlockID property and index of trials in block; minNumberOfTrials (optional) is min number to qualify as a block
           if nargin < 2;minNumberOfTrials=1; end
           for k=1:length(obj)
               obj(k).Trials.setBlockID(minNumberOfTrials);
           end
       end
       function computeNoRespProb(obj)
           for k=1:length(obj)
                obj(k).ProbNoResp = 1 - sum(~isnan([obj(k).Trials.Well_on_ts]))/numel([obj(k).Trials]);
           end
       end
%        function trials = getTrialsMultCrit(obj,cellIn)
%            %trials = [];
%            trials = obj(1).Trials([]);  % this sets trials to be Trials object in case no valid trials are found (only needed for Trials.params search)
%            % search in Trials and in Trials.params struct
%            for k = 1:length(obj)
%                if(isprop(obj(k).Trials(1), param))  % first look in Trials; don't need it if not looking only in trials (this and next 3 lines are all that are needed if no search in Params)
%                     allVals = {obj(k).Trials.(param)};          % make cell array for multiple lenght data
%                     allVals = cellfun(@num2str,allVals,'UniformOutput',0); % convert to char
%                     indx = [strcmpi(allVals, num2str(val))];
%                elseif(isfield(obj(k).Trials(1).Params, param)) % then look in params struct
%                     for j=1:length(obj(k).Trials)
%                         if(isfield(obj(k).Trials(j).Params, param)) % have to re-check because sometimes params may be empty
%                             if(strcmpi(num2str(obj(k).Trials(j).Params.(param)), num2str(val))) % do it the dumb but robust way
%                                 trials = [trials obj(k).Trials(j)];
%                             end
%                         end
%                     end
%                else
%                    
%                    error(['error: field ',param,' not found' ]);
%                end
%                trials = [trials obj(k).Trials([strcmpi(allVals, num2str(val))])]; % string comp is more robust since it will allow multiple digit vals or strings
%            end
%        end
    end
    
end