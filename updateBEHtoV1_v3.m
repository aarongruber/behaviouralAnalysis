function updateBEHtoV1(~)
%
% function updateBEHtoV1(fname)
%
% updateBEHtoV1() converts all files in the current
%
%  Converts older BEH matlab files containing behavioral data from the Gruber 
%  Lab operant boxes to a standard version v1.
%
%  Note: all field names have first letter capitalized
%
%   Aaron Gruber 2011_7_5

% To do to saved structs from beh code:
%   - Dose for every case
%   - have to add probabilities for different trials & rewards (has to
%   account for compound odors

%% define properties of data structure

    if nargin == 0
        files = dir('*.mat');
    elseif nargin == 1
        files = varargin{1};
    end

    for fn = {files.name}
        Session = convertBEH(fn);
        save([char(fn),'.BEHv1'], '-struct', 'Session');
        clear Session;
    end
end

function Session = convertBEH(fn)
    c=load(char(fn));% load
    s = c.session;   % rename for clarity
% Session level
    Session.Rat_id = []; % ID of rat
    Session.Drug = [];    % Name of drug given (incl. saline)
    Session.Dose = [];   % Dose of drug in mg/kg
    Session.Notes = [];  % Notes about behavior
    Session.Program_version = [];    
    Session.Date = [];   
    Session.Time = [];
    Session.Engine = []; % The name of the .m file for the control logic    
    Session.Response= [];% All output from the Arduino controller as text - for debugging
    Session.Epochs = []; % cell array of named epocs (sleep, task, ..) and times [3:20]
    Session.Electrode_placement = []; % placement in {TT*, structName, AP, ML, DV} for each tetrode with spikes or lfp
    Session.BEHtype = []; % More specific name of behavior (sleep, reversal, distract, switch, discrim)
                              %  e.g. {'sleep1',[3433, 4455]; 'discrim', [4466, 4750]; 'sleep2', [4780, 4800]}
    Session.Data = [];   % more specific information about session hardware & parameter configuration
    
% Trial level
    Session.Trial.Trial_start_ts = [];
    Session.Trial.L_light_on_ts = [];
    Session.Trial.L_light_off_ts = [];
    Session.Trial.R_light_on_ts = [];
    Session.Trial.R_light_off_ts = [];
    Session.Trial.Odor1_id = [];           % could make into a vector or cell array
    Session.Trial.Odor2_id = [];
    Session.Trial.Odor_port_on_ts = [];    % poke_on
    Session.Trial.Odor_port_off_ts = [];   % poke off
    Session.Trial.Odor1_on_ts = [];
    Session.Trial.Odor1_off_ts = [];
    Session.Trial.Odor2_on_ts = [];
    Session.Trial.Odor2_off_ts = [];
    Session.Trial.Tone_on_ts = [];
    Session.Trial.Tone_off_ts = [];
    Session.Trial.Well_on_ts = [];
    Session.Trial.Well_off_ts = [];
    Session.Trial.Well_id = [];
    Session.Trial.L_reward_ts = [];
    Session.Trial.R_reward_ts = [];
    Session.Trial.L_lick_ts = [];
    Session.Trial.R_lick_ts = [];
    Session.Trial.Error_id = [];
    Session.Trial.Error_string = [];
    Session.Trial.Error_flg = [];
    Session.Trial.Odor_port_timeout_flg = [];
    Session.Trial.Odor_port_early_go_flg = [];
    Session.Trial.Wrong_choice_flg = [];
    Session.Trial.Well_early_go_flg = [];
    Session.Trial.Stimulation_on_ts = [];
    Session.Trial.Stimulation_off_ts = [];
    Session.Trial.House_light_on_ts = [];
    Session.Trial.House_light_off_ts = [];

    Session.Trial.Params = struct();

   % moved from params
    Session.Trial.Tone_period = [];
    Session.Trial.Rewarded_well = [];
    Session.Trial.L_reward_id = []; 
    Session.Trial.L_reward_delay = []; % could get this from ts
    Session.Trial.L_reward_num = [];   % could get this from length of well ts    
    Session.Trial.R_reward_id = [];
    Session.Trial.R_reward_delay = [];
    Session.Trial.R_reward_num = [];
    Session.Trial.Error_duration = [];
    Session.Trial.No_errors_flg=[];
    Session.Trial.Mode = [];

   % need to add to data in behaivoral scripts for future compatibility 
    Session.Trial.Trial_type_id = [];  % indicate type of trial; probably best post-hoc; BigR, BigL, SmL, LoR, BigR_SmL, SmL_Dist,
    Session.Trial.Odor1_prob = [];     % probability of odor presentation
    Session.Trial.Odor2_prob = [];     
    Session.Trial.L_reward_prob = [];
    Session.Trial.R_reward_prob = [];

   % params
   
   Session.Trial.Params.Odor1_delay = [];
   Session.Trial.Params.Odor1_duration = [];
   Session.Trial.Params.Odor2_delay = [];
   Session.Trial.Params.Odor2_duration = [];
   Session.Trial.Params.Tone_delay = [];
   Session.Trial.Params.Tone_duration = [];
   Session.Trial.Params.Rewarded_well = [];
   Session.Trial.Params.L_reward_duration = [];
   Session.Trial.Params.L_bolus_delay = [];
   Session.Trial.Params.R_reward_duration = [];
   Session.Trial.Params.R_bolus_delay = [];
   Session.Trial.Params.Trial_delay = [];
   Session.Trial.Params.Error_duration = [];
   Session.Trial.Params.Odor_port_timeout = [];
   Session.Trial.Params.Well_timeout = [];
   Session.Trial.Params.Debounce_up = [];
   Session.Trial.Params.Debounce_down_odor_port = [];
   Session.Trial.Params.Debounce_down_well = [];
   Session.Trial.Params.L_prob = [];
   Session.Trial.Params.R_prob = [];
   Session.Trial.Params.E_prob = [];
   Session.Trial.Params.D_prob = [];
   Session.Trial.Params.Error_correction_flg = [];
   Session.Trial.Params.Max_repetitions = [];
   Session.Trial.Params.Incremental = [];
   
    %%
    
    %%
    Session.Data = s.data;                % just copy this 
    
    % for session level data, look in top level of data struct and in .info
    for FieldName = fieldnames(Session)';              % get field names
        % find intersecting field names either with first letter cap or lower
       fName = [intersect(FieldName,fieldnames(s))', intersect(lower(FieldName),fieldnames(s))'];    % find any in top level of saved session struct
       if(~isempty(fName))
           if(~isstruct(s.(char(fName))))                    % if it is not a struct
               fName = fName{1};
               capfName = [upper(fName(1)),fName(2:end)];    %   make first letter capitalized
               Session.(capfName) = s.(fName);               %   copy
           end
       end
       
       fName = [intersect(FieldName,fieldnames(s.info))', intersect(lower(FieldName),fieldnames(s.info))'];    % find any in top level of saved session struct
       if(~isempty(fName))
           if(~isstruct(s.info.(char(fName))))                    % if it is not a struct
               fName = fName{1};
               capfName = [upper(fName(1)),fName(2:end)];    %   make first letter capitalized
               Session.(capfName) = s.info.(fName);              %   copy
           end
       end
    end
    % for trial level stuff, go trial by trial and pull up data from .params when needed
   % FieldName = fieldnames(Session.Trial)
   OKtrialIndx = find(cellfun(@isempty,{s.tr.trial_start_ts})==0); % find trials where start timestamp is not empty: omits empty trials
   
   for FieldName = fieldnames(Session.Trial)'
        %look in tr.params
        fName = [intersect(FieldName,fieldnames(s.tr(1).params))', intersect(lower(FieldName),fieldnames(s.tr(1).params))'];    % find any in top level of saved session struct
        if(~isempty(fName))
            if(~isstruct(s.tr(1).params.(char(fName))))                       % if it is not a struct
                fName = fName{1};
                capfName = [upper(fName(1)),fName(2:end)];  %   make first letter capitalized
                for i=OKtrialIndx
                    Session.Trial(i).(capfName) = s.tr(i).params.(fName);               %   copy
                end
            end
        end
        
        %look in tr
        fName = [intersect(FieldName,fieldnames(s.tr))', intersect(lower(FieldName),fieldnames(s.tr))'];    % find any in top level of saved session struct
        if(~isempty(fName))
            if(~isstruct(s.tr(1).(char(fName))))                       % if it is not a struct
                fName = fName{1};
                capfName = [upper(fName(1)),fName(2:end)];  %   make first letter capitalized
                for i=OKtrialIndx
                    Session.Trial(i).(capfName) = s.tr(i).(fName);               %   copy
                end
            end
        end
       
   end
    
    for FieldName = fieldnames(Session.Trial(1).Params)'
        fName = [intersect(FieldName,fieldnames(s.tr(1).params))', intersect(lower(FieldName),fieldnames(s.tr(1).params))'];    % find any in top level of saved session struct
        if(~isempty(fName))
            if(~isstruct(s.tr(1).params.(char(fName))))                       % if it is not a struct
                fName = fName{1};
                capfName = [upper(fName(1)),fName(2:end)];  %   make first letter capitalized
                for i=OKtrialIndx
                    if(isfield(s.tr(i).params, fName))
                        Session.Trial(i).Params.(capfName) = s.tr(i).params.(fName);               %   copy
                    else
                        Session.Trial(i).Params.(capfName) = NaN;
                    end
                end
            end
        end
        
    end
    'hi'
end
        
    
