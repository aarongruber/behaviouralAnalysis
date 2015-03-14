% 2012-05-17. Leonardo Molina.
% Last modification: 2012-05-25.
% MERGE_SESSIONS merge trials from two or more behavioral files.
%
% session = MERGE_SESSIONS(file1, file2, ...) merges trials from file1, file2, ... into session.
% MERGE_SESSIONS(output_file, file1, file2, ...) merges trials from file1, file2, ... and saves it to output_file.

% Last edit: 2012-05-22.
function session = merge_sessions(varargin)
    % Adjust if necesary. Current settings correspond to behavioral_program_v19.
    fields_to_shift = {'trial_start_ts'; 'L_light_on_ts'; 'L_light_off_ts'; 'R_light_on_ts'; 'R_light_off_ts'; 'odor_port_on_ts'; 'odor_port_off_ts'; 'odor1_on_ts'; 'odor1_off_ts'; 'odor2_on_ts'; 'odor2_off_ts'; 'tone_on_ts'; 'tone_off_ts'; 'well_on_ts'; 'well_off_ts'; 'L_reward_ts'; 'L_quinine_ts'; 'R_reward_ts'; 'R_quinine_ts'; 'L_lick_ts'; 'R_lick_ts'; 'stimulation_on_ts=integer'; 'stimulation_off_ts'; 'house_light_on_ts'; 'house_light_off_ts'; 'odor_port_all_on_ts'; 'odor_port_all_off_ts'; 'R_well_all_on_ts'; 'R_well_all_off_ts'; 'L_well_all_on_ts'; 'L_well_all_off_ts'; 'well_vacuum_on_ts'; 'well_vacuum_off_ts'};
    
    % Load files.
    files = varargin(-nargout+2:end); if numel(files) == 1, return; end
    data = struct;
    for i = 1:numel(files)
        file = files{i};
        try
            load(file); if isempty(session.tr(end).trial_start_ts), session.tr(end) = []; end
            for field = fieldnames(session)'
                data.(sprintf('S%03i', i)).(field{:}) = session.(field{:});
            end
        catch
            fprintf(2, 'Couldn''t load %s... skipping!\n', file);
        end
    end
    clear session;
    session_names = fieldnames(data)'; ns = numel(session_names);
    if ns < 2, fprintf(2, 'Not enough input files!\n'); return; end
    
    % Field names common to all input files.
    common_fields = fieldnames(data.(session_names{1}).tr);
    intertrials = [];
    for i = 2:ns
        common_fields = intersect(common_fields, fieldnames(data.(session_names{i}).tr));
        intertrials = [intertrials data.(session_names{i}).tr.trial_start_ts];
    end
    fields_to_shift = fields_to_shift(ismember(fields_to_shift, common_fields));
    trial_delay = nanmean(intertrials);
        
    % Save information (i.e. data and info fields) from input file #1 only, except for the tr field.
    session = data.(session_names{1});
    session.tr = cell2struct(cell(numel(common_fields),1), common_fields); session.tr(1) = [];
    % Backup all input sessions in the new session.
    session.merged = data;
        
    % Merge.
    shift = 0;
    for i = 1:ns
        % Only keep fields common to all sessions.
        tr = data.(session_names{i}).tr;
        tr_fields = fieldnames(tr);
        rm_fields = common_fields(~ismember(common_fields, tr_fields));
        session.tr = rmfield(session.tr, rm_fields);
        % Report removed fields if any.
        if ~isempty(rm_fields)
            fprintf(2, 'Fields {');
            for rm_field = rm_fields', fprintf(2, ' %s ', rm_field{:}); end
            fprintf(2, '} from input file #%i were not included.\n', i);
        end        
        % Add temporal shift.
        nt = numel(tr);
        for j = 1:numel(fields_to_shift)
            for k = 1:nt
                tr(k).(fields_to_shift{j}) = tr(k).(fields_to_shift{j}) + shift;
            end
        end
        % Actual merging.
        session.tr = [session.tr tr];
        % Update shift.
        if ~isempty(tr(end).trial_start_ts), shift = tr(end).trial_start_ts + trial_delay;
        else shift = shift + trial_delay;
        end
    end
    % Save if user provided an output file.
    if nargout == 0
        save(varargin{1}, 'session');
    end
end