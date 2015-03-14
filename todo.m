% to do
% 
%
% Neural analysis:
%%  overall organization %%
%   - object for units (tet#, depth, struct, grade)
%   - object for lfp (tet#, depth, struct, artifact times)
%   - vid
%   - FSCV
%   - emg
%
% 
% 
%
%% General program & analysis
%
%
% GUI for editing metadata - have vid tracker for move no move?
%        - epochs {rest1, task, ....} -> add to trial 
%
% issue warning if rat_id seems to be misslabled (e.g. unequal names)
%
% improve syntax for getTrials (make like dataset array): always do by
% session to avoid cross session computations.
%
% - bar and error plot: make group enabled like boxplot
%       * write new data parser- improve speed
%
% anova - repeats or not. - make sure that tags are in the correct format
%   - check for sphericity?
%   - check for RepeatMeasures in Trial tag data - burden off users
%
%  computeTrialByTrial -> 
%    option for session level -> rebuild session & compute
%    set range of indicies
%    also take function handle
%
%  computeOverSession (compute on trials object better b/c allows filtering)
%
%  determine if debounce time needs to be added to odor port & wells
%          
%   
% - analysis of response choices (RRL...)
%
% - plots & stats on vectors per session 
%       * have index for vector as tag -> parse & stats like others?
%
% - add property for task type (enginge & mode), or knickname
%
%  - add property (and data times) for beh epoch: rest, task, rest2....
%
% - make 'trials' based on video tracker? (move, rest periods?)
%
% - replace registry with a function that makes registry out of already
% loaded animals (load sessions without predefining registry). Get rid of
% registry at time of construction.
%
%%%%% plots
% ----- box plot -----
%   - group bars
%   - when date is first, just use one per grouping
%   - color for different groups - put in dataset array so all plots are
%   consistent
%
% --- bar&line ----
%   - make a data parser from grpstats 
%   - make plot functs that take dataset arrays & group tags
%   - method to add values & show stats
%
%   - multiple dimensions - not just 2: e.g. group, trial & block (avg over rat)
%   - add easy method for writing rat_id's, data files, group assign on plots
%
%%%%%%% unit analysis %%%%
% - if units, then correct ts in behavioral structure: seperate program
% that does this before this analysis?
%      * code checks for 'events' in same folder & runs registration prog.
% - check that BEH timestamps are in register (range of) neural data
% - link in session & have registry
% - spike object with: ses link, tet#, position, structure,
%       * for analysis group into session & tets for xcorr (same parse
%       funct as for behavior)
%
%

 %% wish list     
 %    - parser for logic calls to merge function calls in 'findTrials' &
 %    'getTrials': use matlab style synax: such as ('Property', val) and
 %    (property = val), (property < val), (property ~= val), (property1 =
 %    val1 & property2 < val). *** should have optional switch '-session' for
 %    doing calc only within session.
 %
 %    - could have program dynamically add properties for each property of
 %    the object being read in. i.e. don't require specification of all properties
 %    in the class. Would add flexibility for other behaviors.
 %
 %   can use dataset array with many non timestamp data as 'nominal'. Would
 %   make it much faster & syntax easier (no need for cellfun) & numerical
 %   & text can be used interchangeably; One field would be pointer to
 %   session. 
 %      - advantages of dataset: easier search syntax,
 %      - dissadvantages: can't copy by reference, not efficient to use
 %                  trials as basic unit-would be session instead
 %      - have custom class that has dataset in it for trials-> allows
 %      encapsulation of methods & protects data
 %
 %  - like a way to compute session-wise stats from the main function (%
 %  responses to L, ....)
 %
 %  - way for users to select rows in spreadsheet for choosing input files
 %  (they just ckeck a column)
 %
 %  - click on a plot and the Anova pops up.
 
 
 %% - NOTES
 % If changed to use strings->nominal; then can use complex logic synax to
 % find trials: NewTrials = trials(tr.prop1=='r' & tr(1:end-2).prop2-tr(2:end).prop3 < 5);
 % combine trials NewTrials = intersect(trials1.Trial_id, trials2.Trial_id)
 
 
 
 %% TO DO in behavioral program
 % - add arduinoID (or com as lazy way to do it)
 % - save files at BEHv1.mat
 % - take out progressive delay in all progs except ones that say 'shaping'
 % - remove delay & odor from matching pennies protocol.
 % - use capital letters for field names
 % - move params to top level
 % - no empty entries - make them nan (either numeric or string 'nan'
 % consistemt with type for property)
 % - save probabilities of reward & odours
 % *** rectify mismatch in param names: 'L_bolus_delay', vs 'R_reward_delay'
 
 
 %%
 
 
 
 % Matlab notes
 evalin
 username
 date
 mfunction name
 .git version number
 
 for over cell (has to be 1xn)
 git
 handles
 painters
 outputting progress and overwriting
 
 boxplot (tough sob)
 
 Arduino integration
 ANOVA - repeated measures