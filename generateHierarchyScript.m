
%% New constructor
%
% This is an example of how to use the new framework of handles class to
% facilitate handling data structures. 
%
% the strategy is to 
%   1) Specify directories that contain matlab BEH sessions and neural data
%   2) Make a registry of all animals that is dynamically populated as each
%      BEH data file is read in; Animal IDs are read in from the session files 
%      and when unique IDs are enountered, new animals are created in the registry 
%   3) Create session objects that will be linked to animals listed in
%      the registry. These have session data and a struct of all trials.
%   4) One of the properties of Trials is a link to the parent Session,
%      which can be used to look up any session or animal data for analysis 

% specify some directories
sessiondir1 = 'C:\Users\aaron.gruber\work\projects\dataStructDev\RatData\ControlRat1\2011_06_13';
sessiondir2 = 'C:\Users\aaron.gruber\work\projects\dataStructDev\RatData\ControlRat1\2011_06_16';
sessiondir3 = 'C:\Users\aaron.gruber\work\projects\dataStructDev\RatData\ControlRat2\2011_04_14';

registry = AnimalRegistry; % this object is a dynamic directory of all rats read in, which come from session data loaded in matlab

% These should be added to animal 1 (because the sesssion.mat data indicate
% a common rat ID.
s1 = Session(sessiondir1,registry);
s2 = Session(sessiondir2,registry);

% This should create a new animal, since we have no animals with the
% corresponding ID
s3 = Session(sessiondir3,registry);

% examples of searching for specific trials
session1trials = s1.getTrials('odor1_id',2)
group1trials = s3.getTrials('odor1_id',2)

% how to use this framework
registry.Animals % is the top level object with instances for each rat
registry.Animals(1) % is the first animal
registry.Animals(1).Sessions(1).getTrials('odor1_id',2) % get trils matching criteria from Session 1
registry.Animals(1).Sessions(:).getTrials('odor1_id',2) % get trials matching criteria from all sessions in Animal 1
registry.Animals.getTrials('odor1_id',2) % get trils from all animals and sessions matching criteria
a1 = registry.Animals(1) % assign a knickname to the handle
leftTR = a1(1).Sessions(:).getTrials('odor1_id',2)  % get trials for odor 2 in all sessions
leftTR(1).Session.date % is the date from the first trial, which is from session 1
leftTR(end).Session.date % is the date from teh last trial, which is session 2 ******BECAUSE EACH TRIAL HAS A REFERENCE TO THE SESSION, WE HAVE ACCESS TO ALL THE DATA FOR THE SESSION IN THE DERIVED OBJECTS THAT MEET SOME CRITERION*****


%% another way if assigning groups to animals at object creation
% g = Group;
% g = [g Group];
% a(1) = Animal(registry,g(1));
% a(1).ID = '2011_001_FB';
% a(2) = Animal(registry,g(1));

%%%%
% load TestTrial
% g = Group;
% g = [g Group];
% a1 = Animal(g(1));
% a2 = Animal(g(1));
% s1 = Session(a1,tr);
% s2 = Session(a2,tr);
% 
% session1trials = s1.getTrials('odorID',5)
% group1trials = g(1).getTrials('odorID',5)