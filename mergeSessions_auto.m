function mergeSessions_auto(dirpath)
%
% Detects multiple files from same animal in a given session (in directory 'dir')  
% merges session data. If no direcotry is given, then program uses pwd.
% original files are deleted, so run only on a copy of data.
%
% Aaron Gruber    2012_5_31


if(nargin<1)
    dirpath = pwd;
end

d = dir([dirpath,filesep,'*.mat']);    % use platform independent filsep command
rd=cellfun(@(x) x(1:16), {d.name}, 'Uniformoutput', false); % get rat/date part
rd_unique = unique(rd);

for i=1:numel(rd_unique)
   indx = ismember(rd,rd_unique(i)); % find indicies of unique rat/date in original data
   if(sum(indx)>1)                   % if more than one, combine
       oldfn = {d(indx).name};         % get full name
       oldfntxt = ['''',oldfn{1}];
       for j=2:numel(oldfn)
          oldfntxt = [oldfntxt,''', ''',oldfn{j}];   % make text string for function 'merge_session' 
       end
       oldfntxt =  [oldfntxt,''''];
       textcmd = ['merge_sessions([''',rd_unique{i},''',','''_merged.mat'']',',',oldfntxt,')'];
       eval(textcmd)
       disp(['Merged :', oldfntxt])
       delete(oldfn{:})                 % delete the source files
   end
end

remove_incomplete_trials_from_file % remove empty trials from all files.
