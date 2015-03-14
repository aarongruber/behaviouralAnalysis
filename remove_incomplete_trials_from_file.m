function remove_incomplete_trials_from_file(fd)
%
% Removes incomplete trials from BEH.mat files in directory 'fd'
% (optional). Otherwise looks in all files in pwd.
%
%  The behavioural control system has been crashing - the resultant BEH
%  file can contain incomplete data; e.g. some of the trials have only
%  partial data, which causes the analysis program to crash.

if(nargin<1)
    fd = pwd;
end

files = dir([fd,'\*.mat']);

for fn = {files.name}
   load(fn{:});
   OKIndx = ones(numel(session.tr),1);
   for i=1:numel(session.tr)
      if(~isfield(session.tr(i).params, 'debounce_down_odor_port'))
          OKIndx(i)=0;
      end
   end
   if(sum(OKIndx)<numel(OKIndx))
     session.tr = session.tr(logical(OKIndx));
     save(fn{:},'session') 
   end
   clear session
end