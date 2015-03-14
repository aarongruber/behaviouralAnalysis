function plot_allBEHinFolder
%
% run session-wise behavioural data plotter
% 'plot_singleSessionBinaryChoice' on all .mat files in directory

d = dir('*.mat');

for k=1:length(d)
   load(d(k).name);
   plot_singleSessionBinaryChoice(session);
   clear session
end