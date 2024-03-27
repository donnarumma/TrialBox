function   [data_trials, out]=removeInactive(data_trials,par)
% function [data_trials, out]=removeInactive(data_trials,par)
execinfo=par.exec;
if ~isempty(execinfo); t=tic; end
InField  =par.InField;
OutField =par.OutField;
hasSpikesBool = (mean([data_trials.(InField)], 2) ~= 0);
for it=1:length(data_trials)
    data_trials(it).(OutField)          = data_trials(it).(InField)(hasSpikesBool,:);
    data_trials(it).(['time' OutField]) = data_trials(it).(['time' InField]); 
end
if ~isempty(execinfo); out.exectime=toc(t); fprintf('Function: %s | Time Elapsed: %.2f s\n',mfilename,out.exectime); end
out.hasSpikeBool=hasSpikesBool;
