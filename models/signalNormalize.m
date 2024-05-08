function   [data_trials, out] = signalNormalize(data_trials,par)
% function [data_trials, out] = signalNormalize(data_trials,par)
execinfo=par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ',mfilename); end

InField     = par.InField;
OutField    = par.OutField;


nTrials =length(data_trials);
for iTrial = 1:nTrials
    signalIn                        = data_trials(iTrial).(InField);
    % rescale       = max(real(signalIn(:)))/rsfactor;       % rescale
    % rescale       = 5*max(real(signalIn(:)));
    rescale                         = sqrt(sum(abs(signalIn(:)).^2));
    % rescale         = 4*std(signalIn(:));
    % rescale       = (max(abs(signalIn(:)) - min(abs(signalIn(:)))))/4;

    data_trials(iTrial).(OutField)  = signalIn/rescale;
    data_trials(iTrial).(['time' OutField])=data_trials(iTrial).(['time' InField]);
end

if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end
out.rescale = rescale;