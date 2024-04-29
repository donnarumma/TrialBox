function   [data_trials,out]=predictNBPW(data_trials,par)
% function [data_trials,out]=predictNBPW(data_trials,par)
execinfo    = par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ', mfilename); end

labs                    = [data_trials.trialType]';
InField                 = par.InField;
OutField                = par.OutField;
SuccessField            = par.SuccessField;

%% Train 
y_pred                  = par.labs_pred; 
Accuracy                = sum(y_pred == labs)/length(labs)*100;
Accuracy_class(1,:)     = accuracy4classes(labs,y_pred);

nTrials     = length(data_trials);
for iTrial=1:nTrials
    data_trials(iTrial).(OutField)          = y_pred(iTrial);
    data_trials(iTrial).(SuccessField)      = y_pred(iTrial)==labs(iTrial);
    timeField                               = data_trials(iTrial).(['time' InField]);
    data_trials(iTrial).(['time' SuccessField])=timeField(end);
    data_trials(iTrial).(['time' OutField])=timeField(end);
end

if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end
%% output save
out.Accuracy        = Accuracy;
out.Accuracy_class  = Accuracy_class;
