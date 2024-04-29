function   [data_trials,out]=predictKAPPA(data_trials,par)
% function [data_trials,out]=predictKAPPA(data_trials,par)
execinfo    = par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ', mfilename); end
labs                    = [data_trials.trialType]';
InField                 = par.InField;
OutField                = par.OutField;
SuccessField            = par.SuccessField;
% nChannels x nTimes x nTrials
data_3d                 = cat(3,data_trials.(InField));
% nTrials x nChannels x nTimes 
data_3d                 = permute(data_3d,[3,1,2]);  
% nTrials x nChannels * nTimes
data_3d                 = reshape(data_3d,size(data_3d,1),size(data_3d,2)*size(data_3d,3)); 

mdl                     = par.mdl;
%% Train 
y_pred                  = predict(mdl, data_3d); 

Cmatrxix = confusionmat(labs, y_pred);
% kappa value
kappavalue = kappaModel(Cmatrxix);

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
out.kappaValue     = kappavalue;