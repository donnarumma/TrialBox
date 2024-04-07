function   [data_trials, out]=pcaDecode(data_trials,par)
% function [data_trials, out]=pcaDecode(data_trials,par)
% pca executed per class
execinfo    = par.exec;
if ~isempty(execinfo); t=tic; end

InField     = par.InField;
OutField    = par.OutField;
Wpca        = par.Wpca;

nTrials     = length(data_trials);
nChannels   = size(par.mu,1);
[~,nTimes]  = size(data_trials(1).(InField));

Z_data      = [data_trials.(InField)];          % data are nChannels x nTimes * nTrials 
murep       = repmat(par.mu,1,size(Z_data,2));
X_data_rec  = Wpca * Z_data + murep;  % nChannels x nTimes*nTrials

X_data_rec  = reshape(X_data_rec,nChannels,nTimes,nTrials);   % nChannels x nTimes x nTrials trials in test

xfld    = 'time';
for iTrial=1:nTrials
    data_trials(iTrial).(OutField)          = X_data_rec(:,:,iTrial);
    data_trials(iTrial).([xfld OutField])   = data_trials(iTrial).([xfld InField]);
end


if ~isempty(execinfo); out.exectime=toc(t); fprintf('Function: %s | Time Elapsed: %.2f s\n',mfilename,out.exectime); end