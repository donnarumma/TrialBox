function   [data_trials, out]=srssdDecode(data_trials,par)
% function [data_trials, out]=srssdDecode(data_trials,par)
% pca executed per class
execinfo    = par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ',mfilename); end

InField     = par.InField;
OutField    = par.OutField;
% nAtoms x nChannels*nTimes
Dsrssd      = par.Dsrssd;
nChannels   = par.nChannels;
nTimes      = par.nTimes;
nTrials     = length(data_trials);
% [nChannels,nTimes]  = size(data_trials(1).(InField));

% nAtoms x nTrials
Zdata       = cat(2,data_trials.(InField));                 
% Xdata_rec   = Wpca * Zdata + murep;                         % nChannels x nTimes*nTrials

% nTrials x nChannels*nTimes  
% Xdata_rec   = Zdata' * Dsrssd;                % 

% nChannels*nTimes x nTrials
Xdata       = Dsrssd'*Zdata;

% nChannels x nTimes x nTrials trials in test
X3d         = reshape(Xdata,nChannels,nTimes,nTrials);  
xfld        = 'time';
for iTrial=1:nTrials
    data_trials(iTrial).(OutField)          = X3d(:,:,iTrial);
    data_trials(iTrial).([xfld OutField])   = data_trials(iTrial).([xfld InField]);
end

if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end