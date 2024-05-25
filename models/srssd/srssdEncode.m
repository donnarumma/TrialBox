function   [data_trials, out]=srssdEncode(data_trials,par)
% function [data_trials, out]=srssdEncode(data_trials,par)
execinfo    = par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ',mfilename); end
nTrials     = length(data_trials);

InField     = par.InField;
OutField    = par.OutField;
xfld        = par.xfld;
Wsrssd      = par.Wsrssd;                   % nAtoms x nChannels*nTimes
X3d         = cat(3,data_trials.(InField)); % nChannels x nTimes x nTrials.
X3d         = permute(X3d,[3,1,2]);         % nTrials x nChannels x nTimes

% set each row as a "synergy" -> [x_1(1), ... x_nChannels(1), x_1(2), ... x_nChannels(2), ... , x_1(nTimes), ... x_nChannels(nTimes),    
Xdata       = reshape(X3d,size(X3d,1),size(X3d,2)*size(X3d,3)); % nTrials x nChannels*nTimes

% Zdata       = Xdata*Wsrssd'; % nTrials x nAtoms
% Zdata       = Zdata';
Zdata       = Wsrssd*Xdata';                % nAtoms x nTrials

for iTrial=1:nTrials
    data_trials(iTrial).(OutField)          = Zdata(:,iTrial);
    time                                    = data_trials(iTrial).([xfld InField]);
    data_trials(iTrial).([xfld OutField])   = time(end);
end

if ~isempty(execinfo); out.exectime=toc(t); fprintf('Time Elapsed: %.2f s\n',out.exectime); end
out.numComponents =size(Wsrssd,1);