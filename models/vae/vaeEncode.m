function   [data_trials, out]=vaeEncode(data_trials,par)
% function [data_trials, out]=vaeEncode(data_trials,par)
execinfo    = par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ',mfilename); end
nTrials     = length(data_trials);

InField     = par.InField;
OutField    = par.OutField;
xfld        = par.xfld;

X4d = cat(4,data_trials.(InField)); % nChannels x nTimes x 1 x nTrials
if canUseGPU
    fprintf('using GPU ');
    X4d = gpuArray(X4d);
end

Z2d = predict(par.netE,X4d)';       % nComponents x 1
for iTrial=1:nTrials
    data_trials(iTrial).(OutField)          = Z2d(:,iTrial);
   
    time                                    = data_trials(iTrial).([xfld InField]);
    data_trials(iTrial).([xfld OutField])   = time(end);
end

if ~isempty(execinfo); out.exectime=toc(t); fprintf('Time Elapsed: %.2f s\n',out.exectime); end

end
