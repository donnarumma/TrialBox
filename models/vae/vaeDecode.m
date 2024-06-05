function   [data_trials, out]=vaeDecode(data_trials,par)
% function [data_trials, out]=vaeDecode(data_trials,par)
execinfo    = par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ',mfilename); end
nTrials     = length(data_trials);

InField     = par.InField;
OutField    = par.OutField;
xfld        = par.xfld;

Z2d         = cat(2,data_trials.(InField))'; % manifold-> nTrials x nComponents  
if canUseGPU
    fprintf('using GPU ');
    Z2d = gpuArray(Z2d);
end
% decode manifold
X3d         = squeeze(predict(par.netD,Z2d));
if isgpuarray(X3d)
    X3d=gather(X3d);
end
for iTrial=1:nTrials
    data_trials(iTrial).(OutField)          = X3d(:,:,iTrial);
    time                                    = data_trials(iTrial).([xfld InField]);
    data_trials(iTrial).([xfld OutField])   = time(end);
end

if ~isempty(execinfo); out.exectime=toc(t); fprintf('Time Elapsed: %.2f s\n',out.exectime); end

end
