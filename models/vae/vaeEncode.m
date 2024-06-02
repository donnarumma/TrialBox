function   [data_trials, out]=vaeEncode(data_trials,par)
% function [data_trials, out]=vaeEncode(data_trials,par)
execinfo    = par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ',mfilename); end
nTrials     = length(data_trials);

InField     = par.InField;
OutField    = par.OutField;
xfld        = par.xfld;
netE        = par.netE;                     % encoder net            
X3d         = cat(3,data_trials.(InField)); % nChannels x nTimes x nTrials (original space)
X4d         = reshape(X3d,size(X3d,1),size(X3d,2),1,size(X3d,3));   % nChannels x nTimes x 1 x nTrials.

dsE         = arrayDatastore(X4d,IterationDimension=4);
mbqE        = minibatchqueue(dsE,                                   ...
                            par.MBnumOutputs,                           ...
                            MiniBatchSize       = par.miniBatchSize,    ...
                            MiniBatchFcn        = @preprocessMiniBatch, ...
                            MiniBatchFormat     = "SSCB");

Zdata            = [];

% Loop over mini-batches.
while hasdata(mbqE)
    Xdata       = next(mbqE);

    % Forward through encoder.
    Zpred       = predict(netE,Xdata);

    % Extract and concatenate predictions.
    Zdata       = cat(2,Zdata,extractdata(Zpred));
end



for iTrial=1:nTrials
    data_trials(iTrial).(OutField)          = Zdata(:,iTrial);
    time                                    = data_trials(iTrial).([xfld InField]);
    data_trials(iTrial).([xfld OutField])   = time(end);
end

if ~isempty(execinfo); out.exectime=toc(t); fprintf('Time Elapsed: %.2f s\n',out.exectime); end
% out.numComponents =size(Wsrssd,1);

end

function   X = preprocessMiniBatch(dataX)
% function X = preprocessMiniBatch(dataX)
    % Concatenate.
    X                       = cat(4,dataX{:});
end