function   [data_trials, out]=vaeDecode(data_trials,par)
% function [data_trials, out]=vaeDecode(data_trials,par)
execinfo    = par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ',mfilename); end
nTrials     = length(data_trials);

InField     = par.InField;
OutField    = par.OutField;
xfld        = par.xfld;
% netD        = par.netD;                     % decoder net
% X3d         = cat(2,data_trials.(InField)); % nCells x nTrials (manifold)
% 
% dsD         = arrayDatastore(X3d,IterationDimension=2);
% mbqD        = minibatchqueue(dsD,                                   ...
%                             par.MBnumOutputs,                           ...
%                             MiniBatchSize       = par.miniBatchSize,    ...
%                             MiniBatchFcn        = @preprocessMiniBatch, ...
%                             MiniBatchFormat     = "CB");

% Xdata            = [];

% Loop over mini-batches.
% while hasdata(mbqD)
%     Zdata       = next(mbqD);
% 
%     % Forward through encoder.
%     Xpred       = predict(netD,Zdata);
% 
%     % Extract and concatenate predictions.
%     Xdata       = cat(4,Xdata,extractdata(Xpred));
% end

Z2d=cat(2,data_trials.(InField))';
X3d=squeeze(predict(par.netD,Z2d));
for iTrial=1:nTrials
    % data_trials(iTrial).(OutField)          = Xdata(:,:,iTrial);
    % Zdata                                   = data_trials(iTrial).(InField);
    % data_trials(iTrial).(OutField)          = predict(par.netD,Zdata');
    data_trials(iTrial).(OutField)          = X3d(:,:,iTrial);
    time                                    = data_trials(iTrial).([xfld InField]);
    data_trials(iTrial).([xfld OutField])   = time(end);
end

if ~isempty(execinfo); out.exectime=toc(t); fprintf('Time Elapsed: %.2f s\n',out.exectime); end
% out.numComponents =size(Wsrssd,1);

end

function   X = preprocessMiniBatch(dataX)
% function X = preprocessMiniBatch(dataX)
    % Concatenate.
    X                       = cat(2,dataX{:});
end