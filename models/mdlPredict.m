function   [data_out,out]=mdlPredict(data_in,par)
% function [data_out,out]=mdlPredict(data_in,par)

execinfo    = par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ', mfilename); end

labs                    = [data_in.trialType]';
InField                 = par.InField;
OutField                = par.OutField;
SuccessField            = par.SuccessField;
ProbField               = par.ProbField;
% nChannels x nTimes x nTrials
data_3d                 = cat(3,data_in.(InField));
% nTrials x nChannels x nTimes 
data_3d                 = permute(data_3d,[3,1,2]);  
% nTrials x nChannels * nTimes
data_3d                 = reshape(data_3d,size(data_3d,1),size(data_3d,2)*size(data_3d,3)); 

mdl                     = par.mdl;
%% Predict 
[y_pred,probs]          = predict(mdl, data_3d); % nTrials x nClasses
Accuracy                = sum(y_pred == labs)/length(labs)*100;
Accuracy_class(1,:)     = accuracy4classes(labs,y_pred);

data_out                = data_in;
nTrials                 = length(data_in);
xfld                    = 'time';
for iTrial=1:nTrials
    data_out(iTrial).(OutField)          = y_pred(iTrial);
    data_out(iTrial).(SuccessField)      = y_pred(iTrial)==labs(iTrial);
    data_out(iTrial).(ProbField)         = probs(iTrial,:)';
    timeField                            = data_in(iTrial).([xfld InField]);
    data_out(iTrial).([xfld SuccessField]) = timeField(end);
    data_out(iTrial).([xfld OutField])  = timeField(end);
    data_out(iTrial).([xfld ProbField]) = timeField(end);
end

if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end
%% output save
out.Accuracy        = Accuracy;
out.Accuracy_class  = Accuracy_class;
