function   [data_out,out]=mdlPredict(data_in,par)
% function [data_out,out]=mdlPredict(data_in,par)

execinfo    = par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ', mfilename); end

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
try
    labs                = [data_in.trialType]';
    % Confusion matrix
    Cmatrxix            = confusionmat(labs, y_pred);
    % kappa value
    kappaValue          = kappaModel(Cmatrxix);
    % accuracy
    Accuracy            = sum(y_pred == labs)/length(labs)*100;
    Accuracy_class(1,:) = accuracy4classes(labs,y_pred);

catch
end
data_out                = data_in;
nTrials                 = length(data_in);
xfld                    = 'time';
for iTrial=1:nTrials
    timeField                               = data_in(iTrial).([xfld InField]);
    % prediction
    data_out(iTrial).(OutField)             = y_pred(iTrial);
    data_out(iTrial).([xfld OutField])      = timeField(end);
    % probabilities
    data_out(iTrial).(ProbField)            = probs(iTrial,:)';
    data_out(iTrial).([xfld ProbField])     = timeField(end);
    % success 
    try
        data_out(iTrial).(SuccessField)         = y_pred(iTrial)==labs(iTrial);
        data_out(iTrial).([xfld SuccessField])  = timeField(end);
    catch
    end
end

if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end
%% output save
try
    out.kappaValue      = kappaValue;
    out.Accuracy        = Accuracy;
    out.Accuracy_class  = Accuracy_class;
    out.Cmatrxix        = Cmatrxix;
catch
end