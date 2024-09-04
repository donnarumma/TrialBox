% function TEST_SHALLOW_LEARN


pms.ifplot      = true;
pms.root_dir    = '~/TESTS/MNIST/deep/';
pms.model_file  = false;
pms.numEpochs   = 150;

ifplot                      = pms.ifplot;
root_dir                    = pms.root_dir;
model_file                  = pms.model_file;
numEpochs                   = pms.numEpochs;
ref                         = datetime('now','Format','yyyyMMddHHmmss');
save_dir                    = sprintf('%s%s/',root_dir,ref);
description                 = mfilename;
out.ref                     = ref; 
out.root_dir                = root_dir;
%--------------------------------------------------------------------------

%% load MNIST Dataset sample in Matlab
[XTrain,YTrain,anglesTrain] = digitTrain4DArrayData;
[XTest, YTest, anglesTest ] = digitTest4DArrayData;

%% arrange data_trials_train and data_trials_test 
% nChannels x nTimes x nTrials 
X3dTrain                    = squeeze(XTrain);
X3dTest                     = squeeze(XTest);
nTimes                      = size(X3dTrain,2); % should be the same as test
time                        = 1:nTimes;
xfld                        = 'time';
InField                     = 'MNIST'; 
% train
nTrials_train               = size(X3dTrain,3);         
labels                      = unique(YTrain);
data_trials_train           = struct;
for iTrial=1:nTrials_train
    data_trials_train(iTrial).(InField)       = X3dTrain(:,:,iTrial);
    data_trials_train(iTrial).([xfld InField])= time;
    data_trials_train(iTrial).trialType       = find(ismember(labels,YTrain(iTrial)))-1;
    data_trials_train(iTrial).trialName       = num2str(data_trials_train(iTrial).trialType);
    data_trials_train(iTrial).behavior        = anglesTrain(iTrial);
    data_trials_train(iTrial).train           = true;
    data_trials_train(iTrial).test            = false;
end

% test
data_trials_test            = struct;
nTrials_test                = size(X3dTest,3);         
labels                      = unique(YTest);
for iTrial=1:nTrials_test
    data_trials_test(iTrial).(InField)       = X3dTest(:,:,iTrial);
    data_trials_test(iTrial).([xfld InField])= time;
    data_trials_test(iTrial).trialType       = find(ismember(labels,YTest(iTrial)))-1;
    data_trials_test(iTrial).trialName       = num2str(data_trials_test(iTrial).trialType);
    data_trials_test(iTrial).behavior        = anglesTrain(iTrial);
    data_trials_test(iTrial).train           = true;
    data_trials_test(iTrial).test            = false;
end

%% Learn Shallow Conv Net
nTrials                         = length(data_trials_train);
seq                             = randperm(nTrials);
par.shallow                     = shallowLearnParams();%vaeModelParams();
%no validation
% par.shallow.log_valid         = [];
% par.shallow.log_train         = 1:5000;
par.shallow.log_valid           = seq(1:1000);          
par.shallow.log_train           = seq(1001:end);
par.shallow.numLatentChannels   = 32; % dimension of the manifold
par.shallow.numEpochs           = 150;
par.shallow.InField             = 'MNIST';

net                             = shallowLearn(data_trials_train,par.shallow);

%% Accuracy on Training
labels_trainvalid       = categorical([data_trials_train.trialType])';
X3d                     = cat(3,data_trials_train.(InField));                     % nCells x nTimes x nTrials
X4d                     = reshape(X3d,size(X3d,1),size(X3d,2),1,size(X3d,3)); % nTrials x nChannels x 1 x nTimes
XTrain                  = X4d(:,:,:,par.shallow.log_train);
TTrain                  = labels_trainvalid(par.shallow.log_train);
YTrain                  = classify(net,XTrain);
accTrain                = mean(YTrain==TTrain);
fprintf('Accuracy on Train: %g\n',accTrain)

%% Accuracy on Validation
if ~isempty(par.shallow.log_valid)
    XValid                  = X4d(:,:,:,par.shallow.log_valid);
    TValid                  = labels_trainvalid(par.shallow.log_valid);
    YValid                  = classify(net,XValid);
    accValid                = mean(YValid==TValid);
    fprintf('Accuracy on Validation: %g\n',accValid)
end

%% Accuracy on Test
labels_test             = categorical([data_trials_test.trialType])';
X3dTest                 = cat(3,data_trials_test.(InField));                     % nCells x nTimes x nTrials
X4dTest                 = reshape(X3dTest,size(X3dTest,1),size(X3dTest,2),1,size(X3dTest,3)); % nTrials x nChannels x 1 x nTimes
XTest                   = X4dTest;
TTest                   = labels_test;
YTest                   = classify(net,XTest);
accTest                = mean(YTest==TTest);
fprintf('Accuracy on Test: %g\n',accTest)