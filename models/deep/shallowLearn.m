function net = shallowLearn(data_trials,par)
InField             = par.InField;
numEpochs           = par.numEpochs;
log_train           = par.log_train;
log_valid           = par.log_valid;
miniBatchSize       = par.miniBatchSize;
validationFrequency = par.validationFrequency;
% log_test            = par.log_test;
% numHiddenUnits      = par.numHiddenUnits;
verbose             = 1;

labels              = categorical([data_trials.trialType])';
numClasses          = length(unique(labels));

X3d                     = cat(3,data_trials.(InField));                     % nCells x nTimes x nTrials
% X3d         = permute(X3d,[3,1,2]);                             % nTrials x nCells x nTimes
% each row is a "synergy" -> [x_1(1), ... x_nChannels(1), x_1(2), ... x_nChannels(2), ... , x_1(nTimes), ... x_nChannels(nTimes),    

X4d                     = reshape(X3d,size(X3d,1),size(X3d,2),1,size(X3d,3)); % nTrials x nChannels x 1 x nTimes
[nCells, nTimes, nChannels, nTrials] = size(X4d);

XTrain              = X4d(:,:,:,log_train);
TTrain              = labels(log_train);

XValid              = X4d(:,:,:,log_valid);
TValid              = labels(log_valid);

% Xdata =cell(length(data_trials),1);
% for it=1:length(data_trials)
%     Xdata{it}=data_trials(it).(InField)'; % transpose in times x nvar
% end

% XTrain              = Xdata(log_train);
% TTrain              = labels(log_train);

% XValid              = Xdata(log_valid);
% TValid              = labels(log_valid);

% XTest               = Xdata(log_test);
% TTest               = labels(log_test);

% inputSize           = size(Xdata{1},1);
% nChannels           = 1;
% [nCells,nTimes]     = size(XTrain{1});
imageSize               = [nCells, nTimes, nChannels];

% layersLSTM = [ ...
%     sequenceInputLayer(nCells)
%     bilstmLayer(par.numLatentChannels,OutputMode="last")
%     fullyConnectedLayer(numClasses)
%     softmaxLayer
%     classificationLayer];
% 
layers = [  imageInputLayer(imageSize,Normalization="none"),                                    ...                              ...
            convolution2dLayer(par.kernsize,  par.numLatentChannels,Padding="same",Stride=2),   ...
            reluLayer,                                                                          ...
            fullyConnectedLayer(numClasses),                                                    ...
            softmaxLayer,                                                                       ...
            classificationLayer];   

% layers = [ ...
%     imageInputLayer(imageSize)
%     convolution2dLayer(par.kernsize,  par.numLatentChannels)
%     reluLayer
%     maxPooling2dLayer(2,'Stride',2)
%     fullyConnectedLayer(numClasses)
%     softmaxLayer
%     classificationLayer];

if isempty (XValid)
    options = trainingOptions("adam", ...
    ExecutionEnvironment="cpu", ...
    GradientThreshold=1, ...
    MaxEpochs=numEpochs, ... 
    MiniBatchSize=miniBatchSize, ...
    SequenceLength="longest", ... Shuffle="never", ...
    Verbose=verbose, ... 
    Plots="training-progress");
else
    options = trainingOptions("adam", ...
        ExecutionEnvironment="cpu", ...
        GradientThreshold=1, ...
        MaxEpochs=numEpochs, ... 
        MiniBatchSize=miniBatchSize, ...
        SequenceLength="longest", ... Shuffle="never", ...
        Verbose=verbose, ... 
        Plots="training-progress", ...
        ValidationData={XValid,TValid}, ...
        ValidationFrequency=validationFrequency);
end
net                     = trainNetwork(XTrain,TTrain,layers,options);