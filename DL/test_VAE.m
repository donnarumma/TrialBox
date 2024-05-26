%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%% load MNIST in Matlab
[XTrain,YTrain,anglesTrain] = digitTrain4DArrayData;
[XTest, YTest, anglesTest ] = digitTest4DArrayData;
% nChannels x nTimes x nTrials 
X3dTrain                    = squeeze(XTrain);
X3dTest                     = squeeze(XTest);

[nChannels, nTimes, RGB, nTrials]  = size(XTrain); 
%%
par.VAE.numEpochs           = 150;
par.VAE.miniBatchSize       = 128;
par.VAE.learnRate           = 1e-3;
par.VAE.numLatentChannels   = 32;
imageSize                   = [nChannels nTimes RGB];

layersE = [
    imageInputLayer(imageSize,Normalization="none") ,   ...
    convolution2dLayer(3,32,Padding="same",Stride=2),   ...
    reluLayer,                                          ...                   
    convolution2dLayer(3,64,Padding="same",Stride=2),   ...
    reluLayer,                                          ...
    fullyConnectedLayer(2*par.VAE.numLatentChannels),   ...
    samplingLayer];

projectionSize              = [7 7 64];
numInputChannels            = imageSize(3);

layersD = [
    featureInputLayer(par.VAE.numLatentChannels)
    projectAndReshapeLayer(projectionSize)
    transposedConv2dLayer(3,64,Cropping="same",Stride=2)
    reluLayer
    transposedConv2dLayer(3,32,Cropping="same",Stride=2)
    reluLayer
    transposedConv2dLayer(3,numInputChannels,Cropping="same")
    sigmoidLayer];

netE        = dlnetwork(layersE);
netD        = dlnetwork(layersD);

dsTrain     = arrayDatastore(XTrain,IterationDimension=4);
numOutputs  = 1;

mbq = minibatchqueue(dsTrain,numOutputs,            ...
    MiniBatchSize       = par.VAE.miniBatchSize,    ...
    MiniBatchFcn        = @preprocessMiniBatch,     ...
    MiniBatchFormat     = "SSCB",                   ...
    PartialMiniBatch    = "discard");

trailingAvgE            = [];
trailingAvgSqE          = [];
trailingAvgD            = [];
trailingAvgSqD          = [];

numIterationsPerEpoch   = ceil(nTrials / par.VAE.miniBatchSize);
numIterations           = par.VAE.numEpochs * numIterationsPerEpoch;

monitor = trainingProgressMonitor(  ...
    Metrics             = "Loss",   ...
    Info                = "Epoch",  ...
    XLabel              = "Iteration");


% learning
epoch                   = 0;
iteration               = 0;
% Loop over epochs.
while epoch < par.VAE.numEpochs && ~monitor.Stop
    epoch = epoch + 1;

    % Shuffle data.
    shuffle(mbq);

    % Loop over mini-batches.
    while hasdata(mbq) && ~monitor.Stop
        iteration = iteration + 1;

        % Read mini-batch of data.
        X = next(mbq);

        % Evaluate loss and gradients.
        [loss,gradientsE,gradientsD] = dlfeval(@modelLoss,netE,netD,X);

        % Update learnable parameters.
        [netE,trailingAvgE,trailingAvgSqE]      = adamupdate(netE, ...
            gradientsE,trailingAvgE,trailingAvgSqE,iteration,par.VAE.learnRate);

        [netD, trailingAvgD, trailingAvgSqD]    = adamupdate(netD, ...
            gradientsD,trailingAvgD,trailingAvgSqD,iteration,par.VAE.learnRate);

        % Update the training progress monitor. 
        recordMetrics(monitor,iteration,Loss=loss);
        updateInfo(monitor,Epoch=epoch + " of " + par.VAE.numEpochs);
        monitor.Progress = 100*iteration/numIterations;
    end
end

%% train
dsTrain                 = arrayDatastore(XTrain,IterationDimension=4);
numOutputs              = 1;

mbqTrain = minibatchqueue(dsTrain,numOutputs,     ...
    MiniBatchSize       = par.VAE.miniBatchSize,...
    MiniBatchFcn        = @preprocessMiniBatch, ...
    MiniBatchFormat     = "SSCB");

XTrain_rec               = modelPredictions(netE,netD,mbqTrain);

%% test
dsTest                  = arrayDatastore(XTest,IterationDimension=4);
numOutputs              = 1;

mbqTest = minibatchqueue(dsTest,numOutputs,     ...
    MiniBatchSize       = par.VAE.miniBatchSize,...
    MiniBatchFcn        = @preprocessMiniBatch, ...
    MiniBatchFormat     = "SSCB");

XTest_rec               = modelPredictions(netE,netD,mbqTest);

%% 
errorsTrain             = errors(XTrain,XTrain_rec);
errorsTest              = errors(XTest,XTest_rec);
errorXTest              = errorsTrain.RMSE;
errorXTrain             = errorsTest.RMSE;
%% plot reconstructions
col             = [0.5,0.5,0.7]; % col borders
hmShow          = 49;
% train
numObservations = size(XTrain,4);
idx             = randperm(numObservations,hmShow);
ImgTrain        = squeeze (XTrain(:,:,:,idx));
ImgTrainRec     = squeeze (XTrain_rec(:,:,:,idx));
subplot(1,2,1)
montage (ImgTrain, 'BorderSize', [1,1], 'BackgroundColor', col);
title('Original');
subplot(1,2,2);
montage (ImgTrainRec, 'BorderSize', [1,1], 'BackgroundColor', col);
title('Reconstructed');
sgtitle(sprintf('Train (RMS: %.4f)',errorXTrain))
%% test
numObservations = size(XTest,4);
idx             = randperm(numObservations,hmShow);
ImgTest         = squeeze (XTest(:,:,:,idx));
ImgTestRec      = squeeze (XTest_rec(:,:,:,idx));
figure;
subplot(1,2,1)
montage (ImgTest, 'BorderSize', [1,1], 'BackgroundColor', col);
title('Original');
subplot(1,2,2);
montage (ImgTestRec, 'BorderSize', [1,1], 'BackgroundColor', col);
title('Reconstructed');
sgtitle(sprintf('Test (RMS: %.4f)',errorXTest))
return

function X = preprocessMiniBatch(dataX)

% Concatenate.
X                       = cat(4,dataX{:});

end

function [loss,gradientsE,gradientsD] = modelLoss(netE,netD,X)

% Forward through encoder.
[Z,mu,logSigmaSq]       = forward(netE,X);

% Forward through decoder.
Y                       = forward(netD,Z);

% Calculate loss and gradients.
loss                    = elboLoss(Y,X,mu,logSigmaSq);
[gradientsE,gradientsD] = dlgradient(loss,netE.Learnables,netD.Learnables);

end

function loss = elboLoss(Y,T,mu,logSigmaSq)

% Reconstruction loss.
reconstructionLoss = mse(Y,T);

% KL divergence.
KL = -0.5 * sum(1 + logSigmaSq - mu.^2 - exp(logSigmaSq),1);
KL = mean(KL);

% Combined loss.
loss = reconstructionLoss + KL;

end

function Y = modelPredictions(netE,netD,mbq)

Y               = [];

% Loop over mini-batches.
while hasdata(mbq)
    X           = next(mbq);

    % Forward through encoder.
    Z           = predict(netE,X);

    % Forward through dencoder.
    XGenerated  = predict(netD,Z);

    % Extract and concatenate predictions.
    Y           = cat(4,Y,extractdata(XGenerated));
end

end