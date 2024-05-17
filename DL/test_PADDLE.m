%--------------------------------------------------------------------------
% Basso, Curzio, et al. "Paddle: proximal algorithm for dual dictionaries learning." Artificial Neural Networks and Machine Learning–ICANN 2011: 21st International Conference on Artificial Neural Networks, Espoo, Finland, June 14-17, 2011, Proceedings, Part I 21. Springer Berlin Heidelberg, 2011.
% Basso, C., Santoro, M., Verri, A., & Villa, S. (2011). 
% Paddle: proximal algorithm for dual dictionaries learning. 
% In Artificial Neural Networks and Machine Learning–ICANN 2011: 
% 21st International Conference on Artificial Neural Networks, 
% Espoo, Finland, June 14-17, 2011, Proceedings, Part I 21 (pp. 379-386). 
% Springer Berlin Heidelberg.
% https://link.springer.com/chapter/10.1007/978-3-642-21735-7_47
%--------------------------------------------------------------------------
% Test settings paddle
%--------------------------------------------------------------------------
nAtoms                      = 100;
eta                         = 1;
tau                         = 0.02;
numIter                     = 500;
rtoll                       = 1E-05;
%--------------------------------------------------------------------------

%% load MNIST in Matlab
[XTrain,YTrain,anglesTrain] = digitTrain4DArrayData;
[XTest, YTest, anglesTest ] = digitTest4DArrayData;
X3dTrain                    = squeeze(XTrain);
X3dTest                     = squeeze(XTest);

%% arrange
[nChannels,nTimes,nTrials]  = size(X3dTrain);
% nChannels*nTimes x nTrials 
XdataTrain                  = reshape(X3dTrain,nChannels*nTimes,nTrials);
XdataTest                   = reshape(X3dTest,nChannels*nTimes,nTrials);
% nTrials x nChannels*nTimes
XdataTrain                  = XdataTrain';
XdataTest                   = XdataTest';
%%
% learnPADDLE_FISTA d x n samples   -> d=nChannels*nTimes, n=nTrials 
% [U, D, C]   = learnPADDLE_FISTA(XdataTrain,  nAtoms, eta, tau, rtoll);
[Manifold, Dpaddle, Wpaddle]= learnPADDLE_FISTA(XdataTrain',  nAtoms, eta, tau, rtoll);
% Dpaddle -> nChannels*nTimes x nAtoms
% Wpaddle -> nAtoms x nChannels*nTimes 

%% manifold construct by simPADDLE     -> see below for an alternative
% nAtoms x nTrials
ZdataTrain                  = simPADDLE(XdataTrain', Dpaddle, tau, rtoll); 
ZdataTest                   = simPADDLE(XdataTest', Dpaddle, tau, rtoll);
% nTrials x nAtoms
ZdataTrain                  = ZdataTrain';  
ZdataTest                   = ZdataTest';
%% manofold construct by Coding Matrix -> see above for an alternative
% nTrials x nAtoms
ZdataTrain                  = XdataTrain*Wpaddle'; % XdataTrainrec   = ((XdataTrain*C')*D');
% nTrials x nAtoms 
ZdataTest                   = XdataTest*Wpaddle';  % XdataTestrec    = ((XdataTest*C')*D');
%% reconstruction
% nTrials x nChannels*nTimes
XdataTrainrec               = ZdataTrain*Dpaddle';
XdataTestrec                = ZdataTest*Dpaddle';
% recontruction error
errorXTrain                 = RMSError(XdataTrainrec,XdataTrain);
errorXTest                  = RMSError(XdataTestrec,XdataTest);
out.errorXTrain             = errorXTrain;
out.errorXTest              = errorXTest;
disp(out)
%% Learn classifier
t=tic;
fprintf('Learning Classifier...')
par.kfold           = 4;
par.numIterations   = 300;
cvp                 = cvpartition(YTrain,'kfold',par.kfold,'Stratify',true);

% fitdiscr input -> nTrials x nAtoms
mdl = fitcdiscr(ZdataTrain, YTrain,...
    'OptimizeHyperparameters','all',...
    'HyperparameterOptimizationOptions', ...
    struct('Optimizer','randomsearch','CVPartition',cvp,'MaxObjectiveEvaluations',par.numIterations, 'AcquisitionFunctionName','expected-improvement-plus','ShowPlots',false,'Verbose',0));
fprintf('Elapsed time %g s\n',toc(t));
%% accuracy on Test
[YTest_pred,probs]          = predict(mdl, ZdataTest); % nTrials x nClasses
out.accuracy = mean(YTest_pred==YTest);
out.accuracyclass = accuracy4classes(YTest_pred,YTest);
disp(out)
%% plot a random number hmShow of test vs reconstructions
numObservations = size(XTest,4);
hmShow          = 49;
idx             = randperm(numObservations,hmShow);
ImTest          = imtile(reshape (XdataTest(idx,:)',nChannels,nTimes,hmShow) ) ;
ImTestRec       = imtile(reshape (XdataTestrec(idx,:)',nChannels,nTimes,hmShow) ) ;
figure
subplot(1,2,1);imshow(ImTest);title('Test Original');
subplot(1,2,2);imshow(ImTestRec);title(sprintf('Reconstructed | Rms Error: %.4f',errorXTest));
%% plot a number hmShow of train vs reconstructions
col             = [0.5,0.5,0.7]; % col borders
hmShow          = 49;
numObservations = size(XTrain,4);
idx             = randperm(numObservations,hmShow);
ImgTrain        = reshape (XdataTrain(idx,:)',nChannels,nTimes,hmShow);
ImgTrainRec     = reshape (XdataTrainrec(idx,:)',nChannels,nTimes,hmShow);
subplot(1,2,1)
montage (ImgTrain, 'BorderSize', [1,1], 'BackgroundColor', col);
title('Original');
subplot(1,2,2);
montage (ImgTrainRec, 'BorderSize', [1,1], 'BackgroundColor', col);
title('Reconstructed');
sgtitle(sprintf('Train (RMS: %.4f)',errorXTrain))
% test
numObservations = size(XTest,4);
idx             = randperm(numObservations,hmShow);
ImgTest         = reshape (XdataTest(idx,:)',nChannels,nTimes,hmShow);
ImgTestRec      = reshape (XdataTestrec(idx,:)',nChannels,nTimes,hmShow);
figure;
subplot(1,2,1)
montage (ImgTest, 'BorderSize', [1,1], 'BackgroundColor', col);
title('Original');
subplot(1,2,2);
montage (ImgTestRec, 'BorderSize', [1,1], 'BackgroundColor', col);
title('Reconstructed');
sgtitle(sprintf('Test (RMS: %.4f)',errorXTest))
%%
% ImTrainRec      = imtile(reshape (XdataTrainrec(idx,:)',nChannels,nTimes,hmShow) ) ;
figure;
subplot(1,2,1);imshow(ImTrain);   title('Train Original');
subplot(1,2,2);imshow(ImTrainRec);title(sprintf('Reconstructed | Rms Error: %.4f',errorXTrain));
%% plot test vs reconstruction (imtile)
numObservations = size(XTest,4);
figure;
idx             = randperm(numObservations,hmShow);

ImTest          = imtile(reshape (XdataTest(idx,:)',nChannels,nTimes,hmShow) ) ;
ImTestRec       = imtile(reshape (XdataTestrec(idx,:)',nChannels,nTimes,hmShow) ) ;
figure
subplot(1,2,1);imshow(ImTest);title('Test Original');
subplot(1,2,2);imshow(ImTestRec);title(sprintf('Reconstructed | Rms Error: %.4f',errorXTest));
%% plot dictionaries (use montage) 
% encoding
figure;
tiledlayout(1,2, 'Padding', 'tight', 'TileSpacing', 'tight');
nexttile;
hold on;
title ('Encoding Dictionary')
IW = Wpaddle';
% IW = rescale(IW,0,1); % uniformly scale
for iAtom = 1:nAtoms
    IW(:,iAtom) = rescale(IW(:,iAtom),0,1);
end
IW=reshape (IW,nChannels,nTimes,nAtoms);
montage (IW, 'BorderSize', [1,1], 'BackgroundColor', col);
% decoding
nexttile;
hold on;
title ('Decoding Dictionary')
ID = Dpaddle;
% ID = rescale(ID,0,1); % uniformly scale
for iAtom = 1:nAtoms
    ID(:,iAtom) = rescale(ID(:,iAtom),0,1);
end
ID=reshape (ID,nChannels,nTimes,nAtoms);
montage (ID, 'BorderSize', [1,1], 'BackgroundColor', col);
return