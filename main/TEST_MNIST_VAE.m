% function TEST_MNIST_VAE
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%% load MNIST in Matlab
[XTrain,YTrain,anglesTrain] = digitTrain4DArrayData;
[XTest, YTest, anglesTest ] = digitTest4DArrayData;

%% arrange data_trials_train and data_trials_test 
% nChannels x nTimes x nTrials 
X3dTrain                    = squeeze(XTrain);
X3dTest                     = squeeze(XTest);
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

%% learn network model: encoder (netE) and decoder (netD)
par.vaeModel                        = vaeModelParams();
par.vaeModel.miniBatchSize          = 128;
% par.vaeModel.numEpochs            = 10;
par.vaeModel.InField                = 'MNIST';
[data_trials_train, out.vaeModel]   = vaeModel(data_trials_train,par.vaeModel);

%% Encode manifold
par.vaeEncode                       = vaeEncodeParams();
par.vaeEncode.netE                  = out.vaeModel.netE;
par.vaeEncode.InField               = 'MNIST';
par.vaeEncode.OutField              = 'vae';
par.vaeEncode.miniBatchSize         = par.vaeModel.miniBatchSize;
par.vaeEncode.MBnumOutputs          = par.vaeModel.MBnumOutputs;
% train
data_trials_train                   = vaeEncode(data_trials_train,par.vaeEncode);
% test
data_trials_test                    = vaeEncode(data_trials_test ,par.vaeEncode);
%% Decode manifold
par.vaeEncode                       = vaeDecodeParams();
par.vaeDecode.netD                  = out.vaeModel.netD;
par.vaeDecode.exec                  = true;
par.vaeDecode.InField               = 'vae';
par.vaeDecode.OutField              = 'MNIST_rec';
par.vaeDecode.xfld                  = 'time';
par.vaeDecode.miniBatchSize         = par.vaeModel.miniBatchSize;
par.vaeDecode.MBnumOutputs          = par.vaeModel.MBnumOutputs;
% train
data_trials_train                   = vaeDecode(data_trials_train,par.vaeDecode);
% test
data_trials_test                    = vaeDecode(data_trials_test ,par.vaeDecode);

%%
%% reconstruction error
% train
% nChannels x nTimes x nTrials
X3dTrain_rec                = cat(3,data_trials_train.MNIST_rec);
% nTrials x nChannels x nTimes
X3dTrain_rec                = permute(X3dTrain_rec,[3,1,2]);
% nTrials x nChannels*nTimes
XdataTrain                  = reshape(X3dTrain    ,size(X3dTrain,1)    ,size(X3dTrain,2)    *size(X3dTrain,3));
XdataTrain_rec              = reshape(X3dTrain_rec,size(X3dTrain_rec,1),size(X3dTrain_rec,2)*size(X3dTrain_rec,3));
errorTrains                 = errors(XdataTrain_rec,XdataTrain);
errorXTrain                 = errorTrains.RMSE;
% test
X3dTest_rec                 = cat(3,data_trials_test.MNIST_rec);
% nTrials x nChannels x nTimes
X3dTest_rec                 = permute(X3dTest_rec,[3,1,2]);
% nTrials x nChannels*nTimes
XdataTest                   = reshape(X3dTest    ,size(X3dTest,1)    ,size(X3dTest,2)    *size(X3dTest,3));
XdataTest_rec               = reshape(X3dTest_rec,size(X3dTest_rec,1),size(X3dTest_rec,2)*size(X3dTest_rec,3));
errorTests                  = errors(XdataTest_rec,XdataTest);
errorXTest                  = errorTests.RMSE;
% errorXTest                  = RMSError(XdataTest_rec,XdataTest);
out.errorXTrain             = errorXTrain;
out.errorXTest              = errorXTest;
disp(out)
%% plot TRAIN
isplot                      = true;
hfg.plot_Montage            = figure('visible',isplot);
par.plot_Montage            = plot_MontageParams;
par.plot_Montage.hfig       = hfg.plot_Montage;   
tcl=tiledlayout(1,2, 'Padding', 'compact', 'TileSpacing', 'tight');
nexttile
par.plot_Montage.InField    ='MNIST';
hfg.plot_Montage            = plot_Montage(data_trials_train,par.plot_Montage); 
title('Original')
nexttile(tcl)
par.plot_Montage.InField    ='MNIST_rec';
hfg.plot_Montage            = plot_Montage(data_trials_train,par.plot_Montage); 
title('Reconstructed')
sgtitle(sprintf('TRAIN: RMS=%.2f',errorXTrain))
%% plot TEST
par.plot_Montage            = plot_MontageParams;
par.plot_Montage.hfig       = figure;   
tcl=tiledlayout(1,2, 'Padding', 'compact', 'TileSpacing', 'tight');
nexttile
par.plot_Montage.InField    ='MNIST';
hfg.plot_Montage            = plot_Montage(data_trials_test,par.plot_Montage); 
title('Original')
nexttile
par.plot_Montage.InField    ='MNIST_rec';
hfg.plot_Montage            = plot_Montage(data_trials_test,par.plot_Montage); 
title('Reconstructed')
sgtitle(sprintf('TEST: RMS=%.2f',errorXTest))
 
%% plot hfg
if ~ifplot
    save_dir                    = '~/TESTS/MNIST/vae/';
    description                 = mfilename;
    par.hfigPrint               = hfigPrintParams();
    par.hfigPrint.pdf_file      = [save_dir filesep description '.pdf'];
    par.hfigPrint.save_dir      = save_dir; 
    hfigPrint(hfg,par.hfigPrint)
end
return