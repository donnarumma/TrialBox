function   TEST_MNIST_VAE(ifplot)
% function TEST_MNIST_VAE(ifplot)
%--------------------------------------------------------------------------
if ~exist('ifplot','var'); ifplot=true; end
%--------------------------------------------------------------------------
%% load MNIST in Matlab
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

%% learn network model: encoder (netE) and decoder (netD)
par.vaeModel                        = vaeModelParams();
par.vaeModel.miniBatchSize          = 128;
par.vaeModel.graphplot              = ifplot;
% par.vaeModel.numEpochs            = 10;
par.vaeModel.InField                = 'MNIST';
[data_trials_train, out.vaeModel]   = vaeModel(data_trials_train,par.vaeModel);

ref         = datetime('now','Format','yyyyMMddHHmmss');
fn          = sprintf('~/TESTS/MNIST/vae/%s_%s.mat',mfilename,ref);
save(fn,'out','par')

%% Encode manifold
par.vaeEncode                       = vaeEncodeParams();
par.vaeEncode.netE                  = out.vaeModel.netE;
par.vaeEncode.InField               = 'MNIST';
par.vaeEncode.OutField              = 'vae';
par.vaeEncode.miniBatchSize         = par.vaeModel.miniBatchSize;
par.vaeEncode.MBnumOutputs          = par.vaeModel.MBnumOutputs;
% train8
data_trials_train                   = vaeEncode(data_trials_train,par.vaeEncode);
% test
data_trials_test                    = vaeEncode(data_trials_test ,par.vaeEncode);
 
%% Decode manifold
par.vaeDecode                       = vaeDecodeParams();
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

%% reconstruction error
par.reconstructionErrors                        = reconstructionErrorsParams();
par.reconstructionErrors.FieldReconstructed     = 'MNIST_rec';
par.reconstructionErrors.FieldTarget            = 'MNIST';
out.errorsTrain                                 = reconstructionErrors(data_trials_train,par.reconstructionErrors);
out.errorsTest                                  = reconstructionErrors(data_trials_test, par.reconstructionErrors);
errorXTrain                                     = out.errorsTrain.RMSE;
errorXTest                                      = out.errorsTest.RMSE;
fprintf('Train error: %g\n',errorXTrain)
fprintf('Test  error: %g\n',errorXTest)

%% plot TRAIN
hfg.trainRec                = figure('visible',ifplot);
hmShow                      = 49;

par.plot_Montage            = plot_MontageParams;
par.plot_Montage.hfig       = hfg.trainRec;
par.plot_Montage.hmShow     = randperm(length(data_trials_train),hmShow);
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
hfg.testRec                 = figure('visible',ifplot);
par.plot_Montage            = plot_MontageParams;
par.plot_Montage.hfig       = hfg.testRec;   
par.plot_Montage.hmShow     = randperm(length(data_trials_test),hmShow);
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