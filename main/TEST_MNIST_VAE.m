function   out = TEST_MNIST_VAE(pms)
% function out = TEST_MNIST_VAE(pms)
%--------------------------------------------------------------------------
if nargin < 1
    pms.ifplot      = true;
    pms.root_dir    = '~/TESTS/MNIST/vae/';
    pms.model_file  = false;
    pms.numEpochs   = 150;
end
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

%% learn autoencoder network model: encoder (netE) and decoder (netD)
if model_file
    lmodel=load(model_file);
    netE    = lmodel.out.vaeModel.netE;
    netD    = lmodel.out.vaeModel.netD;
else
    netE    = [];
    netD    = [];
end

par.vaeModel                        = vaeModelParams();
par.vaeModel.miniBatchSize          = 128;
par.vaeModel.numLatentChannels      = 32; % dimension of the manifold
par.vaeModel.graphplot              = ifplot;
par.vaeModel.numEpochs              = numEpochs;
par.vaeModel.InField                = 'MNIST';
par.vaeModel.netE                   = netE;
par.vaeModel.netD                   = netD;
[data_trials_train, out.vaeModel]   = vaeModel(data_trials_train,par.vaeModel);

%% Encode manifold in numLatentChannels dimension
par.vaeEncode                       = vaeEncodeParams();
par.vaeEncode.netE                  = out.vaeModel.netE;
par.vaeEncode.InField               = 'MNIST';
par.vaeEncode.OutField              = 'vae';
% train
data_trials_train                   = vaeEncode(data_trials_train,par.vaeEncode);
% test
data_trials_test                    = vaeEncode(data_trials_test ,par.vaeEncode);
 
%% Decode manifold into the original space
par.vaeDecode                       = vaeDecodeParams();
par.vaeDecode.netD                  = out.vaeModel.netD;
par.vaeDecode.exec                  = true;
par.vaeDecode.InField               = 'vae';
par.vaeDecode.OutField              = 'MNIST_rec';
par.vaeDecode.xfld                  = 'time';
% train
data_trials_train                   = vaeDecode(data_trials_train,par.vaeDecode);
% test
data_trials_test                    = vaeDecode(data_trials_test ,par.vaeDecode);

%% reconstruction error (train and test)
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
% hfig for train
hfg.trainRec                = figure('visible',ifplot);
hmShow                      = 49;

par.plot_Montage            = plot_MontageParams;
par.plot_Montage.hfig       = hfg.trainRec;
par.plot_Montage.hmShow     = randperm(length(data_trials_train),hmShow);
tcl=tiledlayout(1,2, 'Padding', 'compact', 'TileSpacing', 'tight');
% plot original
nexttile(tcl)
par.plot_Montage.InField    ='MNIST';
plot_Montage(data_trials_train,par.plot_Montage); 
title('Original')
sgtitle(sprintf('TRAIN: RMS=%.2f',errorXTrain))
% plot reconstructed
nexttile(tcl)
par.plot_Montage.InField    ='MNIST_rec';
plot_Montage(data_trials_train,par.plot_Montage); 
title('Reconstructed')
par.plot_Montage.hfig       = [];
%% plot TEST
hfg.testRec                 = figure('visible',ifplot);
par.plot_Montage            = plot_MontageParams;
par.plot_Montage.hfig       = hfg.testRec;   
par.plot_Montage.hmShow     = randperm(length(data_trials_test),hmShow);
tcl=tiledlayout(1,2, 'Padding', 'compact', 'TileSpacing', 'tight');
% plot original
nexttile(tcl)
par.plot_Montage.InField    ='MNIST';
plot_Montage(data_trials_test,par.plot_Montage); 
title('Original')
% plot Reconstructed
nexttile(tcl)
par.plot_Montage.InField    ='MNIST_rec';
plot_Montage(data_trials_test,par.plot_Montage); 
title('Reconstructed')
sgtitle(sprintf('TEST: RMS=%.2f',errorXTest))
par.plot_Montage.hfig       = [];
%% plot hfg
if ~ifplot
    par.hfigPrint               = hfigPrintParams();
    % save pdf with time ref
    par.hfigPrint.pdf_file      = sprintf('%s%s_%s.pdf',save_dir,description,ref);%[save_dir filesep  description '_' ref '.pdf'];
    par.hfigPrint.save_dir      = save_dir; 
    hfigPrint(hfg,par.hfigPrint)
end
%% save learned model
model_file                      = sprintf('%s%s_%s.mat',save_dir,description,ref);
if ~isfolder(save_dir); mkdir(save_dir); end
save(model_file,'out','par')
out.model_file                  = model_file;

return