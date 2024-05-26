%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%% load MNIST in Matlab
[XTrain,YTrain,anglesTrain] = digitTrain4DArrayData;
[XTest, YTest, anglesTest ] = digitTest4DArrayData;
% nChannels x nTimes x nTrials 
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

%% Test settings srssd
%--------------------------------------------------------------------------
logLambdaA          =-7;%[-8,-6];
% lambdaAtom = 10^-5;
% lambdaCoeff = 0.5;
logLambdaA          =-5;%[-8,-6];

lambdaARange        = 10.^logLambdaA;
lambdaCRange        = 1.5;%0.9;%0.5; % 0.02,0.8
nAtoms              = 100;

par.srssdModel      = srssdModelParams();
par.srssdModel.r    = nAtoms; % dictionary size

%%
% 2D struct Jenhatton
mode.rectangle      = true; 
mode.pi4            = true;
par.srssdModel.spG  = get_groups( mode, [nChannels nTimes] );
% srssd
% selectedStruct      = 'ContinuousTemporal';
% structureType  = {'SingleDof','ContinuousTemporal','FixedTemporal','SingleDofCT'};
% selectedStruct          = 'SingleDofCT';
% par.srssdModel.spG      = setActionSPG(nTimes,nChannels,selectedStruct);

par.srssdModel.eta      = 1; % default 0
par.srssdModel.lambda   = 10.^logLambdaA; % lambdaARange -> Regolarization with respect to atoms
par.srssdModel.lambda2  = lambdaCRange;   % lambdaCRange -> Regolarization with respect to Coefficients
InField = 'MNIST'; 
par.srssdModel.InField  = InField;
%% 
time    = 1:nTimes;
xfld    = 'time';
labels  = unique(YTrain);
for iTrial=1:nTrials
    data_trials(iTrial).(InField)       = X3dTrain(:,:,iTrial);
    data_trials(iTrial).([xfld InField])= time;
    data_trials(iTrial).trialType       = find(ismember(labels,YTrain(iTrial)))-1;
    data_trials(iTrial).trialName       = num2str(data_trials(iTrial).trialType);
    data_trials(iTrial).behavior        = anglesTrain(iTrial);
    data_trials(iTrial).train           = true;
    data_trials(iTrial).test            = false;
end

for iTrial=1:size(X3dTest,3)
    data_trials_test(iTrial).(InField)  = X3dTest(:,:,iTrial);
    data_trials_test(iTrial).([xfld InField])= time;
    data_trials_test(iTrial).trialType  = find(ismember(labels,YTest(iTrial)))-1;
    data_trials_test(iTrial).trialName  = num2str(data_trials_test(iTrial).trialType);
    data_trials_test(iTrial).behavior   = anglesTest(iTrial);
    data_trials_test(iTrial).train      = false;
    data_trials_test(iTrial).test       = true;
end

%% model learning
% sr_ssd: n x d -> n=observations d=features -> n=nTrials, d=nChannels*nTimes
% [U, D, C]   = sr_ssd(X, spG, params);
% fprintf('Learning srssd Model...'); t=tic; 
% [Manifold, Dsrssd, Wsrssd]  = sr_ssd(XdataTrain, spG, par.srssdModel);
% fprintf('Elapsed time %.2f s\n',toc(t));

% Dsrssd -> nChannels*nTimes x nAtoms
% Wsrssd -> nAtoms x nChannels*nTimes  ????
%% srssdModel 
[~,out.srssdModel]      = srssdModel(data_trials,par.srssdModel);
%% manifold construction by Coding Matrix with srssdEncode -> see sim_sr_ssd to an alternative
par.srssdEncode             = srssdEncodeParams;
par.srssdEncode.InField     = 'MNIST';
par.srssdEncode.OutField    = 'srssd';
par.srssdEncode.Wsrssd      = out.srssdModel.Wsrssd;           % nAtoms x nChannels*nTimes
data_trials                 = srssdEncode(data_trials,par.srssdEncode);
data_trials_test            = srssdEncode(data_trials_test,par.srssdEncode);

% nAtoms x nTrials
ZdataTrain                  = cat(2,data_trials.srssd); 
ZdataTest                   = cat(2,data_trials_test.srssd); 
% nTrials x nAtoms
ZdataTrain                  = ZdataTrain';
ZdataTest                   = ZdataTest';
%% decode (reconstruction)
par.srssdDecode             = srssdDecodeParams();
par.srssdDecode.Dsrssd      = out.srssdModel.Dsrssd;
par.srssdDecode.InField     ='srssd';
par.srssdDecode.OutField    ='MNIST_rec';
par.srssdDecode.nChannels   = nChannels;
par.srssdDecode.nTimes      = nTimes;
data_trials                 = srssdDecode(data_trials,par.srssdDecode);
data_trials_test            = srssdDecode(data_trials_test,par.srssdDecode);
%% reconstruction error
% Train
% nChannels x nTimes x nTrials
X3dTrain_rec                = cat(3,data_trials.MNIST_rec);
% nTrials x nChannels x nTimes
X3dTrain_rec                = permute(X3dTrain_rec,[3,1,2]);
% nTrials x nChannels*nTimes
XdataTrain_rec              = reshape(X3dTrain_rec,size(X3dTrain_rec,1),size(X3dTrain_rec,2)*size(X3dTrain_rec,3));
errorTrains                 = errors(XdataTrain_rec,XdataTrain);
errorXTrain                 = errorTrains.RMSE;
% errorXTrain                 = RMSError(XdataTrain_rec,XdataTrain);
% test
X3dTest_rec                 = cat(3,data_trials_test.MNIST_rec);
% nTrials x nChannels x nTimes
X3dTest_rec                 = permute(X3dTest_rec,[3,1,2]);
% nTrials x nChannels*nTimes
XdataTest_rec               = reshape(X3dTest_rec,size(X3dTest_rec,1),size(X3dTest_rec,2)*size(X3dTest_rec,3));
errorTests                  = errors(XdataTest_rec,XdataTest);
errorXTest                  = errorTests.RMSE;
% errorXTest                  = RMSError(XdataTest_rec,XdataTest);
out.errorXTrain             = errorXTrain;
out.errorXTest              = errorXTest;
disp(out)
disp(errorXTrain)
disp(errorXTest)
%% Learn classifier
par.qdaModel                = qdaModelParams();
par.qdaModel.numIterations  = 300;
par.qdaModel.InField        = 'srssd';
[data_trials,out.qdaModel]  = qdaModel(data_trials,par.qdaModel);
%% predict on Test
par.mdlPredict              = mdlPredictParams;
par.mdlPredict.InField      = 'srssd'; 
par.mdlPredict.OutField     = 'pred';
par.mdlPredict.mdl          = out.qdaModel.mdl;   
%
[~,res]                     = mdlPredict(data_trials_test,par.mdlPredict);
disp(res);

%% plot TRAIN
isplot                      = true;
hfg.plot_Montage            = figure('visible',isplot);
par.plot_Montage            = plot_MontageParams;
par.plot_Montage.hfig       = hfg.plot_Montage;   
tcl=tiledlayout(1,2, 'Padding', 'compact', 'TileSpacing', 'tight');
nexttile
par.plot_Montage.InField    ='MNIST';
hfg.plot_Montage            = plot_Montage(data_trials,par.plot_Montage); 
title('Original')
nexttile(tcl)
par.plot_Montage.InField    ='MNIST_rec';
hfg.plot_Montage            = plot_Montage(data_trials,par.plot_Montage); 
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
%% plot Dictionaries 
% decoding
par.arrangeDictionary           = arrangeDictionaryParams();
par.arrangeDictionary.nChannels = nChannels;
par.arrangeDictionary.nTimes    = nTimes;
par.arrangeDictionary.dictionary= out.srssdModel.Dsrssd; % decoding
Wsrssd_trials                   = arrangeDictionary([],par.arrangeDictionary);
%
par.plot_Atoms                  = plot_AtomsParams();
par.plot_Atoms.nRows            = 10;
hfg.Decoding                    = plot_Atoms(Wsrssd_trials,par.plot_Atoms);
sgtitle('Dictionary (Decoding)')
% encoding
par.arrangeDictionary.dictionary= out.srssdModel.Wsrssd; % encoding
Dsrssd_trials                   = arrangeDictionary([],par.arrangeDictionary);
%
hfg.Encoding                    = plot_Atoms(Dsrssd_trials,par.plot_Atoms);
sgtitle('Dictionary (Encoding)')
%% plot Manifold
par.meanData                        = meanDataParams;
par.meanData.exec                   = [];
par.meanData.trialTypes             = [data_trials.trialType];
par.meanData.opt                    = [0,0,1]; % mean
par.meanData.InField                = 'srssd';
par.meanData.OutField               = 'usage';
class_usage_trials                  = meanData(data_trials,par.meanData);
% meanData -> get all trial means
% par.meanData                        = meanDataParams;
% par.meanData.exec                   = [];
% par.meanData.opt                    = [0,0,1]; % mean
% par.meanData.trialTypes             = true(size([data_trials_prob_plot.trialType]));
% par.meanData.InField                = SuccessField;
% par.meanData.OutField               = 'accuracy';
% acc_trials                          = meanData(data_trials_prob_plot,par.meanData);   

% plot_EachDimBar
par.plot_EachDimBar                 = plot_EachDimBarParams;
par.plot_EachDimBar.exec            = [];
par.plot_EachDimBar.hfig            = [];
% par.plot_EachDimBar.evaltime        = evaltime;
par.plot_EachDimBar.novariance      = false;
par.plot_EachDimBar.addbar          = false;
par.plot_EachDimBar.cmaps           = linspecer(length(unique([data_trials.trialType])));
par.plot_EachDimBar.legplot         = 1;
par.plot_EachDimBar.InField         = 'usage';
par.plot_EachDimBar.novariance      = false;
par.plot_EachDimBar.keep            = unique([class_usage_trials.trialType]);
par.plot_EachDimBar.nCols           = 1;
par.plot_EachDimBar.ylabel          = '$acc$';
% par.plot_EachDimBar.YLIM            = [0-0.01,1+0.01];
% par.plot_EachDimBar.chanceline      = true;
% GRAPH SUBPLOT TITLES - one graph for each subset
nSets        = 1;%length(channelSets); %% number of graphs
str          = cell(nSets,1);
% for iSet=1:nSets
%     str{iSet}='';
% end
% for iSet=1:nSets
%     channels     =channelSets{iSet};
%     nChannels    =length(channels);
%     % str{iSet}='';
%     for ichannel=1:nChannels
%         str{iSet} = sprintf('%s$${\\mathbf x}_{%d,t}$$',str{iSet},channels(ichannel));
%         if ichannel<nChannels
%              str{iSet} = sprintf('%s,',str{iSet});
%         end 
%     end
%     if ~isempty(explained)
%         str{iSet} = sprintf('%s (%2.1f)',str{iSet},sum(explained(channelSets{iSet})));
%     end
% end
par.plot_EachDimBar.titles{1}   = 'USAGE';
% par.plot_EachDimBar.addbar      = acc_trials;%false;%data_trials_prob_sep(1);
hfg.plot_EachDimBar             = plot_EachDimBar(class_usage_trials,par.plot_EachDimBar);
sgtitle('Usage')
%%