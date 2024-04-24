%% function TEST_GRAZ_CSP_QDA
clear;
ifplot                          = false;
par.data_dir                    = 'D:\D_Ausilio EEG\EEG_FITTS\Data_Graz_Extracted';
par.data_dir                    = '~/DATA/GRAZ/';
par.irng                        = 10;                  % for reproducibility
rng(par.irng);
%% Step 0. Load raw data and arrange trials
par.subject                     = 8;
fsample                         = 250; % Hertz
EEG_trials                      = arrangeGrazTrials(par);
%% Step 1. perform Filterbank-CSP and tsne
signal_name                     = 'eeg';
signal_process                  = 'CSP';
% Filter Bank
par.FilterBankCompute           = FilterBankComputeParams();
par.FilterBankCompute.InField   = signal_name;
par.FilterBankCompute.OutField  = signal_process;
par.FilterBankCompute.FilterBank= 'Nine';
par.FilterBankCompute.fsample   = fsample; 
% cspModel
par.cspModel                    = cspModelParams;
par.cspModel.m                  = 2;
par.cspModel.InField            = signal_process;
par.cspModel.OutField           = signal_process;
% dataSplit
par.dataSplit                   = dataSplitParams;
par.dataSplit.TrainPercentage   = 70;
% execute functions
par.exec.funname                = {'FilterBankCompute','dataSplit'} ;
EEG_trials                      = run_trials(EEG_trials,par);
% cspModel
[~,out.cspModel]                = cspModel(EEG_trials([EEG_trials.train]),par.cspModel);
% cspEncode on train data
par.cspEncode                   = cspEncodeParams;
par.cspEncode.W                 = out.cspModel.W;
par.cspEncode.InField           = signal_process;
par.cspEncode.OutField          = signal_process;
% tsneModel
par.tsneModel                   = tsneModelParams;
par.tsneModel.NumDimensions     = 2;
par.tsneModel.InField           = signal_process;
par.tsneModel.OutField          = 'tsne';
%
par.exec.funname                = {'cspEncode','tsneModel'}; 
EEG_trials                      = run_trials(EEG_trials,par);

%% Step 2. QDA Classification on Filterbank-CSP
% qdaModel
par.qdaModel                    = qdaModelParams();
par.qdaModel.InField            = 'CSP';
[~, out.qdaModel]               = qdaModel(EEG_trials([EEG_trials.train]),par.qdaModel);
disp(out.qdaModel);
% predictQDA
par.mdlPredict                  = mdlPredictParams;
par.mdlPredict.InField          = 'CSP'; 
par.mdlPredict.OutField         = 'pred';
par.mdlPredict.mdl              = out.qdaModel.mdl;   
% just to check test and train
[~,res]                         = mdlPredict(EEG_trials([EEG_trials.train]),par.mdlPredict);
disp(res);
[~,res]                         = mdlPredict(EEG_trials([EEG_trials.test]),par.mdlPredict);
disp(res);
% 
[EEG_trials,res]                = mdlPredict(EEG_trials,par.mdlPredict);

%% plot_AccuracyBars
cmaps                           = linspecer(length(unique([EEG_trials.trialType])));
par.plot_AccuracyBars           = plot_AccuracyBarsParams;
par.plot_AccuracyBars.cmaps     = cmaps;

par.plot_AccuracyBars.hfig      = figure('visible',ifplot);
hfg.AccuracyBarsTrain           = plot_AccuracyBars(EEG_trials([EEG_trials.train]),par.plot_AccuracyBars);
sgtitle('Train')
par.plot_AccuracyBars.hfig      = figure('visible',ifplot);
hfg.AccuracyBarsTest            = plot_AccuracyBars(EEG_trials([EEG_trials.test]),par.plot_AccuracyBars);
sgtitle('Test')


% %% plot Accuracy on Filterbank-CSP
% % plot_EachDimBar
% par.plot_EachDimBar             = plot_EachDimBarParams;
% par.plot_EachDimBar.hfig        = figure('visible',ifplot);
% par.plot_EachDimBar.novariance  = false;
% par.plot_EachDimBar.addbar      = false;
% par.plot_EachDimBar.cmaps       = cmaps;
% par.plot_EachDimBar.legplot     = 1;
% par.plot_EachDimBar.InField     = 'accuracy';
% par.plot_EachDimBar.novariance  = true;
% par.plot_EachDimBar.keep        = 1:4;
% % par.plot_EachDimBar.nCols     = 1;
% par.plot_EachDimBar.ylabel      = '$acc$';
% par.plot_EachDimBar.YLIM        = [0-0.01,1+0.01];
% par.plot_EachDimBar.evaltime    = 1;
% par.plot_EachDimBar.chanceline  = true;
% 
% data_trials                         = EEG_trials([EEG_trials.test]);
% % meanData -> get mean con classes
% par.meanData                        = meanDataParams;
% par.meanData.trialTypes             = [data_trials.trialType];
% par.meanData.opt                    = [0,0,1]; % mean
% par.meanData.InField                = 'success';
% par.meanData.OutField               = 'accuracy';
% class_acc_trials                    = meanData(data_trials,par.meanData);
% % meanData -> get global mean
% par.meanData                        = meanDataParams;
% par.meanData.opt                    = [0,0,1]; % mean
% par.meanData.trialTypes             = true(size([data_trials.trialType]));
% par.meanData.InField                = 'success';
% par.meanData.OutField               = 'accuracy';
% acc_trials                          = meanData(data_trials,par.meanData);   
% par.plot_EachDimBar.addbar          = acc_trials;
% par.plot_EachDimBar.titles{1}       = 'ACCURACY';
% hfg.plot_EachDimBar                 = plot_EachDimBar(class_acc_trials,par.plot_EachDimBar);

%% plot CSP_tsne
data_trials                         = EEG_trials([EEG_trials.train]);
par.plot_scatter                    = plot_scatterParams;
par.plot_scatter.InField            = 'tsne';
par.plot_scatter.cmaps              = linspecer(length(unique([data_trials.trialType])));
par.plot_scatter.hfig               = figure('visible',ifplot);
hfg.tnseScatterTrain                = plot_scatter(EEG_trials([EEG_trials.train]),par.plot_scatter);
title('tsne train')
par.plot_scatter.hfig               = figure('visible',ifplot);
hfg.tnseScatterTest                 = plot_scatter(EEG_trials([EEG_trials.test]),par.plot_scatter);
title('tsne test')

%% pSeparability
par.pSeparability                   = pSeparabilityParams;
par.pSeparability.InField           = par.tsneModel.OutField;
par.pSeparability.OutField          = 'comparisons';
par.pSeparability.exec              = true;
[pVals,pClasses]                    = pSeparability(EEG_trials([EEG_trials.test]),par.pSeparability);

%% pvalue plot per feature
par.plot_pValues                    = plot_pValuesParams;
par.plot_pValues.InField            = par.pSeparability.OutField;
par.plot_pValues.hfig               = figure('visible',ifplot);
par.plot_pValues.nRows              = 1;
par.plot_pValues.decisions          = [];
hfg.pvals                           = plot_pValues(pClasses,par.plot_pValues);
% sgtitle(hfg.pvals,[RatName ' ' daystr ' ' feedstr]);


%% save hfg
if ~ifplot
    save_dir                    = '/TESTS/GRAZ/NSA/';
    par.hfigPrint               = hfigPrintParams();
    par.hfigPrint.pdf_file      = [save_dir mfilename '.pdf'];
    par.hfigPrint.save_dir      = save_dir; 
    hfigPrint(hfg,par.hfigPrint)
end
