% TEST_Graz_OnTest_postDPCA_pipeline6.m

% Dataset eeegplanetion: "Leeb, R., Brunner, C., Müller-Putz, G., Schlögl, A., & Pfurtscheller, G. J. G. U. O. T. (2008). BCI Competition 2008–Graz data set B.
%       Graz University of Technology, Austria, 16, 1-6."
% Link for training data: www.bbci.de/competition/iv/ -> Download of data sets -> agree submit -> Data sets 2a: ‹4-class motor imagery>
% Link for test data: www.bbci.de/competition/iv/ -> News -> Results of the
% BCI Competition IV -> True Labels of Competition's Evaluation Sets -> Data sets 2a:

% Specifics:
% 1:    Left Hand
% 2:    Right Hand
% 3:    Foot
% 4:    Tongue

clear; close all;

par.irng = 10;
rng(par.irng);

subject = 1:9;
itsub1 = [1.4,0.9,1.2,1.0,1.4,1.2,1.2,1.5,0.6];
itsub2 = [1.6,1.1,1.4,1.2,1.6,1.4,1.4,1.9,0.9];
for indsub = subject
    %
    % indsub = 2;
    signal_name                     = 'eeg';
    signal_process                  = 'pca';

    %% Extract and Arrange TRAIN Data
    par.extractGraz.signal_name                  = signal_name;
    par.extractGraz.InField                      = 'train';
    [EEG_train,fsample] = extractGraz(indsub,par.extractGraz);

    StartClass = unique([EEG_train.trialType]);
    % Time Interpolation and selection Trials [0.5;2.5] from CUE (Motor Imagery Interval)
    par.TimeSelect               = TimeSelectParams;
    par.TimeSelect.t1            = 0.5; % in s from ZeroEvent time
    par.TimeSelect.t2            = 2.5; % in s from ZeroEvent time
    par.TimeSelect.InField       = signal_name;
    par.TimeSelect.OutField      = signal_name;
    par.TimeSelect.dt            = 1;

    par.exec.funname ={'TimeSelect'};
    EEG_train1 =run_trials(EEG_train,par);

    par.TimeSelect.t1            = itsub1(indsub); % in s from ZeroEvent time
    par.TimeSelect.t2            = itsub2(indsub); % in s from ZeroEvent time

    EEG_train2=run_trials(EEG_train,par);

    itr1 = par.TimeSelect.t1;
    itr2 = par.TimeSelect.t2;

    par.concatenateEEG.InField = signal_name;
    par.concatenateEEG.fsample = fsample;
    EEG_train = concatenateEEG(EEG_train1,EEG_train2,par.concatenateEEG);

    % Filter Bank
    par.FilterBankCompute            = FilterBankComputeParams();
    par.FilterBankCompute.InField    = signal_name;
    par.FilterBankCompute.OutField   = signal_name;
    par.FilterBankCompute.attenuation = 20;
    par.FilterBankCompute.FilterBank = 'Nine';
    par.FilterBankCompute.fsample    = fsample;

    par.exec.funname ={'FilterBankCompute'};
    EEG_train =run_trials(EEG_train,par);

    %% Extract and Arrange Test Data
    par.extractGraz.signal_name                  = signal_name;
    par.extractGraz.InField                      = 'test';
    [EEG_test,fsample] = extractGraz(indsub,par.extractGraz);

    % Time Interpolation and selection Trials [0.5;2.5] from CUE (Motor Imagery Interval)
    par.TimeSelect               = TimeSelectParams;
    par.TimeSelect.t1            = 0.5; % in s from ZeroEvent time
    par.TimeSelect.t2            = 2.5; % in s from ZeroEvent time
    par.TimeSelect.InField       = signal_name;
    par.TimeSelect.OutField      = signal_name;
    par.TimeSelect.dt            = 1;

    its1 = par.TimeSelect.t1;
    its2 = par.TimeSelect.t2;

    par.exec.funname ={'TimeSelect'};
    EEG_test1 =run_trials(EEG_test,par);

    par.TimeSelect.t1            = itr1; % in s from ZeroEvent time
    par.TimeSelect.t2            = itr2; % in s from ZeroEvent time

    EEG_test2 =run_trials(EEG_test,par);
    EEG_test = concatenateEEG(EEG_test1,EEG_test2,par.concatenateEEG);

    % Filter Bank
    par.FilterBankCompute            = FilterBankComputeParams();
    par.FilterBankCompute.InField    = signal_name;
    par.FilterBankCompute.OutField   = signal_name;
    par.FilterBankCompute.attenuation = 20;
    par.FilterBankCompute.FilterBank = 'Nine';
    par.FilterBankCompute.fsample    = fsample;

    par.exec.funname ={'FilterBankCompute'};
    EEG_test =run_trials(EEG_test,par);

    %% Step 2. perform pca
    % pca Dictionary evaluation on train
    par.pcaSynergy.InField          = signal_name;
    par.pcaSynergy.OutField         = signal_process;
    par.pcaSynergy.numComponents    = 0;
    par.pcaSynergy.perc             = 95;

    par.exec.funname ={'pcaSynergy'};
    [EEG_train,out] =run_trials(EEG_train,par);

    % pca Encode on test data
    par.pcaSynergyEncode.InField = signal_name;
    par.pcaSynergyEncode.OutField = signal_process;
    par.pcaSynergyEncode.Wpca  = out.pcaSynergy.W;
    par.exec.funname ={'pcaSynergyEncode'};
    EEG_test =run_trials(EEG_test,par);

    %% Step 3. Model Classification on pca
    % qdaModel
    par.qdaModel                      = qdaModelParams;
    par.qdaModel.InField              = 'pca';
    par.qdaModel.numIterations        = 100;
    par.qdaModel.kfold                = 5;
    [~, outQDA]               = qdaModel(EEG_train,par.qdaModel);

    % predictQDA
    par.mdlPredict                  = mdlPredictParams;
    par.mdlPredict.InField          = 'pca';
    par.mdlPredict.OutField         = 'QDApred';
    par.mdlPredict.ProbField        = 'QDAProb';
    par.mdlPredict.mdl              = outQDA.mdl;
    [EEG_train,resQDA.train]        = mdlPredict(EEG_train,par.mdlPredict);
    [EEG_test,resQDA.test]          = mdlPredict(EEG_test,par.mdlPredict);

    % knnModel
    par.knnModel                      = knnModelParams;
    par.knnModel.InField              = 'pca';
    par.knnModel.numIterations        = 100;
    par.knnModel.kfold                = 5;
    [~, outKNN]             = knnModel(EEG_train,par.knnModel);

    % predictKNN
    par.mdlPredict                  = mdlPredictParams;
    par.mdlPredict.InField          = 'pca';
    par.mdlPredict.OutField         = 'KNNpred';
    par.mdlPredict.ProbField        = 'KNNProb';
    par.mdlPredict.mdl              = outKNN.mdl;
    [EEG_train,resKNN.train]        = mdlPredict(EEG_train,par.mdlPredict);
    [EEG_test,resKNN.test]          = mdlPredict(EEG_test,par.mdlPredict);


    % Save Result
    % create Tab Result
    params.createStructResult.subj       = indsub;
    params.createStructResult.method     = 'pca';
    params.createStructResult.file       = 'Graz';
    params.createStructResult.train_name = 'AT';
    params.createStructResult.train_tr1  = itr1;
    params.createStructResult.train_tr2  = itr2;
    params.createStructResult.test_name  = 'AE';
    params.createStructResult.test_ts1   = its1;
    params.createStructResult.test_ts2   = its2;
    params.createStructResult.class      = findclass(EEG_train,StartClass);
    params.createStructResult.irng       = par.irng;
    params.createStructResult.method     = signal_process;
    params.createStructResult.Filter     = par.FilterBankCompute.FilterBank;

    % QDA Accuracy save
    [ResultQDA_kappa,ResultQDA_Acc,ResultQDA_class_Acc] = createStructResult(resQDA,params.createStructResult);
    % Update Tab Result
    params.updateTab.dir                = 'D:\TrialBox_Results_excel\Graz_dataset';
    params.updateTab.name               = 'Graz_OnTest_postDPCA_pipe6';
    params.updateTab.sheetnames         = 'QDA';
    updated_Result_tableAccQDA          = updateTab(ResultQDA_Acc,params.updateTab);

    params.updateTab.sheetnames         = 'KappaQDA';
    updated_Result_tableKappaQDA        = updateTab(ResultQDA_kappa,params.updateTab);

    params.updateTab.name               = 'Graz_OnTest_class_postDPCA_pipe6';
    params.updateTab.sheetnames         = 'QDA';
    updated_Resultclass_tableAccQDA     = updateTab(ResultQDA_class_Acc,params.updateTab);

    % KNN Accuracy save
    [ResultKNN_kappa,ResultKNN_Acc,ResultKNN_class_Acc] = createStructResult(resKNN,params.createStructResult);
    % Update Tab Result
    params.updateTab.dir                = 'D:\TrialBox_Results_excel\Graz_dataset';
    params.updateTab.name               = 'Graz_OnTest_postDPCA_pipe6';
    params.updateTab.sheetnames         = 'KNN';
    updated_Result_tableAccKNN          = updateTab(ResultKNN_Acc,params.updateTab);

    params.updateTab.sheetnames         = 'KappaKNN';
    updated_Result_tableKappaKNN        = updateTab(ResultKNN_kappa,params.updateTab);

    params.updateTab.name               = 'Graz_OnTest_class_postDPCA_pipe6';
    params.updateTab.sheetnames         = 'KNN';
    updated_Resultclass_tableAccKNN     = updateTab(ResultKNN_class_Acc,params.updateTab);
end