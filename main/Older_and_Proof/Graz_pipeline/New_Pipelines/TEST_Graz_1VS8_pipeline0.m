% TEST_Graz_1vs8_pipeline0.m

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

signal_name                     = 'eeg';
signal_process                  = 'CSP';
EEG_sub = struct();
for indsub=1:9
    %% Extract and Arrange TRAIN Data
    par.extractGraz.signal_name                  = signal_name;
    par.extractGraz.InField                      = 'train';
    [EEG_sub(indsub).train,fsample] = extractGraz(indsub,par.extractGraz);
end
for indsub=1:9
    %% Extract and Arrange Test Data
    par.extractGraz.signal_name                  = signal_name;
    par.extractGraz.InField                      = 'test';
    [EEG_sub(indsub).test,fsample] = extractGraz(indsub,par.extractGraz);
end

cvall = cvpartition(size(EEG_sub,2),'LeaveOut');
for i=1:size(EEG_sub,2)
    indices = training(cvall,i);
    test = (indices == 0);
    train = ~test;

    EEG_app = EEG_sub(train);
    EEG_train = vertcat(EEG_sub.train);

    for iTr=1:size(EEG_train,1)
        EEG_train(iTr).trialId = iTr;
    end

    EEG_test = EEG_sub(test).test;
    subj_test = find(test==1);

    StartClass = unique([EEG_train.trialType]);
    % Time Interpolation and selection Trials [0.75;2.5] from CUE (Motor Imagery Interval)
    par.TimeSelect               = TimeSelectParams;
    par.TimeSelect.t1            = itsub1(indsub,1) - 0.75; % in s from ZeroEvent time
    par.TimeSelect.t2            = itsub1(indsub,1) + 0.75; % in s from ZeroEvent time
    par.TimeSelect.InField       = signal_name;
    par.TimeSelect.OutField      = signal_name;
    par.TimeSelect.dt            = 1;

    itr1 = par.TimeSelect.t1;
    itr2 = par.TimeSelect.t2;

    % Filter Bank
    par.FilterBankCompute            = FilterBankComputeParams();
    par.FilterBankCompute.InField    = signal_name;
    par.FilterBankCompute.OutField   = signal_name;
    par.FilterBankCompute.attenuation = 10;
    par.FilterBankCompute.FilterBank = 'Nine';
    par.FilterBankCompute.fsample    = fsample;

    par.exec.funname ={'TimeSelect','FilterBankCompute'};
    EEG_train =run_trials(EEG_train,par);

    % Time Interpolation and selection Trials [0.75;2.5] from CUE (Motor Imagery Interval)
    par.TimeSelect               = TimeSelectParams;
    par.TimeSelect.t1            = itsub2(indsub,1) - 0.75; % in s from ZeroEvent time
    par.TimeSelect.t2            = itsub2(indsub,1) + 0.75; % in s from ZeroEvent time
    par.TimeSelect.InField       = signal_name;
    par.TimeSelect.OutField      = signal_name;
    par.TimeSelect.dt            = 1;

    its1 = par.TimeSelect.t1;
    its2 = par.TimeSelect.t2;

    % Filter Bank
    par.FilterBankCompute            = FilterBankComputeParams();
    par.FilterBankCompute.InField    = signal_name;
    par.FilterBankCompute.OutField   = signal_name;
    par.FilterBankCompute.attenuation = 10;
    par.FilterBankCompute.FilterBank = 'Nine';
    par.FilterBankCompute.fsample    = fsample;

    par.exec.funname ={'TimeSelect','FilterBankCompute'};
    EEG_test =run_trials(EEG_test,par);

    %% Step 2. perform CSP
    % CSP Dictionary evaluation on train
    par.cspModel                  = cspModelParams;
    par.cspModel.m                = 2;
    par.cspModel.InField          = signal_name;
    par.cspModel.OutField         = signal_process;

    [~,out.cspModel] = cspModel(EEG_train,par.cspModel);

    % CSP Encode on train and test data
    par.cspEncode                  = cspEncodeParams;
    par.cspEncode.InField          = signal_name;
    par.cspEncode.OutField         = signal_process;
    par.cspEncode.W                = out.cspModel.W;

    par.exec.funname ={'cspEncode'};
    EEG_train =run_trials(EEG_train,par);
    EEG_test =run_trials(EEG_test,par);

    TotalFeatures = size(EEG_test(1).(signal_process),2);

    % Mutual Information
    par.miModel               = miModelParams;
    par.miModel.InField       = signal_process;
    par.miModel.m             = par.cspModel.m;

    [~, out.miModel]=miModel(EEG_train,par.miModel);

    par.miEncode               = miEncodeParams;
    par.miEncode.InField       = signal_process;
    par.miEncode.OutField      = signal_process;
    par.miEncode.IndMI         = out.miModel.IndMI;

    par.exec.funname ={'miEncode'};
    [EEG_train, out]=run_trials(EEG_train,par);
    V_train = out.miEncode.V;
    [EEG_test, out]=run_trials(EEG_test,par);
    V_test = out.miEncode.V;

    %% Step 3. Model Classification on CSP
    % qdaModel
    par.qdaModel                      = qdaModelParams;
    par.qdaModel.InField              = 'CSP';
    par.qdaModel.numIterations        = 100;
    par.qdaModel.kfold                = 5;
    [~, outQDA]               = qdaModel(EEG_train,par.qdaModel);

    % predictQDA
    par.mdlPredict                  = mdlPredictParams;
    par.mdlPredict.InField          = 'CSP';
    par.mdlPredict.OutField         = 'QDApred';
    par.mdlPredict.ProbField        = 'QDAProb';
    par.mdlPredict.mdl              = outQDA.mdl;
    [EEG_train,resQDA.train]        = mdlPredict(EEG_train,par.mdlPredict);
    [EEG_test,resQDA.test]          = mdlPredict(EEG_test,par.mdlPredict);

    % knnModel
    par.knnModel                      = knnModelParams;
    par.knnModel.InField              = 'CSP';
    par.knnModel.numIterations        = 100;
    par.knnModel.kfold                = 5;
    [~, outKNN]                       = knnModel(EEG_train,par.knnModel);

    % predictKNN
    par.mdlPredict                  = mdlPredictParams;
    par.mdlPredict.InField          = 'CSP';
    par.mdlPredict.OutField         = 'KNNpred';
    par.mdlPredict.ProbField        = 'KNNProb';
    par.mdlPredict.mdl              = outKNN.mdl;
    [EEG_train,resKNN.train]        = mdlPredict(EEG_train,par.mdlPredict);
    [EEG_test,resKNN.test]          = mdlPredict(EEG_test,par.mdlPredict);

    %% Naive Bayes Parzen Window
    % fitNBPW
    par.fitNBPW                     = fitNBPWParams;
    par.fitNBPW.labs_train          = [EEG_train.trialType]';
    par.fitNBPW.labs_test           = [EEG_train.trialType]';
    [pred_train,outNBPW.train]      = fitNBPW(V_train,V_train,par.fitNBPW);

    par.fitNBPW                 = fitNBPWParams;
    par.fitNBPW.labs_train      = [EEG_train.trialType]';
    par.fitNBPW.labs_test       = [EEG_test.trialType]';
    [pred_test,outNBPW.test]    = fitNBPW(V_train,V_test,par.fitNBPW);

    % predictNBPW Train
    par.predictNBPW                  = predictNBPWParams;
    par.predictNBPW.InField          = 'CSP';
    par.predictNBPW.OutField         = 'NBPWpred';
    par.predictNBPW.labs_pred        = pred_train;

    [EEG_train,resNBPW.train]     = predictNBPW(EEG_train,par.predictNBPW);
    % KappaValue on Train NBPW
    Cmatrxix_train = confusionmat([EEG_train.trialType]', pred_train);
    resNBPW.train.kappaValue  = kappaModel(Cmatrxix_train);

    % predictNBPW Test
    par.predictNBPW                  = predictNBPWParams;
    par.predictNBPW.InField          = 'CSP';
    par.predictNBPW.OutField         = 'NBPWpred';
    par.predictNBPW.labs_pred        = pred_test;

    [EEG_test,resNBPW.test]       = predictNBPW(EEG_test,par.predictNBPW);
    % KappaValue on Test NBPW
    Cmatrxix_test = confusionmat([EEG_test.trialType]', pred_test);
    resNBPW.test.kappaValue  = kappaModel(Cmatrxix_test);

    % Save Result
    % create Tab Result
    params.createStructResult.subj       = '1VS8';
    params.createStructResult.method     = 'CSP';
    params.createStructResult.file       = 'Graz';
    params.createStructResult.train_name = strcat('AT_ALL');
    params.createStructResult.train_tr1  = itr1;
    params.createStructResult.train_tr2  = itr2;
    params.createStructResult.test_name  = strcat('AE','_',string(subj_test));
    params.createStructResult.test_ts1   = its1;
    params.createStructResult.test_ts2   = its2;
    params.createStructResult.m             = par.cspModel.m;
    params.createStructResult.class         = findclass(EEG_train,StartClass);
    params.createStructResult.irng          = par.irng;
    params.createStructResult.Filter        = par.FilterBankCompute.FilterBank;
    params.createStructResult.n_Features    = size(EEG_train(1).(signal_process),2);
    params.createStructResult.indMi         = par.miModel.k;
    params.createStructResult.attenuation   = par.FilterBankCompute.attenuation;
    params.createStructResult.TotalFeatures = TotalFeatures;
    params.createStructResult.kfold         = 0;

    % QDA Accuracy save
    [ResultQDA_kappa,ResultQDA_Acc,ResultQDA_class_Acc] = createStructResult(resQDA,params.createStructResult);
    % Update Tab Result
    params.updateTab.dir                = 'D:\TrialBox_Results_excel\Graz_dataset';
    params.updateTab.name               = 'Graz_1vs8_pipeline0';
    params.updateTab.sheetnames         = 'QDA';
    updated_Result_tableAccQDA          = updateTab(ResultQDA_Acc,params.updateTab);

    params.updateTab.sheetnames         = 'KappaQDA';
    updated_Result_tableKappaQDA        = updateTab(ResultQDA_kappa,params.updateTab);

    params.updateTab.name               = 'Graz_1vs8_pipeline0_class';
    params.updateTab.sheetnames         = 'QDA';
    updated_Resultclass_tableAccQDA     = updateTab(ResultQDA_class_Acc,params.updateTab);

    % KNN Accuracy save
    [ResultKNN_kappa,ResultKNN_Acc,ResultKNN_class_Acc] = createStructResult(resKNN,params.createStructResult);
    % Update Tab Result
    params.updateTab.dir                = 'D:\TrialBox_Results_excel\Graz_dataset';
    params.updateTab.name               = 'Graz_1vs8_pipeline0';
    params.updateTab.sheetnames         = 'KNN';
    updated_Result_tableAccKNN          = updateTab(ResultKNN_Acc,params.updateTab);

    params.updateTab.sheetnames         = 'KappaKNN';
    updated_Result_tableKappaKNN        = updateTab(ResultKNN_kappa,params.updateTab);

    params.updateTab.name               = 'Graz_1vs8_pipeline0_class';
    params.updateTab.sheetnames         = 'KNN';
    updated_Resultclass_tableAccKNN     = updateTab(ResultKNN_class_Acc,params.updateTab);

    % NBPW Accuracy save
    [ResultNBPW_kappa,ResultNBPW_Acc,ResultNBPW_class_Acc] = createStructResult(resNBPW,params.createStructResult);
    % Update Tab Result
    params.updateTab.dir                = 'D:\TrialBox_Results_excel\Graz_dataset';
    params.updateTab.name               = 'Graz_1vs8_pipeline0';
    params.updateTab.sheetnames         = 'NBPW';
    updated_Result_tableAccNBPW         = updateTab(ResultNBPW_Acc,params.updateTab);

    params.updateTab.sheetnames         = 'KappaNBPW';
    updated_Result_tableKappaNBPW       = updateTab(ResultNBPW_kappa,params.updateTab);

    params.updateTab.name               = 'Graz_1vs8_pipeline0_class';
    params.updateTab.sheetnames         = 'NBPW';
    updated_Resultclass_tableAccNBPW    = updateTab(ResultNBPW_class_Acc,params.updateTab);
end