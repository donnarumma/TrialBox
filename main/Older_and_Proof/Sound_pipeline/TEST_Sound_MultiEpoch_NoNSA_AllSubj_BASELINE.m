% TEST_Baseline_NoNSA_AllSubj.m


% Specifics:
% label indexes
% 1: COHERENT ENGLISH
% 2: COHERENT ITALIAN
% 3: INCOHERENT ENGLISH
% 4: INCOHERENT ITALIAN
% 5: SCUMBLED

clear; close all;

par.irng = 10;
rng(par.irng);


signal_name                     = 'eeg_baseline';
signal_process                  = 'CSP';

EEG_sub = struct();
for indsub=1:9
    %% Extract and Arrange Data
    par.extractSound.signal_name    = signal_name;
    par.extractSound.InField        = 'train';
    par.extractSound.it_end         = 1;
    par.extractSound.multiEpoch     = true;
    [EEG_sub(indsub).trials,fsample]            = extractSound(indsub,par.extractSound);
end
EEG_trials = vertcat(EEG_sub.trials);
for iTr = 1:length(EEG_trials)
    EEG_trials(iTr).trialId = iTr;
end
% remapTypes
par.remapTypes           = remapTypesParams();
par.remapTypes.selection = {1,2};

StartClass = unique([EEG_trials.trialType]);
% % Time Interpolation and selection Trials
% par.TimeSelect               = TimeSelectParams;
% par.TimeSelect.t1            = 0.0; % in s from ZeroEvent time
% par.TimeSelect.t2            = 2.5; % in s from ZeroEvent time
% par.TimeSelect.InField       = signal_name;
% par.TimeSelect.OutField      = signal_name;

% itr1 = par.TimeSelect.t1;
% itr2 = par.TimeSelect.t2;

% Filter Bank
par.FilterBankCompute            = FilterBankComputeParams();
par.FilterBankCompute.InField    = signal_name;
par.FilterBankCompute.OutField   = signal_name;
% par.FilterBankCompute.f_min      = 0.1; % min frequency range in Hz
% par.FilterBankCompute.f_max      = 49.5;
par.FilterBankCompute.FilterBank = 'One';
par.FilterBankCompute.fsample    = fsample;

% %% epochCompute
% par.epochCompute                    = epochComputeParams();
% par.epochCompute.InField            = signal_name;
% par.epochCompute.OutField           = signal_name;
% par.epochCompute.fample             = fsample;
% par.epochCompute.t_epoch            = 0.5; % duration of a single intervals in s
% par.epochCompute.overlap_percent    = 0; % in percentage


par.exec.funname ={'remapTypes','FilterBankCompute'};
[EEG_trials,out.epochCompute] = run_trials(EEG_trials,par);
itr1 = round(out.epochCompute.epochCompute.time_intervals(:,1),2);
itr2 = round(out.epochCompute.epochCompute.time_intervals(:,2),2);

par.multiEEG.Infield = signal_name;
EEG_trials = multiEEG(EEG_trials,par.multiEEG);

% kfold-CrossValidation on the Train dataset
kfoldSplit = 10;
labs = [EEG_trials.trialType]'; %true labels
cvp = cvpartition(labs,'kfold',kfoldSplit,'Stratify',true);

resQDA          = struct();
resKNN          = struct();
resNB         = struct();

outQDA          = struct();
outKNN          = struct();
outNB           = struct();

LabpredictQDA   = struct();
LabpredictKNN   = struct();
LabpredictNB  = struct();
Label_train     = struct();

for i=1:kfoldSplit
    indices = training(cvp,i);
    test = (indices == 0);
    train = ~test;

    EEG_train = EEG_trials(train);
    EEG_test = EEG_trials(test);

    % par.multiEEG.Infield = signal_name;
    % EEG_train = multiEEG(EEG_trials(train),par.multiEEG);
    % EEG_test = multiEEG(EEG_trials(test),par.multiEEG);

    Label_train(i).Iter = [EEG_test.trialType]';
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
    EEG_train = run_trials(EEG_train,par);
    EEG_test = run_trials(EEG_test,par);

    TotalFeatures = size(EEG_test(1).(signal_process),2);

    par.miModel.k = 0;
    % % Mutual Information
    % par.miModel               = miModelParams;
    % par.miModel.InField       = signal_process;
    % par.miModel.m             = par.cspModel.m;
    % % par.miModel.k             = 5;
    % 
    % [~, out.miModel]=miModel(EEG_train,par.miModel);
    % 
    % par.miEncode               = miEncodeParams;
    % par.miEncode.InField       = signal_process;
    % par.miEncode.OutField      = signal_process;
    % par.miEncode.IndMI         = out.miModel.IndMI;
    % 
    % par.exec.funname ={'miEncode'};
    % [EEG_train, out]=run_trials(EEG_train,par);
    % V_train = out.miEncode.V;
    % [EEG_test, out]=run_trials(EEG_test,par);
    % V_test = out.miEncode.V;

    %% Step 3. Model Classification on CSP

    % qdaModel
    par.qdaModel                      = qdaModelParams;
    par.qdaModel.InField              = 'CSP';
    par.qdaModel.numIterations        = 100;
    par.qdaModel.kfold                = 5;
    [~, outQDA(i).Iter]               = qdaModel(EEG_train,par.qdaModel);

    % predictQDA
    par.mdlPredict                  = mdlPredictParams;
    par.mdlPredict.InField          = 'CSP';
    par.mdlPredict.OutField         = 'QDApred';
    par.mdlPredict.ProbField        = 'QDAProb';
    par.mdlPredict.mdl              = outQDA(i).Iter.mdl;
    [EEG_train,resQDA(i).train]     = mdlPredict(EEG_train,par.mdlPredict);
    [EEG_test,resQDA(i).test]       = mdlPredict(EEG_test,par.mdlPredict);

    LabpredictQDA(i).Iter = [EEG_test.QDApred]';

    % knnModel
    par.knnModel                      = knnModelParams;
    par.knnModel.InField              = 'CSP';
    par.knnModel.numIterations        = 100;
    par.knnModel.kfold                = 5;
    [~, outKNN(i).Iter]               = knnModel(EEG_train,par.knnModel);

    % predictKNN
    par.mdlPredict                  = mdlPredictParams;
    par.mdlPredict.InField          = 'CSP';
    par.mdlPredict.OutField         = 'KNNpred';
    par.mdlPredict.ProbField        = 'KNNProb';
    par.mdlPredict.mdl              = outKNN(i).Iter.mdl;
    [EEG_train,resKNN(i).train]     = mdlPredict(EEG_train,par.mdlPredict);
    [EEG_test,resKNN(i).test]       = mdlPredict(EEG_test,par.mdlPredict);

    LabpredictKNN(i).Iter = [EEG_test.KNNpred]';

    % nbModel
    par.nbModel                      = nbModelParams;
    par.nbModel.InField              = 'CSP';
    par.nbModel.numIterations        = 100;
    par.nbModel.kfold                = 5;
    [~, outNB(i).Iter]               = nbModel(EEG_train,par.nbModel);

    % predictNB
    par.mdlPredict                  = mdlPredictParams;
    par.mdlPredict.InField          = 'CSP';
    par.mdlPredict.OutField         = 'NBpred';
    par.mdlPredict.ProbField        = 'NBProb';
    par.mdlPredict.mdl              = outNB(i).Iter.mdl;
    [EEG_train,resNB(i).train]     = mdlPredict(EEG_train,par.mdlPredict);
    [EEG_test,resNB(i).test]       = mdlPredict(EEG_test,par.mdlPredict);

    LabpredictNB(i).Iter = [EEG_test.NBpred]';

end
%% Accuracy and Kappa
label_train = struct2cell(Label_train(:))'; % convert struct in mat
label_train = cell2mat(label_train);

predictQDA_train = struct2cell(LabpredictQDA(:))';
predictQDA_train = cell2mat(predictQDA_train);
predictKNN_train = struct2cell(LabpredictKNN(:))';
predictKNN_train = cell2mat(predictKNN_train);
predictNB_train = struct2cell(LabpredictNB(:))';
predictNB_train = cell2mat(predictNB_train);

AccuracyQDA = sum(predictQDA_train == label_train)/length(label_train)*100;
accuracyQDA_class = accuracy4classes(label_train,predictQDA_train);

AccuracyKNN = sum(predictKNN_train == label_train)/length(label_train)*100;
accuracyKNN_class = accuracy4classes(label_train,predictKNN_train);

AccuracyNB = sum(predictNB_train == label_train)/length(label_train)*100;
accuracyNB_class = accuracy4classes(label_train,predictNB_train);

CmatrxixQDA_train = confusionmat(label_train, predictQDA_train);
kappaQDA = kappaModel(CmatrxixQDA_train);

CmatrxixKNN_train = confusionmat(label_train, predictKNN_train);
kappaKNN = kappaModel(CmatrxixKNN_train);

CmatrxixNB_train = confusionmat(label_train, predictNB_train);
kappaNB= kappaModel(CmatrxixNB_train);

% Save Result
%% create Tab Result
params.createStructResult.subj       = indsub;
params.createStructResult.method     = 'CSP';
params.createStructResult.file       = 'Sound';
params.createStructResult.train_name = 'Strain';
params.createStructResult.train_tr1  = itr1;
params.createStructResult.train_tr2  = itr2;
params.createStructResult.test_name  = 'Strain';
params.createStructResult.test_ts1   = itr1;
params.createStructResult.test_ts2   = itr2;
params.createStructResult.m             = par.cspModel.m;
params.createStructResult.class         = findclass(EEG_train,StartClass);
params.createStructResult.irng          = par.irng;
params.createStructResult.Filter        = par.FilterBankCompute.FilterBank;
params.createStructResult.n_Features    = size(EEG_train(1).(signal_process),2);
params.createStructResult.indMi         = par.miModel.k;
params.createStructResult.attenuation   = par.FilterBankCompute.attenuation;
params.createStructResult.TotalFeatures = TotalFeatures;
params.createStructResult.kfold         = kfoldSplit;


% QDA save result
resultQDA.train.Accuracy = AccuracyQDA;
resultQDA.train.Accuracy_class = accuracyQDA_class;
resultQDA.test.Accuracy = NaN;
resultQDA.test.Accuracy_class = NaN;
resultQDA.train.kappaValue = kappaQDA;
resultQDA.test.kappaValue = NaN;

[ResultQDA_Kappa,ResultQDA_Acc,ResultQDA_class_Acc] = createStructResult(resultQDA,params.createStructResult);

% Update Tab Result
params.updateTab.dir        = 'D:\TrialBox_Results_excel\Sound_dataset';
params.updateTab.name       = 'Baseline_NoNSA_AllSubj';
params.updateTab.sheetnames = 'QDA';

updated_Result_tableAccQDA = updateTab(ResultQDA_Acc,params.updateTab);

params.updateTab.sheetnames = 'KappaQDA';
updated_Result_tableKappaQDA = updateTab(ResultQDA_Kappa,params.updateTab);

params.updateTab.name     = 'Baseline_class_NoNSA_AllSubj';
params.updateTab.sheetnames = 'QDA';
updated_Resultclass_tableAccQDA = updateTab(ResultQDA_class_Acc,params.updateTab);

% KNN save result
resultKNN.train.Accuracy = AccuracyKNN;
resultKNN.train.Accuracy_class = accuracyKNN_class;
resultKNN.test.Accuracy = NaN;
resultKNN.test.Accuracy_class = NaN;
resultKNN.train.kappaValue = kappaKNN;
resultKNN.test.kappaValue = NaN;

[ResultKNN_Kappa,ResultKNN_Acc,ResultKNN_class_Acc] = createStructResult(resultKNN,params.createStructResult);

%% Update Tab Result
params.updateTab.dir        = 'D:\TrialBox_Results_excel\Sound_dataset';
params.updateTab.name       = 'Baseline_NoNSA_AllSubj';
params.updateTab.sheetnames = 'KNN';
updated_Result_tableAccKNN = updateTab(ResultKNN_Acc,params.updateTab);

params.updateTab.sheetnames = 'KappaKNN';
updated_Result_tableKappaKNN = updateTab(ResultKNN_Kappa,params.updateTab);

params.updateTab.name     = 'Baseline_class_NoNSA_AllSubj';
params.updateTab.sheetnames = 'KNN';
updated_Resultclass_tableAccKNN = updateTab(ResultKNN_class_Acc,params.updateTab);

%  Acc
resultNB.train.Accuracy = AccuracyNB;
resultNB.train.Accuracy_class = accuracyNB_class;
resultNB.test.Accuracy = NaN;
resultNB.test.Accuracy_class = NaN;
resultNB.train.kappaValue = kappaNB;
resultNB.test.kappaValue = NaN;

[ResultNB_Kappa,ResultNB_Acc,ResultNB_class_Acc] = createStructResult(resultNB,params.createStructResult);

%% Update Tab Result
params.updateTab.dir        = 'D:\TrialBox_Results_excel\Sound_dataset';
params.updateTab.name       = 'Baseline_NoNSA_AllSubj';
params.updateTab.sheetnames = 'NB';
updated_Result_tableAccNB = updateTab(ResultNB_Acc,params.updateTab);

params.updateTab.sheetnames = 'KappaNB';
updated_Result_tableKappaNB = updateTab(ResultNB_Kappa,params.updateTab);

params.updateTab.name     = 'Baseline_class_NoNSA_AllSubj';
params.updateTab.sheetnames = 'NB';
updated_Resultclass_tableAccNB = updateTab(ResultNB_class_Acc,params.updateTab);
