% TEST_Sound_CrossVal_OnTrain.m


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

indsub = 9;
signal_name                     = 'eeg_sound';
signal_process                  = 'tsne';

%% Extract and Arrange Data
par.extractSound.signal_name    = signal_name;
par.extractSound.InField        = 'train';
[EEG_trials,fsample]            = extractSound(indsub,par.extractSound);

% remapTypes
par.remapTypes           = remapTypesParams();
par.remapTypes.selection = {1,2};

StartClass = unique([EEG_trials.trialType]);
% Time Interpolation and selection Trials [0.5;2.5] from CUE (Motor Imagery Interval)
par.TimeSelect               = TimeSelectParams;
par.TimeSelect.t1            = 0.5; % in s from ZeroEvent time
par.TimeSelect.t2            = 2.5; % in s from ZeroEvent time
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

par.exec.funname ={'remapTypes','TimeSelect','FilterBankCompute'};
EEG_trials =run_trials(EEG_trials,par);


% kfold-CrossValidation on the Train dataset
kfold = 3;
labs = [EEG_trials.trialType]'; %true labels
cvp = cvpartition(labs,'kfold',kfold,'Stratify',true);

resQDA          = struct();
resKNN          = struct();
resNBPW         = struct();
res_kappaNBPW   = struct();

outQDA          = struct();
outKNN          = struct();

LabpredictQDA   = struct();
LabpredictKNN   = struct();
LabpredictNBPW  = struct();
Label_train     = struct();

for i=1:kfold
    indices = training(cvp,i);
    test = (indices == 0);
    train = ~test;
    EEG_train = EEG_trials(train);

    % Bootstrap
    par.BootStrapData               = BootStrapDataParams;
    par.BootStrapData.N             = 100;
    par.BootStrapData.InField       = signal_name;
    par.BootStrapData.OutField      = signal_name;
    EEG_train                       = BootStrapData(EEG_train,par.BootStrapData);
    EEG_train = EEG_train';

    EEG_test = EEG_trials(test);
    Label_train(i).Iter = [EEG_test.trialType]';

    %% Step 2. perform tsne
    par.tsneModel.NumDimensions   = 2;
    par.tsneModel.InField         = signal_name;
    par.tsneModel.OutField        = 'tsne';
    par.tsneModel.xfld            = 'time';
    par.exec.funname ={'tsneModel'};
    EEG_train = run_trials(EEG_train,par);
    EEG_test = run_trials(EEG_test,par);


    %% Step 3. Model Classification on tsne
    % qdaModel
    par.qdaModel                      = qdaModelParams;
    par.qdaModel.InField              = 'tsne';
    par.qdaModel.numIterations        = 100;
    par.qdaModel.kfold                = 5;
    [~, outQDA(i).Iter]               = qdaModel(EEG_train,par.qdaModel);

    % predictQDA
    par.mdlPredict                  = mdlPredictParams;
    par.mdlPredict.InField          = 'tsne';
    par.mdlPredict.OutField         = 'QDApred';
    par.mdlPredict.ProbField        = 'QDAProb';
    par.mdlPredict.mdl              = outQDA(i).Iter.mdl;
    [EEG_train,resQDA(i).train]     = mdlPredict(EEG_train,par.mdlPredict);
    [EEG_test,resQDA(i).test]       = mdlPredict(EEG_test,par.mdlPredict);

    LabpredictQDA(i).Iter = [EEG_test.QDApred]';

    % knnModel
    par.knnModel                      = knnModelParams;
    par.knnModel.InField              = 'tsne';
    par.knnModel.numIterations        = 100;
    par.knnModel.kfold                = 5;
    [~, outKNN(i).Iter]             = knnModel(EEG_train,par.knnModel);

    % predictKNN
    par.mdlPredict                  = mdlPredictParams;
    par.mdlPredict.InField          = 'tsne';
    par.mdlPredict.OutField         = 'KNNpred';
    par.mdlPredict.ProbField        = 'KNNProb';
    par.mdlPredict.mdl              = outKNN(i).Iter.mdl;
    [EEG_train,resKNN(i).train]     = mdlPredict(EEG_train,par.mdlPredict);
    [EEG_test,resKNN(i).test]       = mdlPredict(EEG_test,par.mdlPredict);

    LabpredictKNN(i).Iter = [EEG_test.KNNpred]';
end
%% Accuracy and Kappa
label_train = struct2cell(Label_train(:))'; % convert struct in mat
label_train = cell2mat(label_train);

predictQDA_train = struct2cell(LabpredictQDA(:))';
predictQDA_train = cell2mat(predictQDA_train);
predictKNN_train = struct2cell(LabpredictKNN(:))';
predictKNN_train = cell2mat(predictKNN_train);

AccuracyQDA = sum(predictQDA_train == label_train)/length(label_train)*100;
accuracyQDA_class = accuracy4classes(label_train,predictQDA_train);

AccuracyKNN = sum(predictKNN_train == label_train)/length(label_train)*100;
accuracyKNN_class = accuracy4classes(label_train,predictKNN_train);


CmatrxixQDA_train = confusionmat(label_train, predictQDA_train);
kappaQDA = kappaModel(CmatrxixQDA_train);

CmatrxixKNN_train = confusionmat(label_train, predictKNN_train);
kappaKNN = kappaModel(CmatrxixKNN_train);

% Save Result
%% create Tab Result
params.createStructResult.subj       = indsub;
params.createStructResult.method     = 'tsne';
params.createStructResult.file       = 'Sound';
params.createStructResult.train_name = 'Strain';
params.createStructResult.train_tr1  = itr1;
params.createStructResult.train_tr2  = itr2;
params.createStructResult.test_name  = 'Strain';
params.createStructResult.test_ts1   = itr1;
params.createStructResult.test_ts2   = itr2;
params.createStructResult.class      = findclass(EEG_train,StartClass);
params.createStructResult.irng       = par.irng;
params.createStructResult.method     = signal_process;
params.createStructResult.Filter     = par.FilterBankCompute.FilterBank;

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
params.updateTab.name       = 'Sound_CrossVal_OnTrain_tsne';
params.updateTab.sheetnames = 'QDA';

updated_Result_tableAccQDA = updateTab(ResultQDA_Acc,params.updateTab);

params.updateTab.sheetnames = 'KappaQDA';
updated_Result_tableKappaQDA = updateTab(ResultQDA_Kappa,params.updateTab);

params.updateTab.name     = 'Sound_CrossVal_OnTrain_class_tsne';
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
params.updateTab.name       = 'Sound_CrossVal_OnTrain_tsne';
params.updateTab.sheetnames = 'KNN';
updated_Result_tableAccKNN = updateTab(ResultKNN_Acc,params.updateTab);

params.updateTab.sheetnames = 'KappaKNN';
updated_Result_tableKappaKNN = updateTab(ResultKNN_Kappa,params.updateTab);

params.updateTab.name     = 'Sound_CrossVal_OnTrain_class_tsne';
params.updateTab.sheetnames = 'KNN';
updated_Resultclass_tableAccKNN = updateTab(ResultKNN_class_Acc,params.updateTab);
