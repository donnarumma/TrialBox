% TEST_Graz_CrossVal_OnTest_INTERVAL.m

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

for indsub = 1:9
    %% Step 0. Load raw data and arrange trials
    signal_name                     = 'eeg';
    signal_process                  = 'CSP';

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

    %% Extract and Arrange TEST Data

    Dt = [0.5 4];     % Interval of evaluation (Motor Imagery)
    TW = 2;         % time window
    dt = 1;         %

    num_intervals = floor((Dt(2) - Dt(1) - TW) / dt) + 1;
    outNBPW = struct();
    resQDA = struct();
    resKNN = struct();
    P = struct();
    for n_val = 1:num_intervals

        start_time = Dt(1) + (n_val - 1) * dt;
        end_time = start_time + TW;

        par.extractGraz.signal_name        = signal_name;
        par.extractGraz.InField            = 'test';
        [EEG_test,fsample]                 = extractGraz(indsub,par.extractGraz);

        % Time Interpolation and selection Trials [start_time;end_time] from CUE (Motor Imagery Interval)
        par.TimeSelect               = TimeSelectParams;
        par.TimeSelect.t1            = start_time; % in s from ZeroEvent time
        par.TimeSelect.t2            = end_time; % in s from ZeroEvent time
        par.TimeSelect.InField       = signal_name;
        par.TimeSelect.OutField      = signal_name;
        par.TimeSelect.dt            = 1;

        its1 = par.TimeSelect.t1;
        its2 = par.TimeSelect.t2;

        % Filter Bank
        par.FilterBankCompute            = FilterBankComputeParams();
        par.FilterBankCompute.InField    = signal_name;
        par.FilterBankCompute.OutField   = signal_name;
        par.FilterBankCompute.FilterBank = 'Nine';
        par.FilterBankCompute.fsample    = fsample;

        par.exec.funname ={'TimeSelect','FilterBankCompute'};
        [EEG_test, par.execinfo]=run_trials(EEG_test,par);

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


        % Mutual Information
        par.miModel               = miModelParams;
        par.miModel.InField       = signal_process;

        [~, out.miModel] = miModel(EEG_train,par.miModel);

        par.miEncode               = miEncodeParams;
        par.miEncode.InField       = signal_process;
        par.miEncode.OutField      = signal_process;
        par.miEncode.IndMI         = out.miModel.IndMI;

        par.exec.funname ={'miEncode'};
        [EEG_train, out]=run_trials(EEG_train,par);
        V_train = out.miEncode.V;
        [EEG_test, out]=run_trials(EEG_test,par);
        V_test = out.miEncode.V;

        %% Step 3. QDA Classification on CSP
        % qdaModel
        par.qdaModel                        = qdaModelParams;
        par.qdaModel.InField                = 'CSP';
        par.qdaModel.numIterations          = 100;
        par.qdaModel.kfold                  = 5;
        [~, outQDA]                         = qdaModel(EEG_train,par.qdaModel);

        % predictQDA
        par.mdlPredict                      = mdlPredictParams;
        par.mdlPredict.InField              = 'CSP';
        par.mdlPredict.OutField             = 'QDApred';
        par.mdlPredict.ProbField            = 'QDAProb';
        par.mdlPredict.mdl                  = outQDA.mdl;
        [EEG_train,resQDA(n_val).train]     = mdlPredict(EEG_train,par.mdlPredict);
        [EEG_test,resQDA(n_val).test]       = mdlPredict(EEG_test,par.mdlPredict);

        % knnModel
        par.knnModel                        = knnModelParams;
        par.knnModel.InField                = 'CSP';
        par.knnModel.numIterations          = 100;
        par.knnModel.kfold                  = 5;
        [~, outKNN]                         = knnModel(EEG_train,par.knnModel);

        % predictKNN
        par.mdlPredict                      = mdlPredictParams;
        par.mdlPredict.InField              = 'CSP';
        par.mdlPredict.OutField             = 'KNNpred';
        par.mdlPredict.ProbField            = 'KNNProb';
        par.mdlPredict.mdl                  = outKNN.mdl;
        [EEG_train,resKNN(n_val).train]     = mdlPredict(EEG_train,par.mdlPredict);
        [EEG_test,resKNN(n_val).test]       = mdlPredict(EEG_test,par.mdlPredict);

        for nTr = 1:length(EEG_test)
            P(n_val).QDA.train.prob(nTr,:) = EEG_train(nTr).QDAProb';
            P(n_val).QDA.test.prob(nTr,:) = EEG_test(nTr).QDAProb';
            P(n_val).KNN.train.prob(nTr,:) = EEG_train(nTr).KNNProb';
            P(n_val).KNN.test.prob(nTr,:) = EEG_test(nTr).KNNProb';
        end


        % NBPW
        par.fitNBPW                         = fitNBPWParams;
        par.fitNBPW.labs_train              = [EEG_train.trialType]';
        par.fitNBPW.labs_test               = [EEG_train.trialType]';
        [pred_train,outNBPW(n_val).train]   = fitNBPW(V_train,V_train,par.fitNBPW);

        par.fitNBPW                         = fitNBPWParams;
        par.fitNBPW.labs_train              = [EEG_train.trialType]';
        par.fitNBPW.labs_test               = [EEG_test.trialType]';
        [pred_test,outNBPW(n_val).test]     = fitNBPW(V_train,V_test,par.fitNBPW);

        P(n_val).NBPW.train.prob = outNBPW(n_val).train.prob;
        P(n_val).NBPW.test.prob = outNBPW(n_val).test.prob;
    end

    par.probPredict                 = probPredictParams;
    par.probPredict.InField         = 'QDA';
    par.probPredict.OutField        = 'ResultQDA';
    par.probPredict.label_train     = [EEG_train.trialType]';
    par.probPredict.label_test      = [EEG_test.trialType]';
    ResultQDA                       = probPredict(P,par.probPredict);

    par.probPredict.InField         = 'KNN';
    ResultKNN                       = probPredict(P,par.probPredict);

    par.probPredict.InField         = 'NBPW';
    ResultNBPW                      = probPredict(P,par.probPredict);

    % Create Tab and save Results
    params.createStructInterval.subj       = indsub;
    params.createStructInterval.method     = 'CSP';
    params.createStructInterval.file       = 'Graz';
    params.createStructInterval.train_name = 'AT';
    params.createStructInterval.train_tr1  = itr1;
    params.createStructInterval.train_tr2  = itr2;
    params.createStructInterval.test_name  = 'AE';
    params.createStructInterval.test_ts1   = its1;
    params.createStructInterval.test_ts2   = its2;
    params.createStructInterval.m          = par.cspModel.m;
    params.createStructInterval.class      = findclass(EEG_train,StartClass);
    params.createStructInterval.irng       = par.irng;
    params.createStructInterval.method     = signal_process;

    % QDA Accuracy save
    [ResultQDA_kappa,ResultQDA_Acc]     = createStructInterval(ResultQDA,params.createStructInterval);
    % Update Tab Result
    params.updateTab.dir                = 'D:\TrialBox_Results_excel\Graz_dataset';
    params.updateTab.name               = 'Graz_OnTest_Interval';
    params.updateTab.sheetnames         = 'QDA';
    updated_Result_tableAccQDA          = updateTab(ResultQDA_Acc,params.updateTab);

    params.updateTab.sheetnames         = 'KappaQDA';
    updated_Result_tableKappaQDA        = updateTab(ResultQDA_kappa,params.updateTab);

    % KNN Accuracy save
    [ResultKNN_kappa,ResultKNN_Acc]     = createStructInterval(ResultKNN,params.createStructInterval);
    % Update Tab Result
    params.updateTab.dir                = 'D:\TrialBox_Results_excel\Graz_dataset';
    params.updateTab.name               = 'Graz_OnTest_Interval';
    params.updateTab.sheetnames         = 'KNN';
    updated_Result_tableAccKNN          = updateTab(ResultKNN_Acc,params.updateTab);

    params.updateTab.sheetnames         = 'KappaKNN';
    updated_Result_tableKappaKNN        = updateTab(ResultKNN_kappa,params.updateTab);

    % NBPW Accuracy save
    [ResultNBPW_kappa,ResultNBPW_Acc]   = createStructInterval(ResultNBPW,params.createStructInterval);
    % Update Tab Result
    params.updateTab.dir                = 'D:\TrialBox_Results_excel\Graz_dataset';
    params.updateTab.name               = 'Graz_OnTest_Interval';
    params.updateTab.sheetnames         = 'NBPW';
    updated_Result_tableAccNBPW         = updateTab(ResultNBPW_Acc,params.updateTab);

    params.updateTab.sheetnames         = 'KappaNBPW';
    updated_Result_tableKappaNBPW       = updateTab(ResultNBPW_kappa,params.updateTab);
end