% TEST_Graz_OnTest_NoNSA.m

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

% par.irng = 10;
% rng(par.irng);
subj=1:9;
for i=randi(50,1,5)
    par.irng = i;
rng(par.irng);
    for indsub=2:9

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
                % svcModel
                par.svcModel                      = svcModelParams;
                par.svcModel.InField              = 'CSP';
                par.svcModel.numIterations        = 100;
                par.svcModel.kfold                = 5;
                [~, outSVC]                       = svcModel(EEG_train,par.svcModel);

                % predictSVC
                par.mdlPredict                  = mdlPredictParams;
                par.mdlPredict.InField          = 'CSP';
                par.mdlPredict.OutField         = 'SVCpred';
                par.mdlPredict.ProbField        = 'SVCProb';
                par.mdlPredict.mdl              = outSVC.mdl;
                [EEG_train,resSVC.train]        = mdlPredict(EEG_train,par.mdlPredict);
                [EEG_test,resSVC.test]          = mdlPredict(EEG_test,par.mdlPredict);
% Save Result
                % create Tab Result
                params.createStructResult.subj       = indsub;
                params.createStructResult.method     = 'CSP';
                params.createStructResult.file       = 'Graz';
                params.createStructResult.train_name = strcat('AT','_',string(indsub));
                params.createStructResult.train_tr1  = itr1;
                params.createStructResult.train_tr2  = itr2;
                params.createStructResult.test_name  = strcat('AE','_',string(indsub));
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

                % SVC Accuracy save
                [ResultSVC_kappa,ResultSVC_Acc,ResultSVC_class_Acc] = createStructResult(resSVC,params.createStructResult);
                % Update Tab Result
                params.updateTab.dir                = 'D:\TrialBox_Results_excel\Graz_dataset';
                params.updateTab.name               = 'Graz_OnTest_NoNSA';
                params.updateTab.sheetnames         = 'SVC';
                updated_Result_tableAccSVC          = updateTab(ResultSVC_Acc,params.updateTab);

                params.updateTab.sheetnames         = 'KappaSVC';
                updated_Result_tableKappaSVC        = updateTab(ResultSVC_kappa,params.updateTab);

                params.updateTab.name               = 'Graz_OnTest_NoNSA_class';
                params.updateTab.sheetnames         = 'SVC';
                updated_Resultclass_tableAccSVC     = updateTab(ResultSVC_class_Acc,params.updateTab);
    end
end