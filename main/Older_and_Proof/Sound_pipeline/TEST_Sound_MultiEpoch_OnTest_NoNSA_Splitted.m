% TEST_Graz_OnTest_NoNSA_OnTest.m

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
subj = 1;

for ir=10
    par.irng = ir;
    rng(par.irng);
    for indsub=2:9

        signal_name                     = 'eeg_sound';
        signal_process                  = 'CSP';

        %% Extract and Arrange Data
        par.extractSound.signal_name    = signal_name;
        par.extractSound.InField        = 'train';
        par.extractSound.it_end         = 2.5;
        par.extractSound.multiEpoch     = true;
        [EEG_train,fsample]            = extractSound(indsub,par.extractSound);

        % remapTypes
        par.remapTypes           = remapTypesParams();
        par.remapTypes.selection = {1,2};

        StartClass = unique([EEG_train.trialType]);

        % Time Interpolation and selection Trials [0.5;2.5] from CUE (Motor Imagery Interval)
        par.TimeSelect               = TimeSelectParams;
        par.TimeSelect.t1            = 0.5; % in s from ZeroEvent time
        par.TimeSelect.t2            = 1.5; % in s from ZeroEvent time
        par.TimeSelect.InField       = signal_name;
        par.TimeSelect.OutField      = signal_name;
        % par.TimeSelect.dt            = 1;


        % Filter Bank
        par.FilterBankCompute            = FilterBankComputeParams();
        par.FilterBankCompute.InField    = signal_name;
        par.FilterBankCompute.OutField   = signal_name;
        % par.FilterBankCompute.attenuation = 10;
        par.FilterBankCompute.FilterBank = 'One';
        par.FilterBankCompute.fsample    = fsample;

        %% epochCompute
        par.epochCompute                    = epochComputeParams();
        par.epochCompute.InField            = signal_name;
        par.epochCompute.OutField           = signal_name;
        par.epochCompute.fample             = fsample;
        par.epochCompute.t_epoch            = 0.5; % duration of a single intervals in s
        par.epochCompute.overlap_percent    = 50; % in percentage


        par.exec.funname ={'remapTypes','TimeSelect','FilterBankCompute','epochCompute'};
        [EEG_train,out.epochCompute] =run_trials(EEG_train,par);

        itr1 = round(out.epochCompute.epochCompute.time_intervals(:,1),2);
        itr2 = round(out.epochCompute.epochCompute.time_intervals(:,2),2);

        %% Extract and Arrange Data
        par.extractSound.signal_name    = signal_name;
        par.extractSound.InField        = 'test';
        par.extractSound.it_end         = 2.5;
        par.extractSound.multiEpoch     = true;
        [EEG_test,fsample]            = extractSound(indsub,par.extractSound);

        % remapTypes
        par.remapTypes           = remapTypesParams();
        par.remapTypes.selection = {1,2};

        % Time Interpolation and selection Trials [0.5;2.5] from CUE (Motor Imagery Interval)
        par.TimeSelect               = TimeSelectParams;
        par.TimeSelect.t1            = 1.5; % in s from ZeroEvent time
        par.TimeSelect.t2            = 2.5; % in s from ZeroEvent time
        par.TimeSelect.InField       = signal_name;
        par.TimeSelect.OutField      = signal_name;
        % par.TimeSelect.dt            = 1;



        par.exec.funname ={'remapTypes','TimeSelect','FilterBankCompute','epochCompute'};

        [EEG_test,out.epochCompute] =run_trials(EEG_test,par);

        its1 = round(out.epochCompute.epochCompute.time_intervals(:,1),2);
        its2 = round(out.epochCompute.epochCompute.time_intervals(:,2),2);

        resQDA          = struct();
        resKNN          = struct();
        resNB           = struct();
        resSVC          = struct();

        outQDA          = struct();
        outKNN          = struct();
        outNB           = struct();
        outSVC          = struct();

        LabpredictQDA   = struct();
        LabpredictKNN   = struct();
        LabpredictNB  = struct();
        LabpredictSVC  = struct();

        Label_test     = struct();

        ksplit = 5;
        labs = [EEG_train.trialType]'; %true labels
        cvp = cvpartition(labs,'kfold',ksplit,'Stratify',true);

        for i=1:ksplit
            indices = training(cvp,i);
            test = (indices == 0);
            train = ~test;
            %% multiEEG
            par.multiEEG.Infield = signal_name;
            EEG_trainM = multiEEG(EEG_train(train),par.multiEEG);
            EEG_testM = multiEEG(EEG_test(test),par.multiEEG);

            Label_test(i).Iter = [EEG_testM.trialType]';
            %% Step 2. perform CSP
            % CSP Dictionary evaluation on train
            par.cspModel                  = cspModelParams;
            par.cspModel.m                = 2;
            par.cspModel.InField          = signal_name;
            par.cspModel.OutField         = signal_process;

            [~,out.cspModel] = cspModel(EEG_trainM,par.cspModel);

            % CSP Encode on train and test data
            par.cspEncode                  = cspEncodeParams;
            par.cspEncode.InField          = signal_name;
            par.cspEncode.OutField         = signal_process;
            par.cspEncode.W                = out.cspModel.W;

            par.exec.funname ={'cspEncode'};
            EEG_trainM =run_trials(EEG_trainM,par);
            EEG_testM =run_trials(EEG_testM,par);

            TotalFeatures = size(EEG_testM(1).(signal_process),2);

            % Mutual Information
            par.miModel               = miModelParams;
            par.miModel.InField       = signal_process;
            par.miModel.m             = par.cspModel.m;

            [~, out.miModel]=miModel(EEG_trainM,par.miModel);

            par.miEncode               = miEncodeParams;
            par.miEncode.InField       = signal_process;
            par.miEncode.OutField      = signal_process;
            par.miEncode.IndMI         = out.miModel.IndMI;

            par.exec.funname ={'miEncode'};
            [EEG_trainM, out]=run_trials(EEG_trainM,par);
            V_train = out.miEncode.V;
            [EEG_testM, out]=run_trials(EEG_testM,par);
            V_test = out.miEncode.V;

            %% Step 3. Model Classification on CSP

            % qdaModel
            par.qdaModel                      = qdaModelParams;
            par.qdaModel.InField              = 'CSP';
            par.qdaModel.numIterations        = 100;
            par.qdaModel.kfold                = 5;
            [~, outQDA(i).Iter]               = qdaModel(EEG_trainM,par.qdaModel);

            % predictQDA
            par.mdlPredict                  = mdlPredictParams;
            par.mdlPredict.InField          = 'CSP';
            par.mdlPredict.OutField         = 'QDApred';
            par.mdlPredict.ProbField        = 'QDAProb';
            par.mdlPredict.mdl              = outQDA(i).Iter.mdl;
            [EEG_trainM,resQDA(i).train]     = mdlPredict(EEG_trainM,par.mdlPredict);
            [EEG_testM,resQDA(i).test]       = mdlPredict(EEG_testM,par.mdlPredict);

            LabpredictQDA(i).Iter = [EEG_testM.QDApred]';

            % knnModel
            par.knnModel                      = knnModelParams;
            par.knnModel.InField              = 'CSP';
            par.knnModel.numIterations        = 100;
            par.knnModel.kfold                = 5;
            [~, outKNN(i).Iter]               = knnModel(EEG_trainM,par.knnModel);

            % predictKNN
            par.mdlPredict                  = mdlPredictParams;
            par.mdlPredict.InField          = 'CSP';
            par.mdlPredict.OutField         = 'KNNpred';
            par.mdlPredict.ProbField        = 'KNNProb';
            par.mdlPredict.mdl              = outKNN(i).Iter.mdl;
            [EEG_trainM,resKNN(i).train]     = mdlPredict(EEG_trainM,par.mdlPredict);
            [EEG_testM,resKNN(i).test]       = mdlPredict(EEG_testM,par.mdlPredict);

            LabpredictKNN(i).Iter = [EEG_testM.KNNpred]';

            % svcModel
            par.svcModel                      = svcModelParams;
            par.svcModel.InField              = 'CSP';
            par.svcModel.numIterations        = 100;
            par.svcModel.kfold                = 5;
            [~, outSVC(i).Iter]               = svcModel(EEG_trainM,par.svcModel);

            % predictSVC
            par.mdlPredict                  = mdlPredictParams;
            par.mdlPredict.InField          = 'CSP';
            par.mdlPredict.OutField         = 'SVCpred';
            par.mdlPredict.ProbField        = 'SVCProb';
            par.mdlPredict.mdl              = outSVC(i).Iter.mdl;
            [EEG_trainM,resSVC(i).train]     = mdlPredict(EEG_trainM,par.mdlPredict);
            [EEG_testM,resSVC(i).test]       = mdlPredict(EEG_testM,par.mdlPredict);

            LabpredictSVC(i).Iter = [EEG_testM.SVCpred]';

            % nbModel
            par.nbModel                      = nbModelParams;
            par.nbModel.InField              = 'CSP';
            par.nbModel.numIterations        = 100;
            par.nbModel.kfold                = 5;
            [~, outNB(i).Iter]               = nbModel(EEG_trainM,par.nbModel);

            % predictNB
            par.mdlPredict                  = mdlPredictParams;
            par.mdlPredict.InField          = 'CSP';
            par.mdlPredict.OutField         = 'NBpred';
            par.mdlPredict.ProbField        = 'NBProb';
            par.mdlPredict.mdl              = outNB(i).Iter.mdl;
            [EEG_trainM,resNB(i).train]     = mdlPredict(EEG_trainM,par.mdlPredict);
            [EEG_testM,resNB(i).test]       = mdlPredict(EEG_testM,par.mdlPredict);

            LabpredictNB(i).Iter = [EEG_testM.NBpred]';
        end
        %% Accuracy and Kappa
        label_test = struct2cell(Label_test(:))'; % convert struct in mat
        label_test = cell2mat(label_test);

        predictQDA_test = struct2cell(LabpredictQDA(:))';
        predictQDA_test = cell2mat(predictQDA_test);
        predictKNN_test = struct2cell(LabpredictKNN(:))';
        predictKNN_test = cell2mat(predictKNN_test);
        predictNB_test = struct2cell(LabpredictNB(:))';
        predictNB_test = cell2mat(predictNB_test);

        predictSVC_test = struct2cell(LabpredictSVC(:))';
        predictSVC_test = cell2mat(predictSVC_test);

        AccuracyQDA = sum(predictQDA_test == label_test)/length(label_test)*100;
        accuracyQDA_class = accuracy4classes(label_test,predictQDA_test);

        AccuracyKNN = sum(predictKNN_test == label_test)/length(label_test)*100;
        accuracyKNN_class = accuracy4classes(label_test,predictKNN_test);


        AccuracySVC = sum(predictSVC_test == label_test)/length(label_test)*100;
        accuracySVC_class = accuracy4classes(label_test,predictSVC_test);

        AccuracyNB = sum(predictNB_test == label_test)/length(label_test)*100;
        accuracyNB_class = accuracy4classes(label_test,predictNB_test);

        CmatrxixQDA_test = confusionmat(label_test, predictQDA_test);
        kappaQDA = kappaModel(CmatrxixQDA_test);

        CmatrxixKNN_test = confusionmat(label_test, predictKNN_test);
        kappaKNN = kappaModel(CmatrxixKNN_test);

        CmatrxixNB_test = confusionmat(label_test, predictNB_test);
        kappaNB = kappaModel(CmatrxixNB_test);

        CmatrxixSVC_test = confusionmat(label_test, predictSVC_test);
        kappaSVC = kappaModel(CmatrxixSVC_test);

        %% create Tab Result
        params.createStructResult.subj       = indsub;
        params.createStructResult.method     = 'CSP';
        params.createStructResult.file       = 'Sound';
        params.createStructResult.train_name = 'Strain';
        params.createStructResult.train_tr1  = itr1;
        params.createStructResult.train_tr2  = itr2;
        params.createStructResult.test_name  = 'Stest';
        params.createStructResult.test_ts1   = its1;
        params.createStructResult.test_ts2   = its2;
        params.createStructResult.m             = par.cspModel.m;
        params.createStructResult.class         = findclass(EEG_train,StartClass);
        params.createStructResult.irng          = par.irng;
        params.createStructResult.Filter        = par.FilterBankCompute.FilterBank;
        params.createStructResult.n_Features    = size(EEG_testM(1).(signal_process),2);
        params.createStructResult.indMi         = par.miModel.k;
        params.createStructResult.attenuation   = par.FilterBankCompute.attenuation;
        params.createStructResult.TotalFeatures = TotalFeatures;
        params.createStructResult.kfold         = 0;

        % QDA Accuracy save
        [ResultQDA_kappa,ResultQDA_Acc,ResultQDA_class_Acc] = createStructResult(resQDA,params.createStructResult);
        % Update Tab Result
        params.updateTab.dir        = 'D:\TrialBox_Results_excel\Sound_dataset';
        params.updateTab.name       = 'Sound_OnTest_Split';
        params.updateTab.sheetnames = 'QDA';

        updated_Result_tableAccQDA = updateTab(ResultQDA_Acc,params.updateTab);

        params.updateTab.sheetnames = 'kappaQDA';
        updated_Result_tablekappaQDA = updateTab(ResultQDA_kappa,params.updateTab);

        params.updateTab.name     = 'Sound_Class_OnTest_Split';
        params.updateTab.sheetnames = 'QDA';
        updated_Resultclass_tableAccQDA = updateTab(ResultQDA_class_Acc,params.updateTab);

        % KNN save res
        [ResultKNN_kappa,ResultKNN_Acc,ResultKNN_class_Acc] = createStructResult(resKNN,params.createStructResult);

        %% Update Tab Result
        params.updateTab.dir        = 'D:\TrialBox_Results_excel\Sound_dataset';
        params.updateTab.name       = 'Sound_OnTest_Split';
        params.updateTab.sheetnames = 'KNN';
        updated_Result_tableAccKNN = updateTab(ResultKNN_Acc,params.updateTab);

        params.updateTab.sheetnames = 'kappaKNN';
        updated_Result_tablekappaKNN = updateTab(ResultKNN_kappa,params.updateTab);

        params.updateTab.name     = 'Sound_Class_OnTest_Split';
        params.updateTab.sheetnames = 'KNN';
        updated_Resultclass_tableAccKNN = updateTab(ResultKNN_class_Acc,params.updateTab);
        % SVC Acc
        [ResultSVC_kappa,ResultSVC_Acc,ResultSVC_class_Acc] = createStructResult(resSVC,params.createStructResult);

        %% Update Tab Result
        params.updateTab.dir        = 'D:\TrialBox_Results_excel\Sound_dataset';
        params.updateTab.name       = 'Sound_OnTest_Split';
        params.updateTab.sheetnames = 'SVC';
        updated_Result_tableAccSVC = updateTab(ResultSVC_Acc,params.updateTab);

        params.updateTab.sheetnames = 'kappaSVC';
        updated_Result_tablekappaSVC = updateTab(ResultSVC_kappa,params.updateTab);

        params.updateTab.name     = 'Sound_Class_OnTest_Split';
        params.updateTab.sheetnames = 'SVC';
        updated_Resultclass_tableAccSVC = updateTab(ResultSVC_class_Acc,params.updateTab);

        % NB Acc
        [ResultNB_kappa,ResultNB_Acc,ResultNB_class_Acc] = createStructResult(resNB,params.createStructResult);

        %% Update Tab Result
        params.updateTab.dir        = 'D:\TrialBox_Results_excel\Sound_dataset';
        params.updateTab.name       = 'Sound_OnTest_Split';
        params.updateTab.sheetnames = 'NB';
        updated_Result_tableAccNB = updateTab(ResultNB_Acc,params.updateTab);

        params.updateTab.sheetnames = 'kappaNB';
        updated_Result_tablekappaNB = updateTab(ResultNB_kappa,params.updateTab);

        params.updateTab.name     = 'Sound_Class_OnTest_Split';
        params.updateTab.sheetnames = 'NB';
        updated_Resultclass_tableAccNB = updateTab(ResultNB_class_Acc,params.updateTab);
    end
end