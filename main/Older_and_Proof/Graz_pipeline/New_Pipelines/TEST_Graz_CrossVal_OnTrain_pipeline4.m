% TEST_Graz_CrossVal_OnTrain_pipeline4.m

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

% subj = 1:9;

tsub1 = [0,1.5429,1.3018,1.0606,1.4626,1.2214,1.3822,1.3018,1.5429];

% case 1,4
tsub3 = [0.9952,1.5468,1.3047,1.1488,1.5508,0.9717,1.3004,1.8104,1.5468];

% % case 2
% tsub3 = [0.8998,1.5429,1.3018,1.141,1.5429,0.9802,1.3018,1.8645,0.7391];

% % case 3
% tsub3 = [0.8998,1.5429,1.3018,1.0606,1.4626,1.2214,1.3822,1.3018,1.5429];


for indsub = 2:9

    % T*2cl
    int4sub1 = tsub1(indsub) - 0.75;
    int4sub2 = tsub1(indsub) + 0.75;

    % % tscl
    % int4sub1 = 0.5;
    % int4sub2 = 2.5;

    for deltaT = 1:2

        % Case 1,2
        int4sub3 = tsub3(indsub) - 0.1*deltaT;
        int4sub4 = tsub3(indsub) + 0.1*deltaT;

        % % case 3
        % int4sub3 = tsub3(indsub) - 0.75;
        % int4sub4 = tsub3(indsub) + 0.75;
        %
        % % case 4
        % int4sub3 = 0.5;
        % int4sub4 = 1.65;

        signal_name                     = 'eeg';
        signal_process                  = 'pca';

        %% Extract and Arrange Data
        par.extractGraz.signal_name     = signal_name;
        par.extractGraz.InField         = 'train';
        [EEG_trials,fsample]            = extractGraz(indsub,par.extractGraz);

        StartClass = unique([EEG_trials.trialType]);
        % Time Interpolation and selection Trials [0.75;2.5] from CUE (Motor Imagery Interval)
        par.TimeSelect               = TimeSelectParams;
        par.TimeSelect.t1            = int4sub1; % in s from ZeroEvent time
        par.TimeSelect.t2            = int4sub2; % in s from ZeroEvent time
        par.TimeSelect.InField       = signal_name;
        par.TimeSelect.OutField      = signal_name;
        par.TimeSelect.dt            = 1;


        % Filter Bank
        par.FilterBankCompute               = FilterBankComputeParams();
        par.FilterBankCompute.InField       = signal_name;
        par.FilterBankCompute.OutField      = signal_name;
        par.FilterBankCompute.attenuation   = 10; % one HyperParam
        par.FilterBankCompute.FilterBank    = 'Nine';
        par.FilterBankCompute.fsample       = fsample;

        par.exec.funname ={'TimeSelect','FilterBankCompute'};
        EEG_trials1 =run_trials(EEG_trials,par);

        par.TimeSelect.t1            = int4sub3; % in s from ZeroEvent time
        par.TimeSelect.t2            = int4sub4; % in s from ZeroEvent time

        itr1 = [int4sub1,int4sub2,int4sub3,int4sub4];
        itr2 = [int4sub1,int4sub2,int4sub3,int4sub4];

        EEG_trials2 = run_trials(EEG_trials,par);

        % concatenate Train Trials
        par.concatenateEEG.InField = signal_name;
        par.concatenateEEG.fsample = fsample;
        EEG_trials = concatenateEEG(EEG_trials1,EEG_trials2,par.concatenateEEG);


        % kfold-CrossValidation on the Train dataset
        kfoldSplit   = 10;
        labs    = [EEG_trials.trialType]'; %true labels
        cvp     = cvpartition(labs,'kfold',kfoldSplit,'Stratify',true);

        resQDA          = struct();
        resKNN          = struct();

        outQDA          = struct();
        outKNN          = struct();

        LabpredictQDA   = struct();
        LabpredictKNN   = struct();
        Label_train     = struct();

        for i=1:kfoldSplit
            indices = training(cvp,i);
            test = (indices == 0);
            train = ~test;
            EEG_train = EEG_trials(train);
            EEG_test = EEG_trials(test);

            Label_train(i).Iter = [EEG_test.trialType]';
            %% Step 2. perform pca
            % pca Dictionary evaluation on train
            par.pcaModel.InField          = signal_name;
            par.pcaModel.OutField         = signal_process;
            par.pcaModel.numComponents    = 0;
            par.pcaModel.perc             = 95;

            par.exec.funname ={'pcaModel'};
            [EEG_train,out] =run_trials(EEG_train,par);

            % pca Encode on test data
            par.pcaEncode.InField       = signal_name;
            par.pcaEncode.OutField      = signal_process;
            par.pcaEncode.Wpca          = out.pcaModel.Wpca;
            par.pcaEncode.mu            = out.pcaModel.mu;
            par.pcaEncode.explained     = out.pcaModel.explained;
            par.pcaEncode.numComponents = out.pcaModel.numComponents;
            par.exec.funname = {'pcaEncode'};

            EEG_train  = run_trials(EEG_train,par);
            EEG_test   = run_trials(EEG_test,par);
            
            TotalFeatures = size(EEG_test(1).(signal_process),2);

            %% Step 3. Model Classification on pca

            % qdaModel
            par.qdaModel                      = qdaModelParams;
            par.qdaModel.InField              = 'pca';
            par.qdaModel.numIterations        = 100;
            par.qdaModel.kfold                = 5;
            [~, outQDA(i).Iter]               = qdaModel(EEG_train,par.qdaModel);

            % predictQDA
            par.mdlPredict                  = mdlPredictParams;
            par.mdlPredict.InField          = 'pca';
            par.mdlPredict.OutField         = 'QDApred';
            par.mdlPredict.ProbField        = 'QDAProb';
            par.mdlPredict.mdl              = outQDA(i).Iter.mdl;
            [EEG_train,resQDA(i).train]     = mdlPredict(EEG_train,par.mdlPredict);
            [EEG_test,resQDA(i).test]       = mdlPredict(EEG_test,par.mdlPredict);

            LabpredictQDA(i).Iter = [EEG_test.QDApred]';

            % knnModel
            par.knnModel                      = knnModelParams;
            par.knnModel.InField              = 'pca';
            par.knnModel.numIterations        = 100;
            par.knnModel.kfold                = 5;
            [~, outKNN(i).Iter]             = knnModel(EEG_train,par.knnModel);

            % predictKNN
            par.mdlPredict                  = mdlPredictParams;
            par.mdlPredict.InField          = 'pca';
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
        params.createStructResult               = createStructResultParams();
        params.createStructResult.subj          = indsub;
        params.createStructResult.method        = signal_process;
        params.createStructResult.file          = 'Graz';
        params.createStructResult.train_name    = 'AT';
        params.createStructResult.train_tr1     = itr1;
        params.createStructResult.train_tr2     = itr2;
        params.createStructResult.test_name     = 'AT';
        params.createStructResult.test_ts1      = itr1;
        params.createStructResult.test_ts2      = itr2;
        params.createStructResult.m             = 0;
        params.createStructResult.class         = findclass(EEG_train,StartClass);
        params.createStructResult.irng          = par.irng;
        params.createStructResult.Filter        = par.FilterBankCompute.FilterBank;
        params.createStructResult.n_Features    = size(EEG_train(1).(signal_process),2);
        params.createStructResult.indMi         = 0;
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
        params.updateTab.dir        = 'D:\TrialBox_Results_excel\Graz_dataset';
        params.updateTab.name       = 'Graz_CrossVal_OnTrain_pipeline4';
        params.updateTab.sheetnames = 'QDA';

        updated_Result_tableAccQDA = updateTab(ResultQDA_Acc,params.updateTab);

        params.updateTab.sheetnames = 'KappaQDA';
        updated_Result_tableKappaQDA = updateTab(ResultQDA_Kappa,params.updateTab);

        params.updateTab.name     = 'Graz_CrossVal_OnTrain_class_pipeline4';
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
        params.updateTab.dir        = 'D:\TrialBox_Results_excel\Graz_dataset';
        params.updateTab.name       = 'Graz_CrossVal_OnTrain_pipeline4';
        params.updateTab.sheetnames = 'KNN';
        updated_Result_tableAccKNN = updateTab(ResultKNN_Acc,params.updateTab);

        params.updateTab.sheetnames = 'KappaKNN';
        updated_Result_tableKappaKNN = updateTab(ResultKNN_Kappa,params.updateTab);

        params.updateTab.name     = 'Graz_CrossVal_OnTrain_class_pipeline4';
        params.updateTab.sheetnames = 'KNN';
        updated_Resultclass_tableAccKNN = updateTab(ResultKNN_class_Acc,params.updateTab);
    end
end
