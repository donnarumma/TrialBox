% TEST_Graz_CrossVal_OnTrain_pipeline3.m

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

% select the case: 1 2 3 4

for c=4

    [tsub1,tsub3] = caseSelect(c);

    for indsub = 1:9

        if c==1
            int4sub1 = tsub1(indsub) - 0.75;
            int4sub2 = tsub1(indsub) + 0.75;
            if  int4sub1 <= 0.5
                int4sub1 = 0.5;
            elseif int4sub2 >= 2.5
                int4sub2 = 2.5;
            end
        else
            int4sub1 = 0.5;
            int4sub2 = 2.5;
        end
        for deltaT = 1:2
            if c==3
                int4sub3 = tsub3(indsub) - 0.75;
                int4sub4 = tsub3(indsub) + 0.75;
                if  int4sub3 <= 0.5
                    int4sub3 = 0.5;
                elseif int4sub4 >= 2.5
                    int4sub4 = 2.5;
                end
            else
                int4sub3 = tsub3(indsub) - 0.1*deltaT;
                int4sub4 = tsub3(indsub) + 0.1*deltaT;
            end
        signal_name                     = 'eeg';
        signal_process1                  = 'CSP';
        signal_process2                  = 'pca';

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

        % kfold-CrossValidation on the Train dataset
        kfoldSplit   = 10;
        labs    = [EEG_trials.trialType]'; %true labels
        cvp     = cvpartition(labs,'kfold',kfoldSplit,'Stratify',true);

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

        for i=1:kfoldSplit
            indices = training(cvp,i);
            test = (indices == 0);
            train = ~test;
            EEG_train1 = EEG_trials1(train);
            EEG_test1 = EEG_trials1(test);

            EEG_train2 = EEG_trials2(train);
            EEG_test2 = EEG_trials2(test);

            Label_train(i).Iter = [EEG_test1.trialType]';
            %% Step 2. perform CSP
            par.cspModel                  = cspModelParams;
            par.cspModel.m                = 2;
            par.cspModel.InField          = signal_name;
            par.cspModel.OutField         = signal_process1;

            [~,out.cspModel1] = cspModel(EEG_train1,par.cspModel);
            % CSP Encode on train and test data
            par.cspEncode                  = cspEncodeParams;
            par.cspEncode.InField          = signal_name;
            par.cspEncode.OutField         = signal_process1;
            par.cspEncode.W                = out.cspModel1.W;

            par.exec.funname ={'cspEncode'};
            EEG_train1 =run_trials(EEG_train1,par);
            EEG_test1 =run_trials(EEG_test1,par);

            % Mutual Information
            par.miModel               = miModelParams;
            par.miModel.InField       = signal_process1;
            par.miModel.m             = par.cspModel.m;

            [~, out.miModel]=miModel(EEG_train1,par.miModel);

            par.miEncode               = miEncodeParams;
            par.miEncode.InField       = signal_process1;
            par.miEncode.OutField      = signal_process1;
            par.miEncode.IndMI         = out.miModel.IndMI;

            par.exec.funname ={'miEncode'};
            [EEG_train1, out]=run_trials(EEG_train1,par);
            V_train = out.miEncode.V;
            [EEG_test1, out]=run_trials(EEG_test1,par);
            V_test = out.miEncode.V;

            % pca Dictionary evaluation on train
            par.pcaModel.InField          = signal_name;
            par.pcaModel.OutField         = signal_process2;
            par.pcaModel.numComponents    = 0;
            par.pcaModel.perc             = 95;

            par.exec.funname ={'pcaModel'};
            [EEG_train2,out] =run_trials(EEG_train2,par);

            % pca Encode on test data
            par.pcaEncode.InField = signal_name;
            par.pcaEncode.OutField = signal_process2;
            par.pcaEncode.Wpca  = out.pcaSynergy.W;
            par.exec.funname ={'pcaEncode'};
            EEG_test2 =run_trials(EEG_test2,par);

            EEG_train = EEG_train1;
            EEG_test = EEG_test1;
            featuremix = 'CSP_pca';
            for iTr=1:length(EEG_train)
                EEG_train(iTr).(featuremix) = cat(2,EEG_train1(iTr).(signal_process1),EEG_train2(iTr).(signal_process2));
            end
            for iTr=1:length(EEG_test)
                EEG_test(iTr).(featuremix) = cat(2,EEG_test1(iTr).(signal_process1),EEG_test2(iTr).(signal_process2));
            end

            TotalFeatures = size(EEG_test(1).(signal_process),2);

            %% Step 3. Model Classification on CSP

            % qdaModel
            par.qdaModel                      = qdaModelParams;
            par.qdaModel.InField              = 'CSP_pca';
            par.qdaModel.numIterations        = 100;
            par.qdaModel.kfold                = 5;
            [~, outQDA(i).Iter]               = qdaModel(EEG_train,par.qdaModel);

            % predictQDA
            par.mdlPredict                  = mdlPredictParams;
            par.mdlPredict.InField          = 'CSP_pca';
            par.mdlPredict.OutField         = 'QDApred';
            par.mdlPredict.ProbField        = 'QDAProb';
            par.mdlPredict.mdl              = outQDA(i).Iter.mdl;
            [EEG_train,resQDA(i).train]     = mdlPredict(EEG_train,par.mdlPredict);
            [EEG_test,resQDA(i).test]       = mdlPredict(EEG_test,par.mdlPredict);

            LabpredictQDA(i).Iter = [EEG_test.QDApred]';

            % knnModel
            par.knnModel                      = knnModelParams;
            par.knnModel.InField              = 'CSP_pca';
            par.knnModel.numIterations        = 100;
            par.knnModel.kfold                = 5;
            [~, outKNN(i).Iter]             = knnModel(EEG_train,par.knnModel);

            % predictKNN
            par.mdlPredict                  = mdlPredictParams;
            par.mdlPredict.InField          = 'CSP_pca';
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
        % create Tab Result
        params.createStructResult               = createStructResultParams();
        params.createStructResult.subj          = indsub;
        params.createStructResult.method        = 'CSP_pca';
        params.createStructResult.file          = 'Graz';
        params.createStructResult.train_name    = 'AT';
        params.createStructResult.train_tr1     = itr1;
        params.createStructResult.train_tr2     = itr2;
        params.createStructResult.test_name     = 'AT';
        params.createStructResult.test_ts1      = itr1;
        params.createStructResult.test_ts2      = itr2;
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
        params.updateTab.dir        = 'D:\TrialBox_Results_excel\Graz_dataset';
        params.updateTab.name       = 'Graz_CrossVal_OnTrain_pipeline3';
        params.updateTab.sheetnames = 'QDA';

        updated_Result_tableAccQDA = updateTab(ResultQDA_Acc,params.updateTab);

        params.updateTab.sheetnames = 'KappaQDA';
        updated_Result_tableKappaQDA = updateTab(ResultQDA_Kappa,params.updateTab);

        params.updateTab.name     = 'Graz_CrossVal_OnTrain_class_pipeline3';
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

        % Update Tab Result
        params.updateTab.dir        = 'D:\TrialBox_Results_excel\Graz_dataset';
        params.updateTab.name       = 'Graz_CrossVal_OnTrain_pipeline3';
        params.updateTab.sheetnames = 'KNN';
        updated_Result_tableAccKNN = updateTab(ResultKNN_Acc,params.updateTab);

        params.updateTab.sheetnames = 'KappaKNN';
        updated_Result_tableKappaKNN = updateTab(ResultKNN_Kappa,params.updateTab);

        params.updateTab.name     = 'Graz_CrossVal_OnTrain_class_pipeline3';
        params.updateTab.sheetnames = 'KNN';
        updated_Resultclass_tableAccKNN = updateTab(ResultKNN_class_Acc,params.updateTab);

    end
end