% TEST_Sound_MultiEpoch_NoNSA.m


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

for indsub=1:9
    % indsub = 9;
    % for irng=1:9
    %     par.irng = irng;
    %     rng(par.irng);
        signal_name                     = 'eeg_sound';
        signal_process                  = 'CSP';

        %% Extract and Arrange Data
        par.extractSound.signal_name    = signal_name;
        par.extractSound.InField        = 'train';
        par.extractSound.it_end         = 2.5;
        par.extractSound.multiEpoch     = true;
        [EEG_trials,fsample]            = extractSound(indsub,par.extractSound);

        % remapTypes
        par.remapTypes           = remapTypesParams();
        par.remapTypes.selection = {1,2};

        StartClass = unique([EEG_trials.trialType]);
        % Time Interpolation and selection Trials
        par.TimeSelect               = TimeSelectParams;
        par.TimeSelect.t1            = 0.0; % in s from ZeroEvent time
        par.TimeSelect.t2            = 2.5; % in s from ZeroEvent time
        par.TimeSelect.InField       = signal_name;
        par.TimeSelect.OutField      = signal_name;

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

        %% epochCompute
        par.epochCompute                    = epochComputeParams();
        par.epochCompute.InField            = signal_name;
        par.epochCompute.OutField           = signal_name;
        par.epochCompute.fample             = fsample;
        par.epochCompute.t_epoch            = 0.5; % duration of a single intervals in s
        par.epochCompute.overlap_percent    = 0; % in percentage


        % par.exec.funname ={'remapTypes','TimeSelect','FilterBankCompute','eeg_overlap'};
        par.exec.funname ={'remapTypes','FilterBankCompute','epochCompute'};
        [EEG_trials,out.epochCompute] = run_trials(EEG_trials,par);
        itr1 = round(out.epochCompute.epochCompute.time_intervals(:,1),2);
        itr2 = round(out.epochCompute.epochCompute.time_intervals(:,2),2);

        % kfold-CrossValidation on the Train dataset
        kfoldSplit = 20;
        labs = [EEG_trials.trialType]'; %true labels
        cvp = cvpartition(labs,'kfold',kfoldSplit,'Stratify',true);

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
            EEG_train = repmat(EEG_trials(train),2,1);
            E_train = cat(1,EEG_trials(train).(signal_name));
            for iTrial=1:length(EEG_train)
                EEG_train(iTrial).trialId = iTrial;
                EEG_train(iTrial).(signal_name) = squeeze(E_train(iTrial,:,:,:));
            end

            % % Bootstrap
            % par.BootStrapData               = BootStrapDataParams;
            % par.BootStrapData.N             = 100;
            % par.BootStrapData.InField       = signal_name;
            % par.BootStrapData.OutField      = signal_name;
            % EEG_train                       = BootStrapData(EEG_train,par.BootStrapData);
            % EEG_train = EEG_train';


            EEG_test = repmat(EEG_trials(test),2,1);
            E_test = cat(1,EEG_trials(test).(signal_name));
            for iTrial=1:length(EEG_test)
                EEG_test(iTrial).trialId = iTrial;
                EEG_test(iTrial).(signal_name) = squeeze(E_test(iTrial,:,:,:));
            end
            Label_train(i).Iter = [EEG_test.trialType]';
            %% Step 2. perform CSP
            % CSP Dictionary evaluation on train
            par.cspModel                  = cspModelParams;
            par.cspModel.m                = 14;
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

            % Mutual Information
            par.miModel               = miModelParams;
            par.miModel.InField       = signal_process;
            par.miModel.m             = par.cspModel.m;
            % par.miModel.k             = 5;

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

            [EEG_train,resNBPW(i).train]     = predictNBPW(EEG_train,par.predictNBPW);

            % KappaValue on Train NBPW
            Cmatrxix_train = confusionmat([EEG_train.trialType]', pred_train);
            res_kappaNBPW(i).train.kappaValue  = kappaModel(Cmatrxix_train);

            % predictNBPW Test
            par.predictNBPW                  = predictNBPWParams;
            par.predictNBPW.InField          = 'CSP';
            par.predictNBPW.OutField         = 'NBPWpred';
            par.predictNBPW.labs_pred        = pred_test;

            [EEG_test,resNBPW(i).test]       = predictNBPW(EEG_test,par.predictNBPW);

            LabpredictNBPW(i).Iter = [EEG_test.NBPWpred]';

            % KappaValue on Test NBPW
            Cmatrxix_test = confusionmat([EEG_test.trialType]', pred_test);
            res_kappaNBPW(i).test.kappaValue  = kappaModel(Cmatrxix_test);
        end
        %% Accuracy and Kappa
        label_train = struct2cell(Label_train(:))'; % convert struct in mat
        label_train = cell2mat(label_train);

        predictQDA_train = struct2cell(LabpredictQDA(:))';
        predictQDA_train = cell2mat(predictQDA_train);
        predictKNN_train = struct2cell(LabpredictKNN(:))';
        predictKNN_train = cell2mat(predictKNN_train);
        predictNBPW_train = struct2cell(LabpredictNBPW(:))';
        predictNBPW_train = cell2mat(predictNBPW_train);

        AccuracyQDA = sum(predictQDA_train == label_train)/length(label_train)*100;
        accuracyQDA_class = accuracy4classes(label_train,predictQDA_train);

        AccuracyKNN = sum(predictKNN_train == label_train)/length(label_train)*100;
        accuracyKNN_class = accuracy4classes(label_train,predictKNN_train);

        AccuracyNBPW = sum(predictNBPW_train == label_train)/length(label_train)*100;
        accuracyNBPW_class = accuracy4classes(label_train,predictNBPW_train);

        CmatrxixQDA_train = confusionmat(label_train, predictQDA_train);
        kappaQDA = kappaModel(CmatrxixQDA_train);

        CmatrxixKNN_train = confusionmat(label_train, predictKNN_train);
        kappaKNN = kappaModel(CmatrxixKNN_train);

        CmatrxixNBPW_train = confusionmat(label_train, predictNBPW_train);
        kappaNBPW= kappaModel(CmatrxixNBPW_train);

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
        params.updateTab.name       = 'Sound_MultiEpoch_NoNSA';
        params.updateTab.sheetnames = 'QDA';

        updated_Result_tableAccQDA = updateTab(ResultQDA_Acc,params.updateTab);

        params.updateTab.sheetnames = 'KappaQDA';
        updated_Result_tableKappaQDA = updateTab(ResultQDA_Kappa,params.updateTab);

        params.updateTab.name     = 'Sound_MultiEpoch_class_NoNSA';
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
        params.updateTab.name       = 'Sound_MultiEpoch_NoNSA';
        params.updateTab.sheetnames = 'KNN';
        updated_Result_tableAccKNN = updateTab(ResultKNN_Acc,params.updateTab);

        params.updateTab.sheetnames = 'KappaKNN';
        updated_Result_tableKappaKNN = updateTab(ResultKNN_Kappa,params.updateTab);

        params.updateTab.name     = 'Sound_MultiEpoch_class_NoNSA';
        params.updateTab.sheetnames = 'KNN';
        updated_Resultclass_tableAccKNN = updateTab(ResultKNN_class_Acc,params.updateTab);

        % NBPW Acc
        resultNBPW.train.Accuracy = AccuracyNBPW;
        resultNBPW.train.Accuracy_class = accuracyNBPW_class;
        resultNBPW.test.Accuracy = NaN;
        resultNBPW.test.Accuracy_class = NaN;
        resultNBPW.train.kappaValue = kappaNBPW;
        resultNBPW.test.kappaValue = NaN;

        [ResultNBPW_Kappa,ResultNBPW_Acc,ResultNBPW_class_Acc] = createStructResult(resultNBPW,params.createStructResult);

        %% Update Tab Result
        params.updateTab.dir        = 'D:\TrialBox_Results_excel\Sound_dataset';
        params.updateTab.name       = 'Sound_MultiEpoch_NoNSA';
        params.updateTab.sheetnames = 'NBPW';
        updated_Result_tableAccNBPW = updateTab(ResultNBPW_Acc,params.updateTab);

        params.updateTab.sheetnames = 'KappaNBPW';
        updated_Result_tableKappaNBPW = updateTab(ResultNBPW_Kappa,params.updateTab);

        params.updateTab.name     = 'Sound_MultiEpoch_class_NoNSA';
        params.updateTab.sheetnames = 'NBPW';
        updated_Resultclass_tableAccNBPW = updateTab(ResultNBPW_class_Acc,params.updateTab);
    % end
end