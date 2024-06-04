% TEST_Graz_1VS1_pipeline2.m

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

for c=3:4

    [tsub1,tsub3] = caseSelect(c);

    for indsub1 = 1:2
        for indsub2 = 1:9
            if c==1
                int4sub1 = tsub1(indsub1) - 0.75;
                int4sub2 = tsub1(indsub1) + 0.75;
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
                    int4sub3 = tsub3(indsub1) - 0.75;
                    int4sub4 = tsub3(indsub1) + 0.75;
                    if  int4sub3 <= 0.5
                        int4sub3 = 0.5;
                    elseif int4sub4 >= 2.5
                        int4sub4 = 2.5;
                    end
                else
                    int4sub3 = tsub3(indsub1) - 0.1*deltaT;
                    int4sub4 = tsub3(indsub1) + 0.1*deltaT;
                end
                signal_name                     = 'eeg';
                signal_process                  = 'CSP';

                %% Extract and Arrange TRAIN Data
                par.extractGraz.signal_name                  = signal_name;
                par.extractGraz.InField                      = 'train';
                [EEG_train,fsample] = extractGraz(indsub1,par.extractGraz);

                StartClass = unique([EEG_train.trialType]);
                % Time Interpolation and selection Trials [0.5;2.5] from CUE (Motor Imagery Interval)
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
                par.FilterBankCompute.attenuation   = 10;
                par.FilterBankCompute.FilterBank    = 'Nine';
                par.FilterBankCompute.fsample       = fsample;

                par.exec.funname ={'TimeSelect','FilterBankCompute'};
                EEG_train1 =run_trials(EEG_train,par);

                par.TimeSelect.t1            = int4sub3; % in s from ZeroEvent time
                par.TimeSelect.t2            = int4sub4; % in s from ZeroEvent time

                itr1 = [int4sub1,int4sub2,int4sub3,int4sub4];
                itr2 = [int4sub1,int4sub2,int4sub3,int4sub4];

                EEG_train2=run_trials(EEG_train,par);

                % concatenate Train Trials
                par.concatenateEEG.InField = signal_name;
                par.concatenateEEG.fsample = fsample;
                EEG_train = concatenateEEG(EEG_train1,EEG_train2,par.concatenateEEG);

                %% Extract and Arrange Test Data
                par.extractGraz.signal_name                  = signal_name;
                par.extractGraz.InField                      = 'test';
                [EEG_test,fsample] = extractGraz(indsub2,par.extractGraz);

                % Time Interpolation and selection Trials [0.5;2.5] from CUE (Motor Imagery Interval)
                par.TimeSelect               = TimeSelectParams;
                par.TimeSelect.t1            = int4sub1; % in s from ZeroEvent time
                par.TimeSelect.t2            = int4sub2; % in s from ZeroEvent time
                par.TimeSelect.InField       = signal_name;
                par.TimeSelect.OutField      = signal_name;
                par.TimeSelect.dt            = 1;


                % Filter Bank
                par.FilterBankCompute            = FilterBankComputeParams();
                par.FilterBankCompute.InField    = signal_name;
                par.FilterBankCompute.OutField   = signal_name;
                par.FilterBankCompute.attenuation = 10;
                par.FilterBankCompute.FilterBank = 'Nine';
                par.FilterBankCompute.fsample    = fsample;

                par.exec.funname ={'TimeSelect','FilterBankCompute'};
                EEG_test1 =run_trials(EEG_test,par);

                par.TimeSelect.t1            = int4sub3; % in s from ZeroEvent time
                par.TimeSelect.t2            = int4sub4; % in s from ZeroEvent time

                its1 = [int4sub1,int4sub2,int4sub3,int4sub4];
                its2 = [int4sub1,int4sub2,int4sub3,int4sub4];

                EEG_test2 =run_trials(EEG_test,par);

                EEG_test = concatenateEEG(EEG_test1,EEG_test2,par.concatenateEEG);
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
                par.miModel.k             = 4;

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
                params.createStructResult.subj       = indsub1;
                params.createStructResult.method     = 'CSP';
                params.createStructResult.file       = 'Graz';
                params.createStructResult.train_name = strcat('AT','_',string(indsub1));
                params.createStructResult.train_tr1  = itr1;
                params.createStructResult.train_tr2  = itr2;
                params.createStructResult.test_name  = strcat('AE','_',string(indsub2));
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
                params.updateTab.name               = 'Graz_1VS1_pipeline2';
                params.updateTab.sheetnames         = 'SVC';
                updated_Result_tableAccSVC          = updateTab(ResultSVC_Acc,params.updateTab);

                params.updateTab.sheetnames         = 'KappaSVC';
                updated_Result_tableKappaSVC        = updateTab(ResultSVC_kappa,params.updateTab);

                params.updateTab.name               = 'Graz_1VS1_class_pipeline2';
                params.updateTab.sheetnames         = 'SVC';
                updated_Resultclass_tableAccSVC     = updateTab(ResultSVC_class_Acc,params.updateTab);

               clear out resSVC
            end
        end
    end
end