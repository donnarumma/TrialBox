%% TEST_Graz_2A_SHALLOW_LEARN.m

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

numLatChan = 40;

for indsub=1:9
    
    attenuation = 10;

    if indsub==5
        attenuation = 20;
    end

    signal_name                     = 'eeg';
    signal_process                  = 'SHALLOW';

    %% Extract and Arrange TRAIN Data
    par.extractGraz.signal_name                  = signal_name;
    par.extractGraz.InField                      = 'train';
    [EEG_train,fsample] = extractGraz(indsub,par.extractGraz);

    StartClass = unique([EEG_train.trialType]);
    % Time Interpolation 
    par.TimeSelect               = TimeSelectParams;
    par.TimeSelect.t1            = -4; % in s from ZeroEvent time
    par.TimeSelect.t2            = 4; % in s from ZeroEvent time
    par.TimeSelect.InField       = signal_name;
    par.TimeSelect.OutField      = signal_name;
    par.TimeSelect.dt            = 1;

    EEG_train = TimeSelect(EEG_train,par.TimeSelect);

    % Time Interpolation and selection Trials [0.5;2.5] from CUE (Motor Imagery Interval)
    par.TimeSelect               = TimeSelectParams;
    par.TimeSelect.t1            = 0.5; % in s from ZeroEvent time
    par.TimeSelect.t2            = 2.5; % in s from ZeroEvent time
    par.TimeSelect.InField       = signal_name;
    par.TimeSelect.OutField      = signal_name;

    itr1 = par.TimeSelect.t1;
    itr2 = par.TimeSelect.t2;

    % Filter Bank
    par.FilterBankCompute            = FilterBankComputeParams();
    par.FilterBankCompute.InField    = signal_name;
    par.FilterBankCompute.OutField   = signal_name;
    par.FilterBankCompute.attenuation = attenuation;
    par.FilterBankCompute.FilterBank = 'One';
    par.FilterBankCompute.fsample    = fsample;

    par.exec.funname ={'TimeSelect','FilterBankCompute'};
    EEG_train =run_trials(EEG_train,par);

    %% Extract and Arrange Test Data
    par.extractGraz.signal_name                  = signal_name;
    par.extractGraz.InField                      = 'test';
    [EEG_test,fsample] = extractGraz(indsub,par.extractGraz);

    % Time Interpolation
    par.TimeSelect               = TimeSelectParams;
    par.TimeSelect.t1            = -4; % in s from ZeroEvent time
    par.TimeSelect.t2            = 4; % in s from ZeroEvent time
    par.TimeSelect.InField       = signal_name;
    par.TimeSelect.OutField      = signal_name;
    par.TimeSelect.dt            = 1;

    EEG_test = TimeSelect(EEG_test,par.TimeSelect);

    % Time Interpolation and selection Trials [0.5;2.5] from CUE (Motor Imagery Interval)
    par.TimeSelect               = TimeSelectParams;
    par.TimeSelect.t1            = 0.5; % in s from ZeroEvent time
    par.TimeSelect.t2            = 2.5; % in s from ZeroEvent time
    par.TimeSelect.InField       = signal_name;
    par.TimeSelect.OutField      = signal_name;

    its1 = par.TimeSelect.t1;
    its2 = par.TimeSelect.t2;

    % Filter Bank
    par.FilterBankCompute            = FilterBankComputeParams();
    par.FilterBankCompute.InField    = signal_name;
    par.FilterBankCompute.OutField   = signal_name;
    par.FilterBankCompute.attenuation = attenuation;
    par.FilterBankCompute.FilterBank = 'One';
    par.FilterBankCompute.fsample    = fsample;

    par.exec.funname ={'TimeSelect','FilterBankCompute'};
    EEG_test =run_trials(EEG_test,par);

    %--------------------------------------------------------------------------
    %% Step 2. perform SHALLOW
    ifplot                      = false;
    root_dir                    = 'D:\main_scriptNSA\Older_and_Proof\Main VAE\TEST_SHALLOW';
    model_file                  = false;
    ref                         = datetime('now','Format','yyyyMMddHHmmss');
    save_dir                    = sprintf('%s%s/',root_dir,ref);
    description                 = mfilename;
    out.ref                     = ref;

    numEpochs                   = 100;
    %--------------------------------------------------------------------------
    %% Learn Shallow Conv Net
    nTrials                         = length(EEG_train);
    seq                             = randperm(nTrials);

    par.shallow                     = shallowLearnParams();%vaeModelParams();

    %no validation
    par.shallow.log_valid         = [];
    par.shallow.log_train         = 1:288;

    % si validation
    % par.shallow.log_valid           = seq(1:80);
    % par.shallow.log_train           = seq(81:end);

    par.shallow.numLatentChannels   = numLatChan; % dimension of the manifold
    par.shallow.numEpochs           = numEpochs;
    par.shallow.InField             = signal_name;

    net                             = shallowLearn(EEG_train,par.shallow);
    %% Accuracy on Training
    labels_trainvalid       = categorical([EEG_train.trialType])';
    X3d                     = cat(3,EEG_train.(signal_name));                     % nCells x nTimes x nTrials
    X4d                     = reshape(X3d,size(X3d,1),size(X3d,2),1,size(X3d,3)); % nTrials x nChannels x 1 x nTimes
    XTrain                  = X4d(:,:,:,par.shallow.log_train);
    TTrain                  = labels_trainvalid(par.shallow.log_train);
    YTrain                  = classify(net,XTrain);
    accTrain                = mean(YTrain==TTrain);
    fprintf('Accuracy on Train: %g\n',accTrain)

    %% Accuracy on Validation
    if ~isempty(par.shallow.log_valid)
        XValid                  = X4d(:,:,:,par.shallow.log_valid);
        TValid                  = labels_trainvalid(par.shallow.log_valid);
        YValid                  = classify(net,XValid);
        accValid                = mean(YValid==TValid);
        fprintf('Accuracy on Validation: %g\n',accValid)
    end

    %% Accuracy on Test
    labels_test             = categorical([EEG_test.trialType])';
    X3dTest                 = cat(3,EEG_test.(signal_name));                     % nCells x nTimes x nTrials
    X4dTest                 = reshape(X3dTest,size(X3dTest,1),size(X3dTest,2),1,size(X3dTest,3)); % nTrials x nChannels x 1 x nTimes
    XTest                   = X4dTest;
    TTest                   = labels_test;
    YTest                   = classify(net,XTest);
    accTest                = mean(YTest==TTest);
    fprintf('Accuracy on Test: %g\n',accTest)


    AccShallow.acc_train = accTrain;
    AccShallow.acc_test = accTest;

    try AccShallow.acc_valid = accValid;
    catch
        AccShallow.acc_valid = 0;
    end

    % Save Result
    % create Tab Result
    params.createStructResultShallow.subj          = indsub;
    params.createStructResultShallow.method        = signal_process;
    params.createStructResultShallow.file          = 'Graz';
    params.createStructResultShallow.train_name    = 'AT';
    params.createStructResultShallow.train_tr1     = itr1;
    params.createStructResultShallow.train_tr2     = itr2;
    params.createStructResultShallow.test_name     = 'AE';
    params.createStructResultShallow.test_ts1      = its1;
    params.createStructResultShallow.test_ts2      = its2;
    params.createStructResultShallow.elem_train    = size(par.shallow.log_train,2);
    params.createStructResultShallow.elem_valid    = size(par.shallow.log_valid,2);
    params.createStructResultShallow.elem_test     = size(EEG_test,1);
    params.createStructResultShallow.m             = 0;
    params.createStructResultShallow.class         = findclass(EEG_train,StartClass);
    params.createStructResultShallow.irng          = par.irng;
    params.createStructResultShallow.Filter        = par.FilterBankCompute.FilterBank;
    params.createStructResultShallow.n_Features    = numLatChan;
    params.createStructResultShallow.indMi         = 0;
    params.createStructResultShallow.attenuation   = par.FilterBankCompute.attenuation;
    params.createStructResultShallow.TotalFeatures = numLatChan;
    params.createStructResultShallow.kfold         = 0;

    ResShallow = createStructResultSHALLOW(AccShallow,params.createStructResultShallow);

    % Update Tab Result
    params.updateTab.dir                = 'D:\TrialBox_Results_excel\Graz_dataset2A';
    params.updateTab.name               = 'Graz_2a_SHALLOW_LEARN';
    params.updateTab.sheetnames         = 'SHALLOW';
    updated_Result_tableAccQDA          = updateTab(ResShallow,params.updateTab);

end