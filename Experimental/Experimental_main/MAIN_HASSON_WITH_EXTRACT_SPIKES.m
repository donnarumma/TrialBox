%% MAIN_HASSON_WITH_EXTRACT_SPIKES
%
% Computes Hasson-style directional coupling metrics between
% monkey S and monkey K neural manifolds starting from raw
% bootstrap population activity.
%
% Pipeline
%   1. Extract bootstrap neural populations
%   2. Build PCA neural manifolds
%   3. Construct S and K latent trajectories
%   4. Compute Hasson directional coupling matrices
%   5. Extract S->K and K->S coupling profiles
%   6. Generate temporal coupling plots
%   7. Generate directional difference plots
%
% Inputs
%   - Chamber selection (Frontal / Parietal)
%   - Bootstrap repetitions
%   - Neural analysis interval
%
% Outputs
%   - PCA manifold representations
%   - Hasson coupling matrices
%   - S->K and K->S coupling profiles
%   - K2S-S2K directional difference profiles
%   - PNG, PDF and MAT figures
%
% Author:
%   Mirco Frosolone
%
%
clear; close all;

origState = get(0,'DefaultFigureVisible');

cleanupObj = onCleanup(@() ...
    set(0,'DefaultFigureVisible',origState));

set(0,'DefaultFigureVisible','off');
%% Extraction parameters

Extractionparams.irng            = 10;
Extractionparams.interval        = [-200 600];
Extractionparams.chamber         = 'Frontal';   % 'Frontal' | 'Parietal'
Extractionparams.singleDir_save  = 'SingleSessionBOOT';
Extractionparams.dir_save        = 'TotalSessions';

rng(Extractionparams.irng);

M = 50;

%% =========================================================================
%  Spike Extraction
%  =========================================================================
data_trials = ExtractAllSPIKES_BOOTSTRAP(M,Extractionparams);

cfg = struct();
cfg.chamber = Extractionparams.chamber;
cfg.analysisInterval = [0 500];
cfg.hassonInterval = Extractionparams.interval;
cfg.bootRepetitions = M;

cfg.condNames = { ...
    'Act S - Obs K'
    'Obs S - Act K'
    'Joint'};

switch lower(cfg.chamber)

    case 'frontal'

        cfg.behaviourFile = ...
            'D:\NEW_MainScript_DCM\BEHAVIOUR_DATA\Behaviour_Frontal_TotalSessionsData_0_05s.mat';

        cfg.pcaFolderName = sprintf('PCA_%d_%d', ...
            Extractionparams.interval(1), ...
            Extractionparams.interval(2));

        cfg.numComponents = 12;

    case 'parietal'

        cfg.behaviourFile = ...
            'D:\NEW_MainScript_DCM\BEHAVIOUR_DATA\Behaviour_PARIETAL_TotalSessionsData_0_05s.mat';

        cfg.pcaFolderName = sprintf('PCA_%d_%d', ...
            Extractionparams.interval(1), ...
            Extractionparams.interval(2));

        cfg.numComponents = 6;

    otherwise

        error('Unknown chamber: %s',cfg.chamber);

end

condName        = cfg.condNames;
numComponents   = cfg.numComponents;
analysisWindow  = cfg.analysisInterval;
hassonWindow    = cfg.hassonInterval;

save_root = fullfile( ...
    pwd,...
    'HASSON',...
    'FINAL');
save_dir = fullfile( ...
    save_root,...
    cfg.pcaFolderName,...
    ['t0_t' num2str(cfg.analysisInterval(2))],...
    upper(cfg.chamber),...
    ['BOOT' num2str(cfg.bootRepetitions)]);

fprintf('Loading behavioural data...\n');
behavFile = load(cfg.behaviourFile);
if ~isfield(behavFile,'Data')
    error('Behaviour file does not contain Data structure.');
end

Data = behavFile.Data;

% default params
plotParams = getDefaultSPIKES_PLOTparams;

binWidth                        = plotParams.binWidth;
kernSD                          = plotParams.kernSD;
monkey                          = plotParams.monkey;
% numComponents                   = plotParams.numComponents;
wd2D                            = plotParams.wd2D;
wd3D                            = plotParams.wd3D;
keep                            = plotParams.keep;

fprintf('\n');
fprintf('================ HASSON PARAMETERS ================\n');
fprintf('Chamber                : %s\n', cfg.chamber);
fprintf('PCA Components         : %d\n', cfg.numComponents);
fprintf('Hasson Interval        : [%d %d] ms\n', ...
    cfg.hassonInterval(1), ...
    cfg.hassonInterval(2));
fprintf('Analysis Interval      : [%d %d] ms\n', ...
    cfg.analysisInterval(1), ...
    cfg.analysisInterval(2));
fprintf('Bootstrap repetitions  : %d\n', cfg.bootRepetitions);
fprintf('Behaviour Data         : %s\n', cfg.behaviourFile);
fprintf('===================================================\n\n');

fprintf('Using extracted spike dataset from current session...\n');

data_trials_EXTRACTED = data_trials;
%% =========================================================================
%  PCA Manifold Construction
%  =========================================================================
data_monkey       = struct();
data_monkey_class = struct();
for mnk=1:numel(monkey)
    data_trials                = data_trials_EXTRACTED;
    nTrials                    = length(data_trials);
    for it=1:nTrials
        % data_trials(it).trialTypeDir= data_trials(it).trialType;         % T1-T8
        % data_trials(it).trialTypeCnd= data_trials(it).trialType2;        % joint | obs | solo
        data_trials(it).trialType   = (data_trials(it).trialTypeCond-1)*8+data_trials(it).trialTypeDir;
        data_trials(it).trialName   = typeStringSK(data_trials(it).trialTypeCond,data_trials(it).trialTypeDir);
        data_trials(it).trialNameAll= data_trials(it).trialName;
        data_trials(it).train       = true;
        data_trials(it).valid       = false;
        data_trials(it).test        = false;
        data_trials(it).spikes      = data_trials(it).([monkey{mnk} '_Spikes']);
        data_trials(it).timespikes  = data_trials(it).timeSpikes(1,:);
    end

    %% Preprocess spike trains
    signal_name                     = 'spikes';       % data savead are in 'spikes' field
    signal_process                  = 'dpca';        % data processed are in 'y' field

    % SmoothWindow (moving window smoothing)
    par.SmoothWindow                = SmoothWindowParams;
    par.SmoothWindow.InField        = signal_name;
    par.SmoothWindow.OutField       = signal_process;
    par.SmoothWindow.binWidth       = binWidth;
    % removeInactives (0 mean channels removal)
    par.removeInactive              = removeInactiveParams;
    par.removeInactive.InField      = signal_process;
    par.removeInactive.OutField     = signal_process;

    % function to be execute
    par.exec.funname                = {'SmoothWindow','removeInactive'};
    data_trials                     = run_trials(data_trials,par);

    %% Compute PCA model
    % meanData
    par.meanData                    = meanDataParams;
    par.meanData.trialTypes         = [data_trials.trialType];
    par.meanData.InField            = signal_process;
    par.meanData.OutField           = signal_process;
    % AverageWindow (average window selection with sqrt)
    par.AverageWindow               = AverageWindowParams;
    par.AverageWindow.useSqrt       = true;
    par.AverageWindow.InField       = signal_process;
    par.AverageWindow.OutField      = signal_process;
    par.AverageWindow.binWidth      = binWidth;
    % GaussianSmoother (kernel smoothing)
    par.GaussianSmoother            = GaussianSmootherParams;
    par.GaussianSmoother.InField    = signal_process;
    par.GaussianSmoother.OutField   = signal_process;
    par.GaussianSmoother.kernSD     = kernSD;       % standard deviation of Gaussian kernel, in msec
    par.GaussianSmoother.stepSize   = binWidth;     % time between 2 consecutive datapoints, in msec
    % pcaModel
    par.pcaModel                    = pcaModelParams();
    par.pcaModel.numComponents      = numComponents;
    par.pcaModel.InField            = signal_process;
    par.pcaModel.OutField           = signal_process;
    par.exec.funname                = {'meanData','AverageWindow','GaussianSmoother','pcaModel'};
    [data_trials_class, out]        = run_trials(data_trials,par);
    % pcaEncode
    par.pcaEncode.Wpca              = out.pcaModel.Wpca;
    par.pcaEncode.mu                = out.pcaModel.mu;
    par.pcaEncode.explained         = out.pcaModel.explained;
    par.pcaEncode.InField           = signal_process;
    par.pcaEncode.OutField          = signal_process;
    par.exec.funname                = {'pcaEncode'};
    data_trials_class               = run_trials(data_trials_class,par);

    %% Project single trials
    par.exec.funname                = {'AverageWindow','GaussianSmoother','pcaEncode'};
    data_trials                     = run_trials(data_trials,par);

    data_trials(1).explainedVar     = out.pcaModel.explained;
    %% Save data
    data_monkey.(monkey{mnk})       = data_trials;
    data_monkey_class.(monkey{mnk}) = data_trials_class;
end
%% =========================================================================
%  Behavioural reference metrics
%  =========================================================================
behavMetrics = computeBehaviouralReferenceMetrics(Data);
%% =========================================================================
%  Hasson Directional Coupling Analysis
%  =========================================================================
%% Hasson configuration
cfg.componentGroups = { ...
    1 ...
    2 ...
    3 ...
    4 ...
    1:2 ...
    3:4 ...
    1:4};
nScanH              = numel(cfg.componentGroups);
cfg.hasson.block_size      = 10;%[10 10 10];
cfg.hasson.block_stride    = 0.25; %[0.25 0.4 0.6];
cfg.hasson.sub_block_size   = 3; %[3 3 3];
cfg.hasson.sub_block_stride = 0.5; %[0.5 0.5 0.5];
assert( ...
    numel(cfg.hasson.block_size) == numel(cfg.hasson.block_stride) && ...
    numel(cfg.hasson.block_size) == numel(cfg.hasson.sub_block_size) && ...
    numel(cfg.hasson.block_size) == numel(cfg.hasson.sub_block_stride), ...
    'Hasson parameter vectors must have matching lengths.');
%% Compute Hasson coupling matrices
hassonCfg = struct();
for nconfig = 1:length(cfg.hasson.block_size)
    for iC = 1:nScanH
        currComp = cfg.componentGroups{iC};
        compName = sprintf('Comp%d_%d', currComp(1),currComp(end));
        fprintf('Components: [%s]\n',num2str(currComp));
        fprintf('Output dir: %s\n',compName);
        fprintf('\n=== Hasson scan %d/%d: nComp = [%s] ===\n', iC, nScanH, num2str(currComp));

        hassonCfg.block_size          = cfg.hasson.block_size(nconfig);  % dimensione dei blocchi (lag) (int < trial_length)
        hassonCfg.block_stride        = cfg.hasson.block_stride(nconfig);

        hassonCfg.sub_block_size      = cfg.hasson.sub_block_size(nconfig);  %
        hassonCfg.sub_block_stride    = cfg.hasson.sub_block_stride(nconfig);  %

        hassonCfg.corr_obj            = "manifold";%"hat-hat"  % "manifold" | "hat-obs" | "hat-hat"
        hassonCfg.use_neural          = false;
        hassonCfg.trial_length        =  size(data_monkey.S(1).dpca,2); % lunghezza dei trial (int Fix)

        % trial lenght timespan in ms
        hassonCfg.t_start = hassonWindow(1);
        hassonCfg.t_end   = hassonWindow(2);

        A_struct=data_monkey.S;
        B_struct=data_monkey.K;

        for it=1:numel(A_struct)
            A_struct(it).Manifold       = A_struct(it).dpca(currComp,:);   % S
            A_struct(it).timeManifold   = A_struct(it).timedpca;   % S
            B_struct(it).Manifold       = B_struct(it).dpca(currComp,:);   % K
            B_struct(it).timeManifold   = B_struct(it).timedpca;   % K
        end
        % Y_subject=A_struct;
        % X_subject=B_struct;
        % define which subject
        hassonCfg.y_subject='s';
        hassonCfg.x_subject='k';
        %% Hasson Matrix linear plot
        t_start = analysisWindow(1);
        t_end   = analysisWindow(2);

        win  = 40;   %50
        step = 5;  %20

        t_vec = t_start:step:(t_end-win);
        nT = length(t_vec);
        % Condizioni = 3;
        % Direzioni  = 8;
        % Tempo      = nT;

        S2K_time = nan(3,8,nT);
        K2S_time = nan(3,8,nT);

        for it = 1:nT
            t_in  = t_vec(it);
            t_out = t_in + win;
            H_matrix = struct();
            nh = 1;
            for nCond = 1:3
                for nDir = 1:8
                    hassonCfg.dir  = nDir;
                    hassonCfg.cond = nCond;
                    L_L_matrix = compute_LL_matrix(A_struct,B_struct,hassonCfg);
                    % [H_matrix(nh).c_m, H_matrix(nh).stats] = ...
                    %     LL_detailPC(L_L_matrix,hassonCfg, t_in,t_out, t_in,t_out);
                    [H_matrix(nh).c_m, H_matrix(nh).stats] = ...
                        LL_detailPCadd(L_L_matrix,hassonCfg, t_in,t_out, t_in,t_out);
                    H_matrix(nh).Condition = nCond;
                    H_matrix(nh).Direction = nDir;
                    nh = nh + 1;
                end
            end
            extractPar.InField = 'mean_S2K';
            [c1,c2,c3] = extractHassonMatrix(H_matrix,extractPar);
            S2K_time(:,:,it) = [c1;c2;c3];
            extractPar.InField = 'mean_K2S';
            [c1,c2,c3] = extractHassonMatrix(H_matrix,extractPar);
            K2S_time(:,:,it) = [c1;c2;c3];
        end

        %% Plot S2K vs K2S
        plotHassonProfiles( ...
            S2K_time,...
            K2S_time,...
            t_vec,...
            condName,...
            behavMetrics.RT_mean,...
            behavMetrics.ExT_mean,...
            save_dir,...
            compName);


        %% Plot directional difference
        plotHassonDifference( ...
            S2K_time,...
            K2S_time,...
            t_vec,...
            condName,...
            behavMetrics.RT_mean,...
            behavMetrics.ExT_mean,...
            save_dir,...
            compName);
    end
end
close all