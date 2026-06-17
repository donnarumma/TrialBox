%% TEST_SPREAD_SPARSITY

clear; close all;
par.irng                        = 10;             % for reproducibility
rng(par.irng);
set(0, 'DefaultFigureVisible', 'off');

%% FRONTAL
nf = 'D:\TrialBox\Experimental\Experimental_Data\Total_Frontal_-200_600_BOOT50.mat';
nf_ch = 'PCA_-200_600';
INTERVAL = [-200,500];

chambName = 'FRONTAL';
numComponents = 12;
Def_numComponents = 1:numComponents;

% %% PARIETAL
% nf = 'D:\TrialBox\Experimental\Experimental_Data\Total_Parietal_-200_600_BOOT50.mat';
% nf_ch = 'PCA_-200_600';
% INTERVAL = [-200,500];
% chambName = 'PARIETAL';
% numComponents = 6;
% Def_numComponents = 1:numComponents;

% default params
defParams = getDefaultSPIKES_PLOTparams;
% merge params
params = defParams;

binWidth                        = params.binWidth;
kernSD                          = params.kernSD;
monkey                          = params.monkey;

name_int = num2str(INTERVAL(2));
tokens = regexp(nf, 'BOOT(\d+)', 'tokens');
bootName = str2double(tokens{1}{1});
name_dir = ['DPCA',filesep,'FINAL',filesep,nf_ch,filesep,'t0_t',name_int,filesep,chambName,filesep,'DYNAMICS',filesep,'BOOT',num2str(bootName), ...
    filesep,'pca',num2str(Def_numComponents(1)),'_',num2str(Def_numComponents(end)),'CURVES'];

fprintf('\n');
fprintf('================ SPIKES PLOT PARAMETERS ================\n');
fprintf('PCA number Components   : %d\n', numComponents);
fprintf('Smoothing Kernel (SD)   : %d ms\n', kernSD);
fprintf('Monkey names            : %s\n', strjoin(monkey, ', '));
fprintf('Save Directory          : %s\n', name_dir);
fprintf('========================================================\n\n');

%% Step 0. load data in raw format
fprintf('load data_trials\n');

trials_file                     = load(nf);
% data_trials                   = trials_file.data_trials;
% Treset                          = 0;
data_monkey                     = struct();
data_monkey_class               = struct();
try
    data_trials_EXTRACTED = trials_file.data_trials;
catch
    data_trials_EXTRACTED = trials_file.data_trials_BS;
end
for mnk=1:length(monkey)
    data_trials                 = data_trials_EXTRACTED;
    nTrials                     = length(data_trials);
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

    %% Step 1. prepare data
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

    %% Step 2. perform pca on trials averaged on conditions
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
    %%%%%%%%%%%%%%%% nsa_pca 2-Stage Engine Churchland : kernel smooth + pca ->
    %%%%%%%%%%%%%%%% 'AverageWindow','GaussianSmoother','pcaCompute'
    % see
    % Yu, B. M., Cunningham, J. P., Santhanam, G., Ryu, S., Shenoy, K. V., & Sahani, M. (2008).
    % Gaussian-process factor analysis for low-dimensional single-trial analysis of neural population activity.
    % Advances in neural information processing systems, 21.
    % https:\\journals.physiology.org\doi\full\10.1152\jn.90941.2008
    % Note: perform pca on trials averaged on conditions, similar to demixed pca. see->
    % Kobak, D., Brendel, W., Constantinidis, C., Feierstein, C. E., Kepecs, A., Mainen, Z. F., ... & Machens, C. K. (2016).
    % Demixed principal component analysis of neural population data. elife, 5, e10989.
    % https:\\elifesciences.org\articles\10989
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

    data_trials_class(1).explainedVar  = out.pcaModel.explained;
    %% Step 3. project trials on the pca dictionary found
    % pcaEncode
    % par.pcaEncode.Wpca              = out.pcaModel.Wpca;
    % par.pcaEncode.mu                = out.pcaModel.mu;
    % par.pcaEncode.explained         = out.pcaModel.explained;
    % par.pcaEncode.InField           = signal_process;
    % par.pcaEncode.OutField          = signal_process;
    % functions to be execute
    par.exec.funname                = {'AverageWindow','GaussianSmoother','pcaEncode'};
    data_trials                     = run_trials(data_trials,par);

    data_trials(1).explainedVar  = out.pcaModel.explained;
    data_trials(1).Monkey        = (monkey{mnk});
    %% Save data
    data_monkey.(monkey{mnk}) = data_trials;
    data_monkey_class.(monkey{mnk}) = data_trials_class;
end
pcaField = par.pcaEncode.OutField;

%% MND + Spread evaluation per set di componenti PCA
% Questa versione funziona quando ogni voce di Def_numComponentslist
% e' un vettore di componenti, per esempio:
% {[1 2], [3 4], [1 2 3], 1:numComponents}

data_distanceEval = data_monkey_class;

params.trimTrialsByTime.InField = pcaField;
params.trimTrialsByTime.tStart  = 0; %INTERVAL(1);
params.trimTrialsByTime.tStop   = INTERVAL(2);

numCond = 3;
numDirs = 8;
MDspread = struct();

% ESEMPI:
% Def_numComponentslist = {[1 2], [3 4], [1 2 3], [1 2 3 4]};
% Def_numComponentslist = num2cell(1:numComponents);

Def_numComponentslist = num2cell(1:4);
nSets = length(Def_numComponentslist);

for mnk = 1:length(monkey)

    mName = monkey{mnk};

    data_dist = data_distanceEval.(mName);
    data_dist = getSKconditiondirection(data_dist);
    data_dist = trimTrialsByTime(data_dist, params.trimTrialsByTime);

    timeFieldName = ['time' pcaField];
    timeField = data_dist(1).(timeFieldName);
    [~, i0] = min(abs(timeField));

    for iset = 1:nSets
        compSet = Def_numComponentslist{iset};
        compName = sprintf('DistComp%d_%d', compSet(1), compSet(end));
        for c = 1:numCond
            MV = NaN(1,numDirs);
            iMV = NaN(1,numDirs);
            tMV = NaN(1,numDirs);
            for d = 1:numDirs
                idx = ([data_dist.Condition] == c) & ([data_dist.Direction] == d);
                    if ~any(idx)
                        continue
                    end
                    v = data_dist(idx).(pcaField)(compSet,:);
                    v = v(~isnan(v));
                    [~,iMV(d)] = max(abs(v));
                    MV(d) = v(iMV(d));
                    tMV(d) = data_dist(1).(['time',(pcaField)])(iMV(d));
            end
            
            [vSorted, vDir] = sort(MV, 'ascend');
            dv = diff(vSorted);
            sp = sum(dv.^2);

            MDspread.(mName).val(iset).cond(c).cond     = c;
            MDspread.(mName).val(iset).cond(c).MV       = MV;
            MDspread.(mName).val(iset).cond(c).Tmd      = tMV;           
            MDspread.(mName).val(iset).cond(c).sorted   = vSorted;
            MDspread.(mName).val(iset).cond(c).orderDir = vDir;
            MDspread.(mName).val(iset).cond(c).diff     = dv;
            MDspread.(mName).val(iset).cond(c).spread   = sp;
        end
    end
end


data_distanceEval = data_monkey_class;

params.trimTrialsByTime.InField = pcaField;
params.trimTrialsByTime.tStart  = INTERVAL(1);
params.trimTrialsByTime.tStop   = INTERVAL(2);

Distance = struct();
Mparams = struct();

for iset = 1:nSets
    compSet = Def_numComponentslist{iset};
    compName = sprintf('DistComp%d_%d', compSet(1), compSet(end));

    for mnk = 1:length(monkey)
        monkeyName = monkey{mnk};

        data_dist = data_distanceEval.(monkeyName);
        data_dist = getSKconditiondirection(data_dist);
        data_dist = trimTrialsByTime(data_dist, params.trimTrialsByTime);

        timeField = data_dist(1).(['time' pcaField]);
        [~, i0] = min(abs(timeField));

        parDist = struct();
        parDist.InField = pcaField;
        parDist.wd      = compSet;
        parDist.center  = i0;

        distStruct  = compute_pca_distances_clean(data_dist, parDist);
        statsStruct = compute_pca_metrics_clean(distStruct);

        Distance.(compName).(monkeyName).trial = distStruct;
        Distance.(compName).(monkeyName).stats = statsStruct;
        Distance.(compName).(monkeyName).time  = timeField;
        Distance.(compName).(monkeyName).i0    = i0;
        Distance.(compName).(monkeyName).wd    = compSet;
    end

    S_stats = Distance.(compName).S.stats.condDir;
    K_stats = Distance.(compName).K.stats.condDir;

    % MND
    S_MaxDist = structToMatrix(S_stats, 'MND');
    K_MaxDist = structToMatrix(K_stats, 'MND');

    % MNDT
    S_timeMaxDist = structToMatrix(S_stats, 'MNDT');
    K_timeMaxDist = structToMatrix(K_stats, 'MNDT');

    Mparams.S.val(iset).MND     = S_MaxDist;
    Mparams.K.val(iset).MND     = K_MaxDist;
    Mparams.S.val(iset).MNDT    = S_timeMaxDist;
    Mparams.K.val(iset).MNDT    = K_timeMaxDist;
end

nMonkey = length(monkey);
nComp = nSets;
nDir = 8;
nOut = nMonkey * nComp * nDir;

parPlot = struct();
parPlot.plotPCAwithMDpoints.save_dir    = fullfile(name_dir,'SingleDir');
parPlot.plotPCAwithMDpoints.tStartPlot  = -100;
parPlot.plotPCAwithMDpoints.tEndPlot    = 600;
parPlot.plotPCAwithMDpoints.save        = true;

outPlot = struct();
timeCheckCell = cell(nOut, 1);

k = 0;

for mnk = 1:nMonkey
    monkeyName = monkey{mnk};

    for ncomp = 1:nComp
        for dirIdx = 1:nDir
            k = k + 1;

            outPlot.(monkeyName).comp(ncomp).dir(dirIdx) = ...
                plotPCAwithMDpoints(monkeyName, dirIdx, ncomp, ...
                data_monkey_class, MDspread, Mparams, pcaField, ...
                parPlot.plotPCAwithMDpoints);

            timeCheckCell{k} = outPlot.(monkeyName).comp(ncomp).dir(dirIdx).timeCheck;
        end
    end
end

allTimeCheck = vertcat(timeCheckCell{:});

parPlot.plotPCAwithMDpoints_8Dir.tStartPlot = -100;
parPlot.plotPCAwithMDpoints_8Dir.tEndPlot   = 600;
parPlot.plotPCAwithMDpoints_8Dir.save_dir   = fullfile(name_dir,'AllDir');
parPlot.plotPCAwithMDpoints_8Dir.name_fig   = NaN;
parPlot.plotPCAwithMDpoints_8Dir.save       = true;  

for mnk = 1:nMonkey
    monkeyName = monkey{mnk};
    for ncomp = 1:nComp
        plotPCAwithMDpoints_8Dir(monkeyName, ncomp, data_monkey_class, MDspread, Mparams, pcaField, parPlot.plotPCAwithMDpoints_8Dir);
    end
end

return
%% PLOT DYNAMICS 2 comp
data_dist = data_distanceEval.K;
data_dist = getSKconditiondirection(data_dist);
data_dist = trimTrialsByTime(data_dist, params.trimTrialsByTime);

pca2D = data_dist(1).dpca;
pc3 = pca2D(3,:);
pc4 = pca2D(4,:);
time = data_dist(1).timedpca;

pca2Dplot = [pc3; pc4];
par.compNames = {'PC3','PC4'};
par.mainTitle = 'Monkey K - PC3 PC4';
par.useAbsForPeak = true;
par.doAnimate = false;
par.saveFig = false;
out = plotPCA2Dtime(pca2Dplot, time, par);

%% 2) ANIMAZIONE SEMPLICE DELLA TRAIETTORIA
figure('Color','w','Name','PC3-PC4 animation');
axis equal;
grid on;
hold on;
xlabel('PC3');
ylabel('PC4');
title('Animazione traiettoria PC3-PC4');

xlim([min(pc3)-0.1*range(pc3), max(pc3)+0.1*range(pc3)]);
ylim([min(pc4)-0.1*range(pc4), max(pc4)+0.1*range(pc4)]);

hLine  = plot(nan, nan, 'b-', 'LineWidth', 1.5);
hPoint = plot(nan, nan, 'ro', 'MarkerFaceColor','r', 'MarkerSize',8);

for i = 1:length(time)
    set(hLine,  'XData', pc3(1:i), 'YData', pc4(1:i));
    set(hPoint, 'XData', pc3(i),   'YData', pc4(i));
    title(sprintf('Traiettoria PC3-PC4 - t = %d ms', time(i)));
    drawnow;
    pause(0.03);
end