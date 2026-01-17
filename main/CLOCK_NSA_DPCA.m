% function CLOCK_NSA_DPCA
% PCA Clock, see
% Pezzulo, G., Donnarumma, F., Ferrari-Toniolo, S., Cisek, P., & Battaglia-Mayer, A. (2022). 
% Shared population-level dynamics in monkey premotor cortex during solo action, joint action and action observation. 
% Progress in Neurobiology, 210, 102214.
clear;
par.irng                        = 10;                  % for reproducibility
binWidth                        = 20;
kernSD                          = 30;
rng(par.irng);
%% Step 0. load data in raw format
nf = '~/DATA/SAPIENZA/PROCESSED/all_v5_translated_SK_SUA_t599_all_frontal.mat';                                 % combined sessions
% nf = '~/DATA/SAPIENZA/PROCESSED/all_v5_translated_SK_SUA_t599_all_frequencies20_frontal.mat';                 % alrealdy smoothed % 
% nf = '~/DATA/SAPIENZA/PROCESSED/all_v5_translated_SK_SUA_t599_all_means_frequencies20_cumulative_frontal';    % with mean trials per conditions as first trials 
data_trials                     = load(nf);
data_trials                     = data_trials.dat;
Treset                          = -199;
nTrials                         = length(data_trials);
for it=1:nTrials
    data_trials(it).trialTypeDir= data_trials(it).trialType;         % T1-T8
    data_trials(it).trialTypeCnd= data_trials(it).trialType2;        % joint | obs | solo 
    data_trials(it).trialType   = (data_trials(it).trialTypeCnd-1)*8+data_trials(it).trialTypeDir;
    data_trials(it).trialName   = SapienzaTypeString(data_trials(it).trialTypeCnd,data_trials(it).trialTypeDir);
    data_trials(it).trialNameAll= data_trials(it).trialName;
    data_trials(it).train       = true;
    data_trials(it).valid       = false;
    data_trials(it).test        = false;
    data_trials(it).timespikes  = (1:size(data_trials(it).spikes,2)) + Treset;
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
par.pcaModel.numComponents      = 8;
par.pcaModel.InField            = signal_process;
par.pcaModel.OutField           = signal_process;
%%%%%%%%%%%%%%%% nsa_pca 2-Stage Engine Churchland : kernel smooth + pca -> 
%%%%%%%%%%%%%%%% 'AverageWindow','GaussianSmoother','pcaCompute'
% see 
% Yu, B. M., Cunningham, J. P., Santhanam, G., Ryu, S., Shenoy, K. V., & Sahani, M. (2008).
% Gaussian-process factor analysis for low-dimensional single-trial analysis of neural population activity. 
% Advances in neural information processing systems, 21. 
% https://journals.physiology.org/doi/full/10.1152/jn.90941.2008
% Note: perform pca on trials averaged on conditions, similar to demixed pca. see->
% Kobak, D., Brendel, W., Constantinidis, C., Feierstein, C. E., Kepecs, A., Mainen, Z. F., ... & Machens, C. K. (2016).
% Demixed principal component analysis of neural population data. elife, 5, e10989.
% https://elifesciences.org/articles/10989
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

% Some plots
%% plot_trajectory2D
cmaps                           = linspecer(24);            % 8 directions x 3 conditions
cmaps                           = cmaps([17:24,1:8,9:16],:);% reset color as desired 
timeField                       = data_trials_class(1).(['time' par.pcaEncode.OutField]);
[~,i0]                          = min(abs(timeField));      % select   0  ms % TargetOn
istart                          = find(timeField>=-80,1);   % select -80  ms 
iend                            = find(timeField>=250,1);   % select 250  ms
par.plot_trajectory2D           = plot_trajectory2DParams;
% par.plot_trajectory2D.keep      = 1:8;                     % which conditions (e.g.  Solo)
% par.plot_trajectory2D.keep      = 9:16;                     % which conditions (e.g. obs)
par.plot_trajectory2D.keep      = 17:24;                     % which conditions (e.g. Join)
par.plot_trajectory2D.wd        = [3,4];                    % which components
par.plot_trajectory2D.axlab     = 'x';
par.plot_trajectory2D.cmaps     = cmaps;
par.plot_trajectory2D.cmapslight= lightCmaps(par.plot_trajectory2D.cmaps);
par.plot_trajectory2D.explained = out.pcaModel.explained;
par.plot_trajectory2D.InField   = par.pcaEncode.OutField;  % take pca projections to show
par.plot_trajectory2D.istart    = istart;
par.plot_trajectory2D.center    = i0;
par.plot_trajectory2D.iend      = iend;
hfg.plot_trajectory2D           = plot_trajectory2D(data_trials_class, par.plot_trajectory2D);

%% plot pcas dirs on conditions with plot_Each_mean - showing mean and confidence intervals
% remapTypes
par.remapTypes                  = remapTypesParams;
par.remapTypes.selection        = {  1:8  ,  9:16, 17:24 };
par.remapTypes.names            = {'Joint', 'Obs', 'Solo'};
data_trials_conds               = remapTypes(data_trials,par.remapTypes);
% plot_Each_mean
conditions_m                    = unique([data_trials_conds.trialType]);
nConditions_m                   = length(conditions_m);
cmaps                           = linspecer(nConditions_m); % 
cmaps                           = cmaps([2,1,3],:);         % color map for joint, obs and solo

% plot_EachDimVsTime
par.plot_EachDimVsTime              = plot_EachDimVsTimeParams;
par.plot_EachDimVsTime.cmaps        = cmaps;
par.plot_EachDimVsTime.cmapslight   = lightCmaps(par.plot_EachDimVsTime.cmaps);
par.plot_EachDimVsTime.InField      = par.pcaEncode.OutField;  % take pca projection to show;
par.plot_EachDimVsTime.ylabel       = '$pca$';
par.plot_EachDimVsTime.keep         = 1:nConditions_m;
par.plot_EachDimVsTime.nCols        = 4;
par.plot_EachDimVsTime.legplot      = 2;
% GRAPH SUBPLOT TITLES
explained                           = out.pcaModel.explained;
channels                            = 1:out.pcaModel.numComponents;
nChannels                           = length(channels); %% number of graphs
str                                 = cell(nChannels,1);
for ichannel=1:nChannels
    str{ichannel} = sprintf('$$\\tilde{\\mathbf x}_{%d,t}$$',ichannel);
    if ~isempty(explained)
        str{ichannel} = sprintf('%s (%2.1f)',str{ichannel},explained(ichannel));
    end
end
par.plot_EachDimVsTime.titles       = str;

% meanData
par.meanData                    = meanDataParams;
par.meanData.ifclass            = true;
par.meanData.trialTypes         = [data_trials_conds.trialType];
par.meanData.InField            = par.pcaEncode.OutField;
par.meanData.OutField           = par.pcaEncode.OutField;
par.meanData.P                  = 80;
par.meanData.opt                = [1,0,0];
data_trials_means               = meanData(data_trials_conds,par.meanData);%
%
hfg.pcas                        = plot_EachDimVsTime(data_trials_means,par.plot_EachDimVsTime);

%% plot_scatterGradient
data_trials_scatter             = data_trials;
InField                         = par.pcaEncode.OutField;
nTrials                         = length(data_trials_scatter);
sel                             = 9:16; nn='Obs';% Obs
names                           = unique({data_trials.trialNameAll});
goodtrials                      = false(nTrials,1);
for iTrial=1:nTrials
    time    = data_trials_scatter(iTrial).(['time' InField]);
    nTimes  = length(time);
    tType   = data_trials_scatter(iTrial).trialType;
    data_trials_scatter(iTrial).repTrialType =repmat(tType,1,nTimes);
    data_trials_scatter(iTrial).repTrialName=names(data_trials_scatter(iTrial).repTrialType);
    goodtrials(iTrial) = ismember(data_trials_scatter(iTrial).trialType,sel);
end
data_trials_scatter=data_trials_scatter(goodtrials);

% cmaps                           = linspecer(24);            % 8 directions x 3 conditions
% cmaps                           = cmaps([17:24,1:8,9:16],:);% reset color as desired 
% cmaps                           = cmaps(sel,:);
cmaps                           = linspecer(length(sel));

par.plot_scatterGradient                 = plot_scatterGradientParams();
par.plot_scatterGradient.InField         = InField;

par.plot_scatterGradient.InGradient      = ['time' InField];
par.plot_scatterGradient.lats            = [3,4];              % directions to be plot     
%par.plot_scatterGradient.lats            = [2,3,1];              % directions to be plot     
% start gradient color for each class
par.plot_scatterGradient.cmaps           = cmaps;
par.plot_scatterGradient.fine            = 1;
par.plot_scatterGradient.cmapslight      = lightCmaps(par.plot_scatterGradient.cmaps);
% 
hfg.plot_scatterGradient                 = plot_scatterGradient(data_trials_scatter,par.plot_scatterGradient);
