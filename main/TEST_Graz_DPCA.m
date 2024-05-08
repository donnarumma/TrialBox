% TEST_Graz_DPCA.m
clear;

ifplot                          = false;

par.irng                        = 10;                  % for reproducibility
rng(par.irng);
addpath(genpath('D:\D_Ausilio EEG\EEG_FITTS\NSA_FITTS'))
addpath(genpath('D:\D_Ausilio EEG\NSA'))

binWidth                        = 30;
kernSD                          = 150;
% binWidth                        = 20;
% kernSD                          = 30;
P                               = 95;

indsub = 6;

%% Step 0. Load raw data and arrange trials
signal_name                     = 'eeg';

par.extractGraz.signal_name     = signal_name;
par.extractGraz.InField         = 'train';
[EEG_trials,fsample]            = extractGraz(indsub,par.extractGraz);

% Time interpolation and selection from Cue
par.TimeSelect               = TimeSelectParams;
par.TimeSelect.t1            = 0.0; % in s from ZeroEvent time
par.TimeSelect.t2            = 4.0; % in s from ZeroEvent time
par.TimeSelect.InField       = signal_name;
par.TimeSelect.OutField      = signal_name;
par.TimeSelect.dt            = 1;

par.exec.funname ={'TimeSelect'};
[EEG_trials, par.execinfo]=run_trials(EEG_trials,par);

%% Step 2. Smooth for dpca processing
signal_process                  = 'dpca';        % data processed are in 'y' field

% signalNormalize
par.signalNormalize             = signalNormalizeParams;
par.signalNormalize.InField     = signal_name;
par.signalNormalize.OutField    = signal_process;

% SmoothWindow (moving window smoothing)
par.SmoothWindow                = SmoothWindowParams;
par.SmoothWindow.InField        = signal_process;
par.SmoothWindow.OutField       = signal_process;
par.SmoothWindow.binWidth       = binWidth;
% removeInactives (0 mean channels removal)
par.removeInactive              = removeInactiveParams;
par.removeInactive.InField      = signal_process;
par.removeInactive.OutField     = signal_process;
% dataSplit
par.dataSplit                   = dataSplitParams;
par.dataSplit.TrainPercentage   = 80;
% function to be execute
par.exec.funname                = {'signalNormalize','SmoothWindow','removeInactive','dataSplit'};
[data_trials, outStep2]         = run_trials(EEG_trials,par);

%% Step 3. Perform pca on trials averaged on conditions  
% meanData
par.meanData                    = meanDataParams;
par.meanData.trialTypes         = [data_trials([data_trials.train]).trialType];
par.meanData.InField            = signal_process;
par.meanData.OutField           = signal_process;
% AverageWindow (average window selection with sqrt)
par.AverageWindow               = AverageWindowParams;
par.AverageWindow.useSqrt       = false;
par.AverageWindow.InField       = signal_process;
par.AverageWindow.OutField      = signal_process;
par.AverageWindow.binWidth      = binWidth;
% GaussianSmoother (kernel smoothing)
par.GaussianSmoother            = GaussianSmootherParams;
par.GaussianSmoother.InField    = signal_process;
par.GaussianSmoother.OutField   = signal_process;
par.GaussianSmoother.kernSD     = kernSD;       % standard deviation of Gaussian kernel, in msec
par.GaussianSmoother.stepSize   = binWidth;     % time between 2 consecutive datapoints, in msec
% pcaModel or pcaModel
par.pcaModel                    = pcaModelParams();
% par.pcaModel.numComponents    = 8;
par.pcaModel.numComponents      = 0;
par.pcaModel.perc               = 99;
par.pcaModel.InField            = signal_process;
par.pcaModel.OutField           = signal_process;
%%%%%%%%%%%%% nsa_pca 2-Stage Engine Churchland : kernel smooth + pca -> %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% 'AverageWindow','GaussianSmoother','pcaModel'            %%%%%%%%%%%%%%%%%%
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
[data_trials_class, out]        = run_trials(data_trials([data_trials.train]),par);

% pcaEncode
par.pcaEncode.Wpca              = out.pcaModel.Wpca;
par.pcaEncode.mu                = out.pcaModel.mu;
par.pcaEncode.explained         = out.pcaModel.explained;
par.pcaEncode.InField           = signal_process;
par.pcaEncode.OutField          = signal_process;
par.exec.funname                = {'pcaEncode'};
data_trials_class               = run_trials(data_trials_class,par);


%% Step 4. Project trials on the dictionary found 
% pcaEncode or Project
par.pcaEncode.Wpca              = out.pcaModel.Wpca;
par.pcaEncode.mu                = out.pcaModel.mu;
par.pcaEncode.explained         = out.pcaModel.explained;
par.pcaEncode.InField           = signal_process;
par.pcaEncode.OutField          = signal_process;
par.exec.funname                = {'AverageWindow','GaussianSmoother','pcaEncode'};
% par.exec.funname                = {'TimeSelect','AverageWindow','GaussianSmoother','pcaEncode'};
[data_trials, out2]             = run_trials(data_trials,par);

classes                         = unique([data_trials.trialType]);
nClasses                        = length(classes);
cmaps                           = linspecer(nClasses); % 
timeField                       = data_trials_class(1).(['time' par.pcaEncode.OutField]);

par.dpca_pCompute.it_start = par.TimeSelect.t1 + 0.25;
par.dpca_pCompute.it_stop = par.TimeSelect.t2 - 0.25;

valmin = nan(size(data_trials(1).dpca,1),1);
t_point = nan(size(data_trials(1).dpca,1),1);

for icomp = 1:size(data_trials(1).dpca,1)
[valmin(icomp),t_point(icomp)] = dpca_pCompute(icomp,data_trials,par.dpca_pCompute);
end

[pmin,dpca_comp] = min(valmin);
ctime = t_point(dpca_comp);

%% Plot Step
%% Plot 3D Trajectories  
[~,i0]                          = min(abs(timeField));      % select   0  ms % TargetOn
par.plot_trajectory3D           = plot_trajectory3DParams;
hfig                            = figure;
par.plot_trajectory3D.hfig      = hfig;
par.plot_trajectory3D.keep      = 1:nClasses;               % which conditions u-path, shortcut, deadend
par.plot_trajectory3D.wd        = 1:3;                      % which components
par.plot_trajectory3D.axlab     = 'x';
par.plot_trajectory3D.cmaps     = cmaps;
par.plot_trajectory3D.cmapslight= lightCmaps(par.plot_trajectory3D.cmaps);
par.plot_trajectory3D.explained = out.pcaModel.explained;
par.plot_trajectory3D.InField   = par.pcaEncode.OutField;  % take pca projections to show
par.plot_trajectory3D.istart    = 1;
par.plot_trajectory3D.center    = i0;
par.plot_trajectory3D.iend      = length(timeField);
hfg.trajectory3D                = plot_trajectory3D(data_trials_class,par.plot_trajectory3D);

%%
% 1:    Left Hand
% 2:    Right Hand
% 3:    Foot
% 4:    Tongue
% remapTypes
par.remapTypes                  = remapTypesParams;
par.remapTypes.selection        = { 1  , 2};
nSels                           = length(par.remapTypes.selection);
cmaps_LR                        = nan(nSels,3);
for isel=1:length(par.remapTypes.selection)
    cmaps_LR(isel,:)            = mean(cmaps(par.remapTypes.selection{isel},:),1);
end 
par.remapTypes.names            = {'Left Hand', 'Right Hand'};
data_trials_LR                  = remapTypes(data_trials,par.remapTypes);
%% meanData
par.meanData                    = meanDataParams;
par.meanData.ifclass            = true;
par.meanData.trialTypes         = [data_trials_LR.trialType];
par.meanData.InField            = par.pcaEncode.OutField;
par.meanData.OutField           = par.pcaEncode.OutField;
par.meanData.P                  = P;
par.meanData.opt                = [1,0,0];
data_trials_means_LR            = meanData(data_trials_LR,par.meanData);%
%%
% remapTypes
par.remapTypes                  = remapTypesParams;
par.remapTypes.selection        = {  3  ,  4 };
nSels                           = length(par.remapTypes.selection);
cmaps_FT                        = nan(nSels,3);
for isel=1:nSels
    cmaps_FT(isel,:)            = mean(cmaps(par.remapTypes.selection{isel},:),1);
end 
par.remapTypes.names            = {'Foot', 'Tongue'};
data_trials_FT                  = remapTypes(data_trials,par.remapTypes);
% meanData
par.meanData                    = meanDataParams;
par.meanData.ifclass            = true;
par.meanData.trialTypes         = [data_trials_FT.trialType];
par.meanData.InField            = par.pcaEncode.OutField;
par.meanData.OutField           = par.pcaEncode.OutField;
par.meanData.P                  = P;
par.meanData.opt                = [1,0,0];
data_trials_means_FT            = meanData(data_trials_FT,par.meanData);%


%%
% remapTypes
par.remapTypes                  = remapTypesParams;
par.remapTypes.selection        = { 1,2, 3 , 4 };
nSels                           = length(par.remapTypes.selection);
cmaps_LRFT                      = nan(nSels,3);
for isel=1:nSels
    cmaps_LRFT(isel,:)            = mean(cmaps(par.remapTypes.selection{isel},:),1);
end 
par.remapTypes.names            = {'Left Hand', 'Right Hand','Foot','Tongue'};
data_trials_LRFT                  = remapTypes(data_trials,par.remapTypes);
% meanData
par.meanData                    = meanDataParams;
par.meanData.ifclass            = true;
par.meanData.trialTypes         = [data_trials_LRFT.trialType];
par.meanData.InField            = par.pcaEncode.OutField;
par.meanData.OutField           = par.pcaEncode.OutField;
par.meanData.P                  = P;
par.meanData.opt                = [1,0,0];
data_trials_means_LRFT            = meanData(data_trials_LRFT,par.meanData);%


%% plot pcas dirs on conditions with plot_EachDimVsTime - showing mean and confidence intervals
% plot_Each_mean
classes                          = unique([data_trials_means_LR.trialType]);
nClasses                         = length(classes);
% plot_EachDimVsTime
par.plot_EachDimVsTime              = plot_EachDimVsTimeParams;
par.plot_EachDimVsTime.cmaps        = cmaps_LR;
par.plot_EachDimVsTime.cmapslight   = lightCmaps(par.plot_EachDimVsTime.cmaps);
par.plot_EachDimVsTime.InField      = par.pcaEncode.OutField;  % take pca projection to show;
par.plot_EachDimVsTime.ylabel       = '$pca$';
par.plot_EachDimVsTime.keep         = 1:nClasses;
par.plot_EachDimVsTime.nCols        = 3;
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

hfg.pcas                            = plot_EachDimVsTime(data_trials_means_LR,par.plot_EachDimVsTime);

%% plot pcas dirs on conditions with plot_EachDimVsTime - showing mean and confidence intervals
% plot_Each_mean
classes                          = unique([data_trials_means_FT.trialType]);
nClasses                         = length(classes);
% plot_EachDimVsTime
par.plot_EachDimVsTime              = plot_EachDimVsTimeParams;
par.plot_EachDimVsTime.cmaps        = cmaps_FT;
par.plot_EachDimVsTime.cmapslight   = lightCmaps(par.plot_EachDimVsTime.cmaps);
par.plot_EachDimVsTime.InField      = par.pcaEncode.OutField;  % take pca projection to show;
par.plot_EachDimVsTime.ylabel       = '$pca$';
par.plot_EachDimVsTime.keep         = 1:nClasses;
par.plot_EachDimVsTime.nCols        = 3;
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

hfg.pcas                       = plot_EachDimVsTime(data_trials_means_FT,par.plot_EachDimVsTime);

% 
%% plot pcas dirs on conditions with plot_EachDimVsTime - showing mean and confidence intervals
% plot_Each_mean
classes                          = unique([data_trials_means_LRFT.trialType]);
nClasses                         = length(classes);
% plot_EachDimVsTime
par.plot_EachDimVsTime              = plot_EachDimVsTimeParams;
par.plot_EachDimVsTime.cmaps        = cmaps_LRFT;
par.plot_EachDimVsTime.cmapslight   = lightCmaps(par.plot_EachDimVsTime.cmaps);
par.plot_EachDimVsTime.InField      = par.pcaEncode.OutField;  % take pca projection to show;
par.plot_EachDimVsTime.ylabel       = '$pca$';
par.plot_EachDimVsTime.keep         = 1:nClasses;
par.plot_EachDimVsTime.nCols        = 3;
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

hfg.pcas                       = plot_EachDimVsTime(data_trials_means_LRFT,par.plot_EachDimVsTime);


%% pSeparability
par.pSeparability                   = pSeparabilityParams;
par.pSeparability.InField           = par.pcaModel.OutField;
par.pSeparability.OutField          = 'comparisons';
par.pSeparability.exec              = true;
[pVals,pClasses]                    = pSeparability(data_trials,par.pSeparability);

%% pvalue plot per feature
par.plot_pValues                    = plot_pValuesParams;
par.plot_pValues.InField            = par.pSeparability.OutField;
par.plot_pValues.hfig               = figure('visible',ifplot);
par.plot_pValues.nRows              = 1;
par.plot_pValues.decisions          = [];
hfg.pvals                           = plot_pValues(pClasses,par.plot_pValues);
% sgtitle(hfg.pvals,[RatName ' ' daystr ' ' feedstr]);


%% save hfg
if ~ifplot
    save_dir                    = 'D:\TrialBox_Results_excel\NSA\Graz\Prova';
    save_dir                    = strcat(save_dir,datestr(datetime('now'), 'dd_mm_yy_HH_MM'),'\');
    par.hfigPrint               = hfigPrintParams();
    par.hfigPrint.pdf_file      = [save_dir mfilename '.pdf'];
    par.hfigPrint.save_dir      = save_dir; 
    hfigPrint(hfg,par.hfigPrint)
end

return