function   par = TEST_SHORTCUT_NSA(ifplot,root_dir,AB,RatName,GO,evalinterval,irng)
% function par = TEST_SHORTCUT_NSA(ifplot)
% SHORTCUT DECISION, for methods see DPCA  
% Pezzulo, G., Donnarumma, F., Ferrari-Toniolo, S., Cisek, P., & Battaglia-Mayer, A. (2022). 
% Shared population-level dynamics in monkey premotor cortex during solo action, joint action and action observation. 
% Progress in Neurobiology, 210, 102214.
%%
if ~exist('irng','var')
    irng    = 10;
end
par.irng                        = irng;  % for reproducibility
rng(par.irng);
if ~exist('ifplot','var')
    ifplot  = true;
end
if ~exist('AB','var')
    AB      = true;    % trials AB or BA;
end
if ~exist('RatName','var')
    % RatName ='Rat067'; % subject
    RatName ='Rat066'; % subject
    % RatName ='Rat068'; % subject
end
if ~exist('root_dir','var')
    root_dir='';
end
par.irng                = 10;                  % for reproducibility
rng(par.irng);
%% Step 0. Load data in raw format 
runs=11;             % selected options (see params);
[S_cell,params]      = RatNSA_main(RatName,runs);

% params.iDecision     = iDecision;
Trials               = S_cell.Trials;
Data                 = S_cell.Data;
ratname              = Data.ratname;
iday                 = Data.day;
phase                = Data.phase;
kernSD               = params.kernSD;

% combine trial Decision 1 with opposite direction (decision 3)
Decisions  = [ 1,  3];  ABs = [AB,~AB];
% Decisions  = 1;         ABs = AB;

% params.AB            = AB;
gap_after            = params.data.gap_after;
gap_before           = params.data.gap_before;
questionable_cells   = params.data.load_questionable_cells;

endTR                = gap_after+gap_before; % trial duration in ms
if AB
    ABstr='AB';
else
    ABstr='BA';
end

description          = sprintf('%s_day%g_ph%g_D%g_%s_Q%g_b%ga%g_KRN%g',ratname,           ...
                                                                    iday,              ...
                                                                    phase,             ...
                                                                    Decisions(1),         ...
                                                                    ABstr,             ...
                                                                    questionable_cells,...
                                                                    gap_before,        ...
                                                                    gap_after,                                                                                      ...
                                                                    kernSD);
disp(description);

interval            = [];
labelDecision       = [];
condPerTrial        = [];
for iab=1:length(ABs)
    iDecision       = Decisions(iab);
    Traj_labels     = Trials.Traj_labels(Trials.ABlabels==ABs(iab),:);
    % Traj_edges      = Trials.Traj_edges(Trials.ABlabels==ABs(iab),:);
    % Traj_indexes    = Trials.Traj_indexes(Trials.ABlabels==ABs(iab),:);
    Traj_Decision   = Trials.Traj_Decision{iDecision}(Trials.ABlabels==ABs(iab),:);
    % T_Traj_Decision = Trials.T_Traj_Decision{iDecision}(Trials.ABlabels==ABs(iab),:);
    
    ABlabels        = Trials.ABlabels(Trials.ABlabels==ABs(iab),:);
    dinterval       = Data.T(Traj_Decision(:,[1,3]));
    interval        = [interval; dinterval];
    labelDecision   = [labelDecision; ones(size(dinterval,1),1)*iDecision];
    condPerTrial    = [condPerTrial;Traj_labels];
end

conditions           = unique(condPerTrial);
% NConditions          = length(conditions);
NeuronLabels         = Data.Spikes.label;
nNeurons             = length(NeuronLabels);
Tasklabels           = Trials.labels;

NTrialsOriginal     = size(interval,1);


% interval            = Data.T(Traj_Decision(:,[1,3]));
exclude_trials      = ~(abs(1000*diff(interval')-endTR)<exp(-10));
nTrials             = sum(~exclude_trials);
goodTrials          = find(~exclude_trials);
labelDecision       = labelDecision(goodTrials,:);
fprintf('Good Trials %g/%g\n',nTrials,NTrialsOriginal);
if nTrials<NTrialsOriginal
    Excluded=find(exclude_trials);
    for ie=1:length(Excluded)
        fprintf('Excluded: Trial %g\n',Excluded(ie));
    end
end
% create binned Spikes Struct
Spikes               = struct;
Spikes.label         = NeuronLabels;
Spikes.interval      = interval;
% extraction params
% freq_interval_in_bins=20;                           % freq on 20 bins
dt_in_s              = 1/1000;                      % 1ms bin
% Treset               = -gap_before;                 % trial zero time  

Spikes.t             = cell(nTrials,nNeurons);
Spikes.interval(exclude_trials,:)=[];
condPerTrial(exclude_trials)=[];
% loop on trials to select the chosen time window around the decision point [-gap_before,gap_after]
for iTrial=1:nTrials   
    for iNeuron=1:nNeurons
        cur_interval       =Spikes.interval(iTrial,:);
        spike_times        =Data.Spikes.t{iNeuron};
        sel_logical        =spike_times>cur_interval(1) & spike_times< cur_interval(2);            
        Spikes.t{iTrial,iNeuron}=spike_times(sel_logical)';
    end
end
% spikes
dat      =prepare_trials(Spikes,condPerTrial,[Spikes.label],dt_in_s); % bin spikes in dt_in_s s
dat      =setDatLabels(dat,Tasklabels);
dat      =setDatId(dat,goodTrials);
%% Step 1. Arrange to data_trials
data_trials                             = dat;
nTrials                                 = length(data_trials);
xfld                                    = 'time';
for it=1:nTrials
    data_trials(it).train               = true;
    data_trials(it).valid               = false;
    data_trials(it).xField              = xfld;
    data_trials(it).xFieldDesc          = 'ms';
    data_trials(it).test                = false;
    data_trials(it).labelDecision       = labelDecision(it);
    data_trials(it).([xfld 'spikes'])   = -gap_before:(dt_in_s*1000):gap_after; %ms
end
if AB
    feedstr =sprintf('feeder A->B');
else 
    feedstr =sprintf('feeder B->A');
end
daystr      = sprintf('Day%g',params.data.day);
titlestr    = [RatName ' ' daystr ' ' feedstr];

%% Step 2. Smooth for dpca processing. Compute dpca - CISEK pca
signal_name                     = 'spikes';       % data savead are in 'spikes' field
signal_process                  = 'dpca';        % data processed are in 'y' field
% SmoothWindow (moving window smoothing)
par.SmoothWindow                = SmoothWindowParams;
par.SmoothWindow.InField        = signal_name;
par.SmoothWindow.OutField       = signal_process;
par.SmoothWindow.binWidth       = params.binWidth;
% removeInactives (0 mean channels removal)
par.removeInactive              = removeInactiveParams;
par.removeInactive.InField      = signal_process;
par.removeInactive.OutField     = signal_process;
% function to be execute
par.exec.funname                = {'SmoothWindow','removeInactive'};
data_trials                     = run_trials(data_trials,par);

%% Step 3. Perform pca on trials averaged on conditions 
% meanData
par.meanData                    = meanDataParams;
% par.meanData.trialTypes         = [data_trials.trialType];
par.meanData.InField            = signal_process;
par.meanData.OutField           = signal_process;
% AverageWindow (average window selection with sqrt)
par.AverageWindow               = AverageWindowParams;
par.AverageWindow.useSqrt       = true;
par.AverageWindow.InField       = signal_process;
par.AverageWindow.OutField      = signal_process;
par.AverageWindow.binWidth      = params.binWidth;
% GaussianSmoother (kernel smoothing)
par.GaussianSmoother            = GaussianSmootherParams;
par.GaussianSmoother.InField    = signal_process;
par.GaussianSmoother.OutField   = signal_process;
par.GaussianSmoother.kernSD     = params.kernSD;       % standard deviation of Gaussian kernel, in msec
par.GaussianSmoother.stepSize   = params.binWidth;     % time between 2 consecutive datapoints, in msec
% pcaModel
par.pcaModel                    = pcaModelParams();
% par.pcaModel.numComponents    = 8;
par.pcaModel.numComponents      = 0;
par.pcaModel.perc               = 90;
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

% TimeSelect, compute pca in an interval
% par.TimeSelect                  = TimeSelectParams;
% par.TimeSelect.t1               = -1000;
% par.TimeSelect.t2               = 500;
% par.TimeSelect.InField          = signal_process;
% par.TimeSelect.OutField         = signal_process;
% par.exec.funname                = {'TimeSelect','meanData','AverageWindow','GaussianSmoother','pcaModel'};
par.exec.funname                = {'meanData','AverageWindow','GaussianSmoother','pcaModel'};
% par.exec.funname                = {'AverageWindow','GaussianSmoother','pcaModel'};
[data_trials_class, out]        = run_trials(data_trials,par);

% pcaProject
par.pcaEncode.Wpca              = out.pcaModel.Wpca;
par.pcaEncode.mu                = out.pcaModel.mu;
par.pcaEncode.explained         = out.pcaModel.explained;
par.pcaEncode.InField           = signal_process;
par.pcaEncode.OutField          = signal_process;
par.exec.funname                = {'pcaEncode'};
% par.exec.funname                = {'pcaEncode','meanData'};
data_trials_class               = run_trials(data_trials_class,par);


%% Step 4. Project trials on the dictionary found 
% pcaProject
% par.pcaProject.Wpca             = out.pcaModel.W;
% par.pcaProject.mu               = out.pcaModel.mu;
% par.pcaProject.explained        = out.pcaModel.explained;
% par.pcaProject.InField          = signal_process;
% par.pcaProject.OutField         = signal_process;
% % TimeSelect
% par.TimeSelect                  = TimeSelectParams;
% par.TimeSelect.t1               = -1000;
% par.TimeSelect.t2               = 1000;
% par.TimeSelect.InField          = signal_process;
% par.TimeSelect.OutField         = signal_process;
% functions to be execute
par.exec.funname                = {'AverageWindow','GaussianSmoother','pcaEncode'};
% par.exec.funname                = {'TimeSelect','AverageWindow','GaussianSmoother','pcaProject'};
data_trials                     = run_trials(data_trials,par);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Some plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmaps                           = linspecer(length(unique([data_trials_class.trialType]))); % u-path, shortcut, deadend
timeField                       = data_trials_class(1).(['time' par.pcaEncode.OutField]);
[~,i0]                          = min(abs(timeField));      % select   0  ms % TargetOn

% separate GO and RETURN
% Decision Point Decisions(1) trials in the forward direction
% (e.g. first decision in A->B direction)
data_trials_GO                  = data_trials([data_trials.labelDecision]==Decisions(1));
if length(Decisions)>1
    % Decision point Decisions(2) is the in reverse direction 
    % (e.g. third decision in B->A direction
    data_trials_RETURN          = data_trials([data_trials.labelDecision]==Decisions(2));
end
if GO
    data_trials                 = data_trials_GO;
else
    data_trials                 = data_trials_RETURN;
end
keep                            = unique([data_trials.trialType]);


%% Bootstrap
par.BootStrapData               = BootStrapDataParams;
par.BootStrapData.TimeLagged    = 1;
par.BootStrapData.N             = 100;
par.BootStrapData.InField       = signal_process;
par.BootStrapData.OutField      = signal_process;
bootdata_trials                 = BootStrapData(data_trials,par.BootStrapData);

%% plot 2D Trajectories
% i500                            = find(timeField>=-500,1);   % select -250  ms 
% plot_trajectory2D
par.plot_trajectory2D           = plot_trajectory2DParams;

par.plot_trajectory2D.keep      = keep;                      % which conditions u-path, shortcut, deadend
par.plot_trajectory2D.wd        = [1,2];                    % which components
par.plot_trajectory2D.axlab     = 'x';
par.plot_trajectory2D.cmaps     = cmaps;
par.plot_trajectory2D.cmapslight= lightCmaps(par.plot_trajectory2D.cmaps);
par.plot_trajectory2D.explained = out.pcaModel.explained;
par.plot_trajectory2D.InField   = par.pcaEncode.OutField;  % take pca projections to show
par.plot_trajectory2D.istart    = 1;
par.plot_trajectory2D.center    = i0;
par.plot_trajectory2D.iend      = length(timeField);

% par.plot_trajectory2D.hfig      = figure('visible',ifplot);
% hfg.trajectory2D_class          = plot_trajectory2D(data_trials_class, par.plot_trajectory2D);
% par.plot_trajectory2D.hfig      = figure('visible',ifplot);
% hfg.trajectory2D_meanboot       = plot_trajectory2D(meanData(bootdata_trials,par.meanData), par.plot_trajectory2D);
par.plot_trajectory2D.hfig      = figure('visible',ifplot);
hfg.trajectory2D_mean           = plot_trajectory2D(meanData(data_trials,par.meanData), par.plot_trajectory2D);
title(titlestr);

%% Plot 3D Trajectories  
par.plot_trajectory3D           = plot_trajectory3DParams;
par.plot_trajectory3D.hfig      = figure('visible',ifplot);
par.plot_trajectory3D.keep      = keep;                      % which conditions u-path, shortcut, deadend
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
title(titlestr);
% view([12 25]);


%% plot pcas dirs on conditions with plot_Each_mean - showing mean and confidence intervals
% conditions                      = unique([data_trials_class.trialType]);
% nConditions                     = length(conditions);
% cmaps                           = linspecer(nConditions); % 
% % plot_Each_mean
% par.plot_Each_mean              = nsa_plot_params();
% par.plot_Each_mean.opt          = [1,0,0];                  % use bootstrap
% par.plot_Each_mean.InField      = par.pcaEncode.OutField;  % take pca projection to show
% par.plot_Each_mean.keep         = conditions;
% par.plot_Each_mean.cmaps        = cmaps;
% par.plot_Each_mean.cmapslight   = lightCmaps(par.plot_Each_mean.cmaps);
% par.plot_Each_mean.legplot      = 2;
% par.plot_Each_mean.P            = 68;%68;
% par.plot_Each_mean.explained    = out.pcaModel.explained;
% par.plot_Each_mean.nCols        = min(4,out.pcaModel.numComponents);
% hfg.pcas                        = plot_Each_mean(data_trials, par.plot_Each_mean);

%% plot_Ellipsoid2D
par.plot_ellipsoid2D            = plot_ellipsoid2DParams;
par.plot_ellipsoid2D.hfig       = figure('visible',ifplot);
par.plot_ellipsoid2D.InField    = par.pcaEncode.OutField;
par.plot_ellipsoid2D.P          = 68;
par.plot_ellipsoid2D.opt        = [0,1,0];                  % bootstrap
par.plot_ellipsoid2D.wd         = [1,2];                    % which components

selected                        = keep;%1:3;
T                               = length(timeField);
dt                              = 5;
indsT                           = 1:dt:T;
par.plot_ellipsoid2D.explained  = out.pcaModel.explained;
try
    hfg.ellipsoid2D                 = plot_ellipsoid2D(data_trials,selected,indsT,cmaps,par.plot_ellipsoid2D);
catch
    hfg.ellipsoid2DBoot             = plot_ellipsoid2D(bootdata_trials,selected,indsT,cmaps,par.plot_ellipsoid2D);
end
nCols=2;
%% plot pcas dirs on conditions with plot_Each_mean - showing mean and confidence intervals
par.plot_EachDimVsTime              = plot_EachDimVsTimeParams;
par.plot_EachDimVsTime.hfig         = figure('visible',ifplot);
par.plot_EachDimVsTime.cmaps        = cmaps;
par.plot_EachDimVsTime.cmapslight   = lightCmaps(par.plot_EachDimVsTime.cmaps);
par.plot_EachDimVsTime.legplot      = 2;
par.plot_EachDimVsTime.keep         = 1:3;
par.plot_EachDimVsTime.InField      = signal_process;
par.plot_EachDimVsTime.ylabel       = '$pca$';
par.plot_EachDimVsTime.nCols        = nCols;
% GRAPH SUBPLOT TITLES
explained       = out.pcaModel.explained;
channels        = 1:out.pcaModel.numComponents;
nChannels       = length(channels); %% number of graphs
toleg           = cell(nChannels,1);
for ichannel=1:nChannels
    toleg{ichannel} = sprintf('$$\\tilde{\\mathbf x}_{%d,t}$$',ichannel);
    if ~isempty(explained)
        toleg{ichannel} = sprintf('%s (%2.1f)',toleg{ichannel},explained(ichannel));
    end
end
par.plot_EachDimVsTime.titles       = toleg;

% mean on pca
par.meanData                        = meanDataParams;
par.meanData.trialTypes             = [data_trials.trialType];
par.meanData.ifclass                = true;
par.meanData.P                      = 68;
par.meanData.opt                    = [1,0,0];
par.meanData.InField                = signal_process;
par.meanData.OutField               = signal_process;
par.meanData.SE                     = 0;
hfg.pcasVsTime                      = plot_EachDimVsTime(meanData(data_trials,par.meanData),par.plot_EachDimVsTime);
%%

%% pSeparability
par.pSeparability                   = pSeparabilityParams;
par.pSeparability.InField           = signal_process;
par.pSeparability.OutField          = 'comparisons';
par.pSeparability.exec              = true;
pdata_trials                        = bootdata_trials; % data_trials
[pVals,pClasses]                    = pSeparability(pdata_trials,par.pSeparability);

%% pvalue plot per feature
par.plot_pValues                    = plot_pValuesParams;
par.plot_pValues.InField            = par.pSeparability.OutField;
par.plot_pValues.hfig               = figure('visible',ifplot);
par.plot_pValues.exec               = 1;
par.plot_pValues.xfld               = 'time';
par.plot_pValues.dt                 = 50;
par.plot_pValues.nCols              = nCols;
par.plot_pValues.explained          = explained;
hfg.pvals   = plot_pValues(pClasses,par.plot_pValues);
sgtitle(hfg.pvals,titlestr);

%% BayesianInferenceClassificator, cumulated pca components
% time select: [-500,0] window
% evalinterval                    = [-1000,0];
evaltime                        = evalinterval(end);
par.TimeSelect                  = TimeSelectParams;
par.TimeSelect.t1               = evalinterval(1);%-1000;%-500;
par.TimeSelect.t2               = evalinterval(end);%0;
par.TimeSelect.InField          = signal_process;
par.TimeSelect.OutField         = signal_process;
% BayesianInferenceClassification
par.BayesianInferenceClassification            = BayesianInferenceClassificationParams();
par.BayesianInferenceClassification.InField    = signal_process;
par.BayesianInferenceClassification.exec       = true;
par.BayesianInferenceClassification.channelSets= {1:out.pcaModel.numComponents};        % all pcas together
% par.BayesianInferenceClassification.channelSets=num2cell(1:out.pcaModel.numComponents); % all pcas separated
% par.BayesianInferenceClassification.channelSets={1,2,1:3,4:6};                            % some pca combinations

%% train classifier for each bin and save in train field
data_to_class                                  = TimeSelect(bootdata_trials,par.TimeSelect);
par.BayesianInferenceClassification.train      = [data_to_class.train];
[~,out.BayesianInferenceClassification]        = BayesianInferenceClassification(data_to_class,par.BayesianInferenceClassification);
% par.BayesianInferenceClassification.train      = res.timemdl1;% 0;
par.BayesianInferenceClassification.train      = 0;
for iset=1:length(par.BayesianInferenceClassification.channelSets)
    tfld = ['timemdl' num2str(iset)];   
    par.BayesianInferenceClassification.(tfld)= out.BayesianInferenceClassification.(tfld);% 0;
end

%% classify original trials
% retrain if necessary
% par.BayesianInferenceClassification.train     = [data_trials.train];
data_trials_prob_raw                           = BayesianInferenceClassification(TimeSelect(data_trials,par.TimeSelect),par.BayesianInferenceClassification);
data_trials_prob                               = data_trials_prob_raw;
%% classify bootstrapped trials
% retrain if necessary
% par.BayesianInferenceClassification.train      = [data_to_class.train];
bootdata_trials_prob                           = BayesianInferenceClassification(data_to_class,par.BayesianInferenceClassification);
% data_trials_prob=bootdata_trials_prob;

%% some plots

%% plot classification probabilites in time
% plot_EachDimVsTime
par.plot_EachDimVsTime              = plot_EachDimVsTimeParams;
par.plot_EachDimVsTime.hfig         = figure('visible',ifplot);
par.plot_EachDimVsTime.cmaps        = cmaps;
par.plot_EachDimVsTime.cmapslight   = lightCmaps(par.plot_EachDimVsTime.cmaps);
par.plot_EachDimVsTime.InField      = 'prob';
par.plot_EachDimVsTime.ylabel       = '$P(C|o_{t})$';
par.plot_EachDimVsTime.keep         = 1:3;
par.plot_EachDimVsTime.nCols        = 1;
par.plot_EachDimVsTime.YLIM         = [0-0.01,1+0.01];
par.plot_EachDimVsTime.legplot      = 2;

% GRAPH SUBPLOT TITLES
explained    =out.pcaModel.explained;
channelSets  =par.BayesianInferenceClassification.channelSets;
nSets        =length(channelSets); %% number of graphs
toleg=cell(nSets,1);
for iSet=1:nSets
    channels     =channelSets{iSet};
    nChannels    =length(channels);
    toleg{iSet}='';
    for ichannel=1:nChannels
        toleg{iSet} = sprintf('%s$${\\mathbf x}_{%d,t}$$',toleg{iSet},channels(ichannel));
        if ichannel<nChannels
             toleg{iSet} = sprintf('%s,',toleg{iSet});
        end 
    end
    if ~isempty(explained)
        toleg{iSet} = sprintf('%s (%2.1f)',toleg{iSet},sum(explained(channelSets{iSet})));
    end
end
par.plot_EachDimVsTime.titles=toleg;

par.meanData                    = meanDataParams;
par.meanData.ifclass            = true;
par.meanData.trialTypes         = [data_trials_prob.trialType];
par.meanData.InField            = 'prob';
par.meanData.OutField           = 'prob';
par.meanData.P                  = 68;
par.meanData.opt                = [0,0,1];

hfg.CumPcaClassificationInTime  = plot_EachDimVsTime(meanData(data_trials_prob,par.meanData),par.plot_EachDimVsTime);
sgtitle(titlestr);

%% accuracy in trials
nBins                               = 8;
nTrials                             = length(data_trials_prob_raw);
DT                                  = nTrials-nBins+1;
LastWindow                          = nTrials-DT+1;
% nWindows                            = length(1:LastWindow);
% timeTrials                          = nan(DT,nWindows);
% evaltime                            = 0; % seconds
try
    exp_ind                         = getExpIndex(RatName,AB);
catch
    exp_ind=1;
end
save_dir                    = [root_dir description filesep];
accuracy_filename           = [root_dir mfilename '_accuracy'];
class_filename              = [root_dir mfilename '_class_accuracy'];
try 
    fprintf('Try loading %s\n',accuracy_filename);
    accuracy_file   = load(accuracy_filename);
    accuracy_window = accuracy_file.accuracy_window;

    fprintf('Try loading %s\n',class_filename);
    class_file      = load(class_filename);
    class_window    = class_file.class_window;
catch
    fprintf('WARNING: Failed loading %s\n',accuracy_filename);
    fprintf('WARNING: Failed loading %s\n',class_filename);
    pause;
end
% accuracy_window                     = nan(1,nBins);
for iClass=1:length(unique([data_trials_prob.trialType]))
    class_window{iClass}(exp_ind,:)=nan(1,nBins);
end
for iWindow=1:LastWindow
    idx_sel=iWindow:(iWindow+DT-1);
    SuccessField                        = 'success';
    data_trials_for_accuracy            = data_trials_prob_raw(idx_sel);
    
    % meanData -> get all trial means
    pms.meanData                        = meanDataParams;
    % pms.meanData.exec                   = [];
    pms.meanData.opt                    = [0,0,1]; % mean
    pms.meanData.trialTypes             = true(size([data_trials_for_accuracy.trialType]));
    pms.meanData.InField                = SuccessField;
    pms.meanData.OutField               = 'accuracy';
    acc_trials                          = meanData(data_trials_for_accuracy,pms.meanData);
    
    timextkl                            = acc_trials.([pms.meanData.xfld pms.meanData.OutField]);
    [~,indextime]                       = min(abs(timextkl-evaltime));
    accuracy_window(exp_ind,iWindow)    = acc_trials(1).(pms.meanData.OutField)(indextime);
    
    % for each class
    pms.meanData                        = meanDataParams;
    % pms.meanData.exec                   = [];
    pms.meanData.trialTypes             = [data_trials_for_accuracy.trialType];
    pms.meanData.opt                    = [0,0,1]; % mean
    pms.meanData.InField                = SuccessField;
    pms.meanData.OutField               = 'accuracy';
    class_acc_trials                    = meanData(data_trials_for_accuracy,pms.meanData);
    for idxclass=1:length(class_acc_trials)
        iClass = class_acc_trials(idxclass).trialType;
        class_window{iClass}(exp_ind,iWindow)=class_acc_trials(idxclass).(pms.meanData.OutField)(indextime);
    end
    
    % accuracy_window(iWindow)            = class_acc_trials(2).(pms.meanData.OutField)(indextime);
end


save(accuracy_filename,'accuracy_window');
save(class_filename,'class_window');

hfg.accuracyWindows=figure('visible',ifplot);hold on; grid on; box on;
plot(accuracy_window(exp_ind,:),'DisplayName','mean','linewidth',4,'Color',0.8*ones(1,3))
for idxclass=1:length(class_acc_trials)
    iClass = class_acc_trials(idxclass).trialType;
    plot(class_window{iClass}(exp_ind,:),'DisplayName',class_acc_trials(idxclass).trialName, 'linewidth',4,'color',cmaps(iClass,:));
end
title(titlestr)
ylabel(sprintf('Accuracy at t=%g ms',evaltime))
xlabel('bin');
legend;

%% plot hfg
if ~ifplot
    % save_dir                    = ['~/TESTS/DARTMOUTH/SHORTCUT/NSA/' description '/'];
    par.hfigPrint               = hfigPrintParams();
    par.hfigPrint.pdf_file      = [root_dir mfilename '_' description '.pdf'];
    par.hfigPrint.save_dir      = save_dir; 
    hfigPrint(hfg,par.hfigPrint)
end

return
%%
figure;plot(acc_trials(1).timeaccuracy,acc_trials(1).accuracy,'linewidth',4)
preds                               = squeeze(cat(3,data_trials_for_accuracy.success))';


%%
par.plot_AccuracyBars               = plot_AccuracyBarsParams;
par.plot_AccuracyBars.explained     = out.pcaModel.explained;
par.plot_AccuracyBars.evaltime      = evaltime;
par.plot_AccuracyBars.cmaps         = cmaps;
par.plot_AccuracyBars.channelSets   = par.BayesianInferenceClassification.channelSets;%{1:out.pcaModel.numComponents};
par.plot_AccuracyBars.nCols         = 1;

for iWindow=1:LastWindow
    idx_sel=iWindow:(iWindow+DT-1);
    % [data_trials_sel,res]                   = BayesianInferenceClassification(data_trials(idx_sel),par.BayesianInferenceClassification);
    hfig.AllComponentsAccuracyTrials(iWindow)= plot_AccuracyBars(data_trials_prob_raw(idx_sel),par.plot_AccuracyBars);
    % hfig.AllComponentsAccuracyTrials(iTrial)= plot_AccuracyBars(data_trials_sel,par.plot_AccuracyBars);
    sgtitle(['Window ', num2str(iWindow) ', t=' num2str(evaltime)]);
    timeTrials(:,iWindow)=idx_sel;
    % data_trials_sels{iTrial}=data_trials_sel;
end
%% BayesianInferenceClassifications separated components
par.BayesianInferenceClassificationSeparated            = BayesianInferenceClassificationParams();
par.BayesianInferenceClassificationSeparated.InField    = signal_process;
par.BayesianInferenceClassificationSeparated.exec       = true;
par.BayesianInferenceClassificationSeparated.channelSets= num2cell(1:out.pcaModel.numComponents); 
[data_trials_prob_sep,res]= BayesianInferenceClassification(data_trials,par.BayesianInferenceClassificationSeparated);

%% plot PROBABILITY - bar plot each dimension 
par.plot_EachDimBar                 = plot_EachDimBarParams;
par.plot_EachDimBar.novariance      = false;
par.plot_EachDimBar.addbar          = false;
par.plot_EachDimBar.cmaps           = cmaps;
par.plot_EachDimBar.legplot         = 1;
par.plot_EachDimBar.InField         = 'prob';
par.plot_EachDimBar.keep            = 1:3;
par.plot_EachDimBar.nCols           = 4;
par.plot_EachDimBar.ylabel          = '$P(C|x_{t})$';
par.plot_EachDimBar.YLIM            = [0-0.01,1+0.01];
par.plot_EachDimBar.evaltime        = 0;%i0+1;
par.plot_EachDimBar.chanceline      = true;
% GRAPH SUBPLOT TITLES - one graph for each subset
explained                           = out.pcaModel.explained;
channelSets                         = par.BayesianInferenceClassificationSeparated.channelSets;
nSets                               = length(channelSets); %% number of graphs
toleg=cell(nSets,1);
for iSet=1:nSets
    channels     =channelSets{iSet};
    nChannels    =length(channels);
    toleg{iSet}='';
    for ichannel=1:nChannels
        toleg{iSet} = sprintf('%s$${\\mathbf x}_{%d,t}$$',toleg{iSet},channels(ichannel));
        if ichannel<nChannels
             toleg{iSet} = sprintf('%s,',toleg{iSet});
        end 
    end
    if ~isempty(explained)
        toleg{iSet} = sprintf('%s (%2.1f)',toleg{iSet},sum(explained(channelSets{iSet})));
    end
end
par.plot_EachDimBar.titles          = toleg;

% meanData -> get class means
par.meanData                    = meanDataParams;
par.meanData.trialTypes         = [data_trials_prob_sep.trialType];
par.meanData.P                  = 90;
par.meanData.opt                = [1,0,0]; % bootstrap
par.meanData.InField            = 'prob';
par.meanData.OutField           = 'prob';
class_mean_trials               = meanData(data_trials_prob_sep,par.meanData);
% meanData -> get all trial means
par.meanData                    = meanDataParams;
par.meanData.P                  = 90;
par.meanData.opt                = [1,0,0]; % bootstrap
par.meanData.trialTypes         = true(size([data_trials_prob_sep.trialType]));
par.meanData.InField            = 'prob';
par.meanData.OutField           = 'prob';
mean_trials                     = meanData(data_trials_prob_sep,par.meanData);   
par.plot_EachDimBar.addbar      = mean_trials;
hfig.ComponentsProbBars         = plot_EachDimBar(class_mean_trials,par.plot_EachDimBar);

%% plot ACCURACY with plot_EachDimBar
par.plot_EachDimBar                 = plot_EachDimBarParams;
par.plot_EachDimBar.novariance      = false;
par.plot_EachDimBar.addbar          = false;
par.plot_EachDimBar.cmaps           = cmaps;
par.plot_EachDimBar.legplot         = 1;
par.plot_EachDimBar.InField         = 'accuracy';
par.plot_EachDimBar.novariance      = true;
par.plot_EachDimBar.keep            = 1:3;
par.plot_EachDimBar.nCols           = 4;
par.plot_EachDimBar.ylabel          = '$acc$';
par.plot_EachDimBar.YLIM            = [0-0.01,1+0.01];
par.plot_EachDimBar.evaltime        = i0+1;
par.plot_EachDimBar.chanceline      = true;
% GRAPH SUBPLOT TITLES - one graph for each subset
explained    =out.pcaModel.explained;
channelSets  =par.BayesianInferenceClassificationSeparated.channelSets;
nSets        =length(channelSets); %% number of graphs
toleg=cell(nSets,1);
for iSet=1:nSets
    channels     =channelSets{iSet};
    nChannels    =length(channels);
    toleg{iSet}='';
    for ichannel=1:nChannels
        toleg{iSet} = sprintf('%s$${\\mathbf x}_{%d,t}$$',toleg{iSet},channels(ichannel));
        if ichannel<nChannels
             toleg{iSet} = sprintf('%s,',toleg{iSet});
        end 
    end
    if ~isempty(explained)
        toleg{iSet} = sprintf('%s (%2.1f)',toleg{iSet},sum(explained(channelSets{iSet})));
    end
end
par.plot_EachDimBar.titles=toleg;

% meanData -> get class means
par.meanData                    = meanDataParams;
par.meanData.trialTypes         = [data_trials_prob_sep.trialType];
par.meanData.opt                = [0,0,1]; % mean
par.meanData.InField            = 'success';
par.meanData.OutField           = 'accuracy';
class_acc_trials                = meanData(data_trials_prob_sep,par.meanData);
% meanData -> get all trial means
par.meanData                    = meanDataParams;
par.meanData.opt                = [0,0,1]; % mean
par.meanData.trialTypes         = true(size([data_trials_prob_sep.trialType]));
par.meanData.InField            = 'success';
par.meanData.OutField           = 'accuracy';
acc_trials                      = meanData(data_trials_prob_sep,par.meanData);   
par.plot_EachDimBar.addbar      = acc_trials;%false;%data_trials_prob_sep(1);
hfig.ComponentsAccuracyBars     = plot_EachDimBar(class_acc_trials,par.plot_EachDimBar);

%%
evaltime                            = -900; % seconds

par.plot_AccuracyBars               = plot_AccuracyBarsParams;
par.plot_AccuracyBars.evaltime      = 1;
par.plot_AccuracyBars.explained     = out.pcaModel.explained;
par.plot_AccuracyBars.channelSets   = {1:out.pcaModel.numComponents};
par.plot_AccuracyBars.cmaps         = cmaps;
hfig.AllComponentsAccuracyBars      = plot_AccuracyBars(data_trials_prob(1:10),par.plot_AccuracyBars);

%% accuracy in trials, separated pcas
nBins                               = 5;
nTrials                             = length(data_trials);
DT                                  = nTrials-nBins+1;
LastWindow                          = nTrials-DT+1;
nWindows                            = length(1:LastWindow);
timeTrials                          = nan(DT,nWindows);
evaltime                            = -500; % seconds
par.plot_AccuracyBars               = plot_AccuracyBarsParams;
par.plot_AccuracyBars.evaltime      = evaltime;
par.plot_AccuracyBars.cmaps         = cmaps;
par.plot_AccuracyBars.channelSets   = par.BayesianInferenceClassificationSeparated.channelSets;% num2cell(1:out.pcaModel.numComponents);
par.plot_AccuracyBars.nCols         = min(4,out.pcaModel.numComponents);

for iTrial=1:LastWindow
    idx_sel=iTrial:(iTrial+DT-1);
    hfig.AllComponentsAccuracyTrials(iTrial) = plot_AccuracyBars(data_trials_prob_sep(idx_sel),par.plot_AccuracyBars);
    % [data_trials_sel,res]              = BayesianInferenceClassification(data_trials(idx_sel),par.BayesianInferenceClassification);
    sgtitle(['Trials ', num2str(iTrial) ', t=' num2str(evaltime)]);
    timeTrials(:,iTrial)=idx_sel;
    % data_trials_sels{iTrial}=data_trials_sel;
end
%%

%% QDA
par.TimeSelect                  = TimeSelectParams;
par.TimeSelect.t1               = -40;%-520;
par.TimeSelect.t2               = -20;%-500;
par.TimeSelect.InField          = signal_process;
par.TimeSelect.OutField         = signal_process;
data_trials_W1                  = TimeSelect(bootdata_trials,par.TimeSelect);
% AverageWindow
par.AverageWindow               = AverageWindowParams;
par.AverageWindow.useSqrt       = false;
par.AverageWindow.InField       = signal_process;
par.AverageWindow.OutField      = signal_process;
par.AverageWindow.binWidth      = size(data_trials_W1(1).(['time' signal_process]),2);
data_trials_W1                  = AverageWindow(data_trials_W1,par.AverageWindow);
% dataSplit
% par.dataSplit                   = dataSplitParams;
% par.dataSplit.TrainPercentage   = 60;
% data_trials_W1               = dataSplit(data_trials_W1,par.dataSplit);
% fitQDA
par.fitQDA                      = fitQDAParams;
par.fitQDA.InField              = signal_process;
% [~,out.fitQDA]                  = fitQDA(data_trials_qda([data_trials_qda.train]),par.fitQDA);
[~,out.fitQDA]                  = fitQDA(data_trials_W1,par.fitQDA);
disp(out.fitQDA);
% TimeSelect W2
par.TimeSelect                  = TimeSelectParams;
par.TimeSelect.t1               = -20;%-500;
par.TimeSelect.t2               = 0;%-480;
par.TimeSelect.InField          = signal_process;
par.TimeSelect.OutField         = signal_process;
% data_trials_W2                  = TimeSelect(bootdata_trials,par.TimeSelect);
data_trials_W2                  = TimeSelect(data_trials,par.TimeSelect);
% AverageWindow
par.AverageWindow               = AverageWindowParams;
par.AverageWindow.useSqrt       = false;
par.AverageWindow.InField       = signal_process;
par.AverageWindow.OutField      = signal_process;
par.AverageWindow.binWidth      = size(data_trials_W2(1).(['time' signal_process]),2);
data_trials_W2                  = AverageWindow(data_trials_W2,par.AverageWindow);
% predictQDA
par.predictQDA                  = predictQDAParams;
par.predictQDA.InField          = signal_process;
par.predictQDA.OutField         = 'QDA';
par.predictQDA.mdl              = out.fitQDA.mdl;
[data_trials_qda_test,res]      = predictQDA(data_trials_W2,par.predictQDA);
disp(res);
% TimeSelect
par.TimeSelect                  = TimeSelectParams;
par.TimeSelect.t1               = -1000;
par.TimeSelect.t2               = 0;
par.TimeSelect.InField          = signal_process;
par.TimeSelect.OutField         = signal_process;
% functions to execute
par.exec={'TimeSelect','AverageWindow','fitQDA'};

data_trials_to_fit              = bootdata_trials;
[~,res]                         = fitQDA(data_trials_to_fit,par.fitQDA);

% TimeSelect
par.TimeSelect                  = TimeSelectParams;
par.TimeSelect.t1               = -500;
par.TimeSelect.t2               = 0;
par.TimeSelect.InField          = signal_process;
par.TimeSelect.OutField         = signal_process;
% functions to execute
par.exec={'TimeSelect','AverageWindow','predictQDA'};

% predictQDA
par.predictQDA.exec             = true;                 
par.predictQDA.mdl              = res.mdl; 
par.fitQDA.InField              = signal_process;
par.fitQDA.OutField             = 'QDA';
[data_trials_to_fit,res]        = predQDA(data_trials_to_pred,par.fitQDA);



