% function TEST_SHORTCUT_NSA()
% SHORTCUT DECISION, for methods see DPCA  
% Pezzulo, G., Donnarumma, F., Ferrari-Toniolo, S., Cisek, P., & Battaglia-Mayer, A. (2022). 
% Shared population-level dynamics in monkey premotor cortex during solo action, joint action and action observation. 
% Progress in Neurobiology, 210, 102214.
%%
clear;
par.irng                        = 10;                  % for reproducibility
rng(par.irng);
%% Step 0. Load data in raw format 
% RatName              ='Rat067'; % subject
RatName              ='Rat066'; % subject
% RatName              ='Rat068'; % subject
runs=11;             % selected options (see params);
[S_cell,params]      = RatNSA_main(RatName,runs);

params.AB            = true;

Trials               = S_cell.Trials;
Data                 = S_cell.Data;
Encoder              = S_cell.Encoder;
ratname              = Data.ratname;
iday                 = Data.day;
phase                = Data.phase;
kernSD               = params.kernSD;
iDecision            = params.iDecision;
prolong              = params.prolong;
condPerTrial         = Trials.Traj_labels(Trials.ABlabels==params.AB,:);

conditions           = unique(condPerTrial);
NConditions          = length(conditions);
NeuronLabels         = Data.Spikes.label;
NNeurons             = length(NeuronLabels);
Tasklabels           = Trials.labels;
gap_after            = params.data.gap_after;
gap_before           = params.data.gap_before;
questionable_cells   = params.data.load_questionable_cells;

Traj_Decision        = Trials.Traj_Decision{iDecision}(Trials.ABlabels==params.AB,:);
Trials.Traj_edges    = Trials.Traj_edges(Trials.ABlabels==params.AB,:);

Trials.Traj_indexes  = Trials.Traj_indexes(Trials.ABlabels==params.AB,:);
Trials.Traj_labels   = Trials.Traj_labels(Trials.ABlabels==params.AB,:);
for iD=1:length(Trials.Traj_Decision)
    Trials.Traj_Decision{iD}   = Trials.Traj_Decision{iD}(Trials.ABlabels==params.AB,:);
    Trials.T_Traj_Decision{iD} = Trials.T_Traj_Decision{iD}(Trials.ABlabels==params.AB,:);
end
Trials.ABlabels      = Trials.ABlabels(Trials.ABlabels==params.AB,:);
NTrialsOriginal      = length(Traj_Decision);


endTR                = gap_after+gap_before; % trial duration in ms
if params.AB
    ABstr='AB';
else
    ABstr='BA';
end
description          = sprintf('%s_day%g_ph%g_D%g_%s_Q%g_b%ga%g_KRN%g',ratname,           ...
                                                                    iday,              ...
                                                                    phase,             ...
                                                                    iDecision,         ...
                                                                    ABstr,             ...
                                                                    questionable_cells,...
                                                                    gap_before,        ...
                                                                    gap_after,                                                                                      ...
                                                                    kernSD);
disp(description);
runIdx              = params.runIdx;
dat_name            = description;
interval            = Data.T(Traj_Decision(:,[1,3]));
exclude_trials      = ~(abs(1000*diff(interval')-endTR)<exp(-10));
nTrials             = sum(~exclude_trials);
goodTrials          = find(~exclude_trials);
        
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
freq_interval_in_bins=20;                           % freq on 20 bins
dt_in_s              = 1/1000;                      % 1ms bin
Treset               = -gap_before;                 % trial zero time  

Spikes.t             = cell(nTrials,NNeurons);
Spikes.interval(exclude_trials,:)=[];
condPerTrial(exclude_trials)=[];
% loop on trials to select the chosen time window around the decision point [-gap_before,gap_after]
for iTrial=1:nTrials   
    for iNeuron=1:NNeurons
        interval           =Spikes.interval(iTrial,:);
        spike_times        =Data.Spikes.t{iNeuron};
        sel_logical        =spike_times>interval(1) & spike_times< interval(2);            
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
    data_trials(it).([xfld 'spikes'])   = -gap_before:(dt_in_s*1000):gap_after; %ms
end

%% Step 2. Smooth for dpca processing
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
[data_trials, outStep2]         = run_trials(data_trials,par);

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
% pcaCompute
par.pcaModel                    = pcaModelParams();
% par.pcaCompute.numComponents    = 8;
par.pcaModel.numComponents      = 0;
par.pcaModel.perc               = 90;
par.pcaModel.InField            = signal_process;
par.pcaModel.OutField           = signal_process;
%%%%%%%%%%%%% nsa_pca 2-Stage Engine Churchland : kernel smooth + pca -> %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% 'AverageWindow','GaussianSmoother','pcaCompute'            %%%%%%%%%%%%%%%%%%
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
% par.exec.funname                = {'TimeSelect','meanData','AverageWindow','GaussianSmoother','pcaCompute'};
par.exec.funname                = {'meanData','AverageWindow','GaussianSmoother','pcaModel'};
[data_trials_class, out]        = run_trials(data_trials,par);

% pcaProject
par.pcaEncode.Wpca              = out.pcaModel.Wpca;
par.pcaEncode.mu                = out.pcaModel.mu;
par.pcaEncode.explained         = out.pcaModel.explained;
par.pcaEncode.InField           = signal_process;
par.pcaEncode.OutField          = signal_process;
par.exec.funname                = {'pcaEncode'};
data_trials_class               = run_trials(data_trials_class,par);

%% Step 4. Project trials on the dictionary found 
% pcaProject
% par.pcaProject.Wpca             = out.pcaCompute.W;
% par.pcaProject.mu               = out.pcaCompute.mu;
% par.pcaProject.explained        = out.pcaCompute.explained;
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
[data_trials, out2]             = run_trials(data_trials,par);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Some plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmaps                           = linspecer(length(unique([data_trials_class.trialType]))); % u-path, shortcut, deadend
timeField                       = data_trials_class(1).(['time' par.pcaEncode.OutField]);
[~,i0]                          = min(abs(timeField));      % select   0  ms % TargetOn

%% Bootstrap
par.BootStrapData               = BootStrapDataParams;
par.BootStrapData.N             = 100;
par.BootStrapData.InField       = signal_process;
par.BootStrapData.OutField      = signal_process;
bootdata_trials                 = BootStrapData(data_trials,par.BootStrapData);

%% plot 2D Trajectories
i500                            = find(timeField>=-500,1);   % select -250  ms 
% plot_trajectory2D
par.plot_trajectory2D           = plot_trajectory2DParams;
par.plot_trajectory2D.keep      = 1:3;                      % which conditions u-path, shortcut, deadend
par.plot_trajectory2D.wd        = [1,2];                    % which components
par.plot_trajectory2D.axlab     = 'x';
par.plot_trajectory2D.cmaps     = cmaps;
par.plot_trajectory2D.cmapslight= lightCmaps(par.plot_trajectory2D.cmaps);
par.plot_trajectory2D.explained = out.pcaModel.explained;
par.plot_trajectory2D.InField   = par.pcaEncode.OutField;  % take pca projections to show
par.plot_trajectory2D.istart    = 1;
par.plot_trajectory2D.center    = i0;
par.plot_trajectory2D.iend      = length(timeField);
hfg.trajectory2D                = plot_trajectory2D(data_trials_class, par.plot_trajectory2D);
hfg.trajectory2D                = plot_trajectory2D(meanData(bootdata_trials,par.meanData), par.plot_trajectory2D);
hfg.trajectory2D                = plot_trajectory2D(meanData(data_trials,par.meanData), par.plot_trajectory2D);

%% Plot 3D Trajectories  
par.plot_trajectory3D           = plot_trajectory3DParams;
par.plot_trajectory3D.keep      = 1:3;                      % which conditions u-path, shortcut, deadend
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
par.plot_ellipsoid2D            = CI_computeParams;
par.plot_ellipsoid2D.InField    = par.pcaEncode.OutField;
par.plot_ellipsoid2D.P          = 68;
par.plot_ellipsoid2D.opt        = [0,1,0];                  % bootstrap
par.plot_ellipsoid2D.wd         = [1,2];                    % which components

selected                        = 1:3;
T                               = length(timeField);
dt                              = 5;
indsT                           = 1:dt:T;
par.plot_ellipsoid2D.explained  = out.pcaModel.explained;
hfg.ellipsoid2D                 = plot_ellipsoid2D(data_trials,selected,indsT,cmaps,par.plot_ellipsoid2D);

%% plot pcas dirs on conditions with plot_Each_mean - showing mean and confidence intervals
par.plot_EachDimVsTime              = plot_EachDimVsTimeParams;
par.plot_EachDimVsTime.cmaps        = cmaps;
par.plot_EachDimVsTime.cmapslight   = lightCmaps(par.plot_EachDimVsTime.cmaps);
par.plot_EachDimVsTime.legplot      = 2;
par.plot_EachDimVsTime.keep         = 1:3;
par.plot_EachDimVsTime.InField      = signal_process;
par.plot_EachDimVsTime.ylabel       = '$pca$';
par.plot_EachDimVsTime.nCols        = 4;
% GRAPH SUBPLOT TITLES
explained       =out.pcaModel.explained;
channels        =1:out.pcaModel.numComponents;
nChannels       =length(channels); %% number of graphs
toleg             =cell(nChannels,1);
for ichannel=1:nChannels
    toleg{ichannel} = sprintf('$$\\tilde{\\mathbf x}_{%d,t}$$',ichannel);
    if ~isempty(explained)
        toleg{ichannel} = sprintf('%s (%2.1f)',toleg{ichannel},explained(ichannel));
    end
end
par.plot_EachDimVsTime.titles=toleg;

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
par.pSeparability.exec              = true;
pdata_trials                        = bootdata_trials; % data_trials
[pVals,pComps]                      = pSeparability(pdata_trials,par.pSeparability);
%% pvalue plot per feature
par.plot_pValues.InField     = 'comparisons';
par.plot_pValues.xfld        = 'time';
par.plot_pValues.dt          = 50;

hfg.pvals       = plot_pValues(pComps,par.plot_pValues);
return
par.InField     = 'comparisons';
par.xfld        = 'time';
close all;
hfg.pvals       = figure; 
xfld            = par.xfld;
time            = pdata_trials(1).([xfld signal_process]);
dt              = 50; % ms
inds            = time>time(1)+dt & time<time(end)-dt;
nClasses        = length(pComps);
InField         = par.InField;
newk            = 1:nChannels;
nRows           = 2;
nCols           = 2;%3;
legplot         = 2;
lw              = 4;
cc              = [];    
if legplot==2
    hSub        = subplot(nRows+1, 1, 1);
    set(gca,'xticklabel','');set(gca,'yticklabel','');
    axis off;
    newk        = indexesForSubplotsWithLegend2(nRows,nCols);
    nRows_mod   = nRows+1;
    nCols_mod   = nCols;
elseif legplot==1
    hSub    = subplot(1, nCols+1, 1);
    set(gca,'xticklabel','');set(gca,'yticklabel','');
    axis off;
    newk=indexesForSubplotsWithLegend(nRows,nCols);
    nRows_mod=nRows;
    nCols_mod=nCols+1;
else
    nRows_mod=nRows;
    nCols_mod=nCols;
end
selY=newk(1):nCols_mod:max(newk);
selX=newk(end-nCols+1):newk(end);
nComps  = nchoosek(nClasses,2);
% cmaps   = linspecer(nComps);
cmaps   = linspecer(nClasses+nComps+1);
cmaps   = cmaps(nClasses+1:end,:);
limval  = [inf,-inf];
totit   = cell(1,nChannels);
for ichannel=1:nChannels
    totit{ichannel} = sprintf('$$\\tilde{\\mathbf x}_{%d,t}$$',ichannel);
    if ~isempty(explained)
        totit{ichannel} = sprintf('%s (%2.1f)',totit{ichannel},explained(ichannel));
    end
end
for iChannel=1:nChannels
    subplot(nRows_mod,nCols_mod,newk(iChannel))
    hold on; box on; grid on;
    set(gca,'fontsize',13)
    cc=[];
    toleg=cell(0,0);
    indcol  = 0;
    xline(0,'label','Decision Point','LabelVerticalAlignment','bottom','LineWidth',5)
    yline(0.05,'LineStyle','--');
    yline(0.1,'LineStyle','--');
    yline(0.5,'LineStyle','--');
    if false
        cc(end+1)=plot(time,pVals(iChannel,:),'Color',cmaps(end,:),'LineWidth',lw);
        trialNames = {pComps.trialName};
        toleg{end+1}=trialNames{1};
        for itn=2:length(trialNames)
            toleg{end}=sprintf('%s vs %s',toleg{end},trialNames{itn});
        end
        % toleg{end+1}=sprintf('%s vs %s',pComps(iClass).trialName,pComps(ic).trialName);
        limval(1)=min([limval(1);pVals(:)]);
        limval(2)=max([limval(2);pVals(:)]);
    else
        for iClass=1:nClasses
            for ic=1:nClasses
                if ic<=iClass
                    continue
                else
                    indcol=indcol+1;
                    % col = cmaps(iClass,:)+cmaps(ic,:);
                    % if any(col>1)
                    %     col=col/max(col);
                    % end
                    % col=col*0.8;
                    % cc(end+1)=plot(time,comp,'Color',col,'LineWidth',lw);
                    comp        = squeeze(pComps(iClass).(InField)(iChannel,:,ic));
                    comp        = smooth(comp);
                    cc(end+1)   = plot(time(inds),comp(inds),'Color',cmaps(indcol,:),'LineWidth',lw);
                    % pVals
                    toleg{end+1}=sprintf('%s vs %s',pComps(iClass).trialName,pComps(ic).trialName);
                    limval(1)=min([limval(1);comp(:)]);
                    limval(2)=max([limval(2);comp(:)]);
        
                end
            end
        end
    end
    title(totit{iChannel},'Interpreter','latex')
    % set(gca,'yscale','log')

    if ismember(newk(iChannel),selX)
        xlabel('time [ms]');
    end
    if ismember(newk(iChannel),selY)
        % ylabel(par.ylabel,'interpreter','latex');
        ylabel('p (log scale)')
    end  
    
    yticks([0,0.01,0.05,0.1,0.5,1]);
    % legend(cc,toleg);
end
for iChannel=1:nChannels
    subplot(nRows_mod,nCols_mod,newk(iChannel))
    hold on; box on; grid on;
    set(gca,'yscale','log');
    ylim(limval);
    % 
end

if legplot>0
  hLegend   = legend(cc,toleg);
  set(hLegend, 'position', get(hSub, 'position'));
  uistack(hLegend,'bottom');
end
if params.AB
    feedstr=sprintf('feeder A->B');
else 
    feedstr=sprintf('feeder B->A');
end
daystr = sprintf('Day%g',params.data.day);
sgtitle([RatName ' ' daystr ' ' feedstr]);
p=[0,0,1500,600];
set(hfg.pvals,'Position',p);

%% plot hfg

return
%% BayesianInferenceClassificator, cumulated pca components 
par.BayesianInferenceClassification            = BayesianInferenceClassificationParams();
par.BayesianInferenceClassification.InField    = signal_process;
par.BayesianInferenceClassification.exec       = true;
par.BayesianInferenceClassification.channelSets= {1:out.pcaCompute.numComponents};        % all pcas together
% par.BayesianInferenceClassification.channelSets=num2cell(1:out.pcaCompute.numComponents); % all pcas separated
% par.BayesianInferenceClassification.channelSets={1,2,1:3,4:6};                            % some pca combinations
[data_trials_prob,res]                         = BayesianInferenceClassification(TimeSelect(data_trials,par.TimeSelect),par.BayesianInferenceClassification);

%%
data_to_class                                  = TimeSelect(bootdata_trials,par.TimeSelect);
par.BayesianInferenceClassification.train      = [data_to_class.train];
[data_trials_prob_mdl,res]                     = BayesianInferenceClassification(data_to_class,par.BayesianInferenceClassification);
par.BayesianInferenceClassification.train      = 0;
[data_trials_prob,res]                         = BayesianInferenceClassification(data_trials_prob_mdl,par.BayesianInferenceClassification);
%%
% par.BayesianInferenceClassification.train       = [data_trials.train];
% [data_trials_prob_mdl,res]                      = BayesianInferenceClassification(data_trials,par.BayesianInferenceClassification);
% par.BayesianInferenceClassification.train       = 0;
% [data_trials_prob,res]                          = BayesianInferenceClassification(data_trials_prob_mdl,par.BayesianInferenceClassification);

%% plot classification probabilites in time
% plot_EachDimVsTime
par.plot_EachDimVsTime              = plot_EachDimVsTimeParams;
par.plot_EachDimVsTime.cmaps        = cmaps;
par.plot_EachDimVsTime.cmapslight   = lightCmaps(par.plot_EachDimVsTime.cmaps);
par.plot_EachDimVsTime.InField      = 'prob';
par.plot_EachDimVsTime.ylabel       = '$P(C|o_{t})$';
par.plot_EachDimVsTime.keep         = 1:3;
par.plot_EachDimVsTime.nCols        = 1;
par.plot_EachDimVsTime.YLIM         = [0-0.01,1+0.01];
par.plot_EachDimVsTime.legplot      = 2;

% GRAPH SUBPLOT TITLES
explained    =out.pcaCompute.explained;
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

hfig.CumPcaClassificationInTime = plot_EachDimVsTime(meanData(data_trials_prob,par.meanData),par.plot_EachDimVsTime);

%% BayesianInferenceClassifications separated components
par.BayesianInferenceClassificationSeparated            = BayesianInferenceClassificationParams();
par.BayesianInferenceClassificationSeparated.InField    = signal_process;
par.BayesianInferenceClassificationSeparated.exec       = true;
par.BayesianInferenceClassificationSeparated.channelSets= num2cell(1:out.pcaCompute.numComponents); 
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
explained                           = out.pcaCompute.explained;
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
explained    =out.pcaCompute.explained;
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
par.plot_AccuracyBars.explained     = out.pcaCompute.explained;
par.plot_AccuracyBars.channelSets   = {1:out.pcaCompute.numComponents};
par.plot_AccuracyBars.cmaps         = cmaps;
hfig.AllComponentsAccuracyBars      = plot_AccuracyBars(data_trials_prob(1:10),par.plot_AccuracyBars);

%% accuracy in trials, separates pcas
nBins                               = 5;
nTrials                             = length(data_trials);
DT                                  = nTrials-nBins+1;
LastWindow                          = nTrials-DT+1;
nWindows                            = length(1:LastWindow);
timeTrials                          = nan(DT,nWindows);
evaltime                            = -500; % seconds
par.plot_AccuracyBars               = plot_AccuracyBarsParams;
par.plot_AccuracyBars.evaltime      = evaltime;
par.plot_AccuracyBars.channelSets   = par.BayesianInferenceClassificationSeparated.channelSets;% num2cell(1:out.pcaCompute.numComponents);
par.plot_AccuracyBars.nCols         = min(4,out.pcaCompute.numComponents);

for iTrial=1:LastWindow
    idx_sel=iTrial:(iTrial+DT-1);
    hfig.AllComponentsAccuracyTrials(iTrial) = plot_AccuracyBars(data_trials_prob_sep(idx_sel),par.plot_AccuracyBars);
    % [data_trials_sel,res]              = BayesianInferenceClassification(data_trials(idx_sel),par.BayesianInferenceClassification);
    sgtitle(['Trials ', num2str(iTrial) ', t=' num2str(evaltime)]);
    timeTrials(:,iTrial)=idx_sel;
    % data_trials_sels{iTrial}=data_trials_sel;
end
%%

%% accuracy in trials
nBins                               = 8;
nTrials                             = length(data_trials);
DT                                  = nTrials-nBins+1;
LastWindow                          = nTrials-DT+1;
nWindows                            = length(1:LastWindow);
timeTrials                          = nan(DT,nWindows);
evaltime                            = 0; % seconds
par.plot_AccuracyBars               = plot_AccuracyBarsParams;
par.plot_AccuracyBars.explained     = out.pcaCompute.explained;
par.plot_AccuracyBars.evaltime      = evaltime;
par.plot_AccuracyBars.channelSets   = par.BayesianInferenceClassification.channelSets;%{1:out.pcaCompute.numComponents};
par.plot_AccuracyBars.nCols         = 1;

for iTrial=1:LastWindow
    idx_sel=iTrial:(iTrial+DT-1);
    % [data_trials_sel,res]                   = BayesianInferenceClassification(data_trials(idx_sel),par.BayesianInferenceClassification);
    hfig.AllComponentsAccuracyTrials(iTrial)= plot_AccuracyBars(data_trials_prob(idx_sel),par.plot_AccuracyBars);
    % hfig.AllComponentsAccuracyTrials(iTrial)= plot_AccuracyBars(data_trials_sel,par.plot_AccuracyBars);
    sgtitle(['Trials ', num2str(iTrial) ', t=' num2str(evaltime)]);
    timeTrials(:,iTrial)=idx_sel;
    % data_trials_sels{iTrial}=data_trials_sel;
end
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



