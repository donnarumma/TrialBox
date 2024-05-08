function par = TEST_SHORTCUT_TRAJ(ifplot,root_dir,AB,RatName,irng)
% function TEST_SHORTCUT_TRAJ(ifplot)
% SHORTCUT DECISION, for methods see DPCA  
% Pezzulo, G., Donnarumma, F., Ferrari-Toniolo, S., Cisek, P., & Battaglia-Mayer, A. (2022). 
% Shared population-level dynamics in monkey premotor cortex during solo action, joint action and action observation. 
% Progress in Neurobiology, 210, 102214.
%%
% clear;
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
%% Step 0. Load data in raw format 
runs=11;             % selected options (see params);
[S_cell,params]      = RatNSA_main(RatName,runs);

params.AB            = AB;
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

% conditions           = unique(condPerTrial);
% NConditions          = length(conditions);
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
if ~exist('root_dir','var')
    root_dir = '';
end
save_dir          = [root_dir description filesep];

% runIdx              = params.runIdx;
% dat_name            = description;
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
% freq_interval_in_bins=20;                           % freq on 20 bins
dt_in_s              = 1/1000;                      % 1ms bin
% Treset               = -gap_before;                 % trial zero time  

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

%% plot 2D reference
params.plots.hfig           = figure('visible',ifplot);
hfg.linearised2Dreference   = params.plots.linearised2Dreference(S_cell.Data,S_cell.Encoder,params.plots);

%% Plot Summary Trial Conditions
hfig                = figure('visible',ifplot);
params.plots.hfig   = hfig;
hfg.TrialConditions = plot_TrialConditions(data_trials,params.plots);
sgtitle(strrep(description,'_',' '),'fontsize',16);

%% Plot Trial Trajectories
for iTrial=1:nTrials
    itraj               = goodTrials(iTrial);
    hfig                = figure('visible',ifplot);
    params.plots.hfig   = hfig;
    hfg.(['TrialTrajectory' num2str(iTrial)])=params.plots.TrialTrajectory(Data,Trials,[itraj,0],params.plots); % 0 overplot the decision points
end

%% Plot Decision Trajectories
idp=params.iDecision;
    % params for "fast" (swr) decoder
params.decoder.t_intervals_width        = 0.025;
params.decoder.t_intervals_overlap      = 0.005;

for iTrial=1:nTrials
    itraj               = goodTrials(iTrial);
    hfig                = figure('visible',ifplot);
    hfg.(['SpikeOrderedsvsMotorDecodedTraj' num2str(iTrial)]) = hfig;
    params.plots.hfig   = hfig;        
    plot_SpikeOrderedsvsMotorDecoded(Data,Encoder,Trials,[itraj,idp,prolong],params.decoder,params.plots);
    
    labh=strsplit(hfig.Children(1).String,':');
    splb=strsplit(labh{1},'Traj');
    lab = [splb{1} 'Trial' num2str(iTrial) ':' labh{2}];% 'Decision'];
    set(hfig, 'PaperPositionMode','auto');
    set(hfig, 'PaperPositionMode','auto');
    sgtitle(lab)
end

%% plot hfg
if ~ifplot
    par.hfigPrint               = hfigPrintParams();
    par.hfigPrint.pdf_file      = [root_dir mfilename '_' description '.pdf'];
    par.hfigPrint.wm            = 170;
    par.hfigPrint.save_dir      = save_dir; 
    hfigPrint(hfg,par.hfigPrint)
end