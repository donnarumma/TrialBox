function par=TEST_CLOCK_CEBRA_NSA(ifplot,irng)
% see PCA Clock, 
% Pezzulo, G., Donnarumma, F., Ferrari-Toniolo, S., Cisek, P., & Battaglia-Mayer, A. (2022). 
% Shared population-level dynamics in monkey premotor cortex during solo action, joint action and action observation. 
% Progress in Neurobiology, 210, 102214.
% alternate approach in CEBRA
% clear;
try 
    irng;
catch
    irng=10;
end
try
    ifplot;
catch
    ifplot                      = true;             
end
reload                          = false;
max_iterations                  = 10000;
output_dimension                = 3;
par.irng                        = irng;   % for reproducibility
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
signal_process                  = 'neural';       % data processed are in 'y' field
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

%%
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

% functions to be execute
par.exec.funname                = {'AverageWindow','GaussianSmoother'};
[data_trials, out2]             = run_trials(data_trials,par);

%% Step 2. perform cebra model
cebra_codes_dir                     = '~/tools/TrialBox/pywraps/'; 
model_name                          = 'cebra_model.pkl'; 
if reload
    % model_dir     = '/TESTS/CEBRA/SAPIENZA/cebra20240420080115/'; output_dimension = 3;
    model_dir       = '/TESTS/CEBRA/SAPIENZA/cebra20240420083547/'; output_dimension = 8;

    model_filename  = [model_dir, model_name];
else
    par.cebraModel                      = cebraModelParams();
    par.cebraModel.seed                 = par.irng;
    nTimes                              = length(data_trials(1).(['time'  signal_process])); 
    for iTrial=1:nTrials
        data_trials(iTrial).behavior    = data_trials(iTrial).trialType*ones(1,nTimes);
    end
    strdate                             = datetime('now', 'Format', 'yyyyMMddHHmmss');
    script_rundir                       = ['/tmp/cebra' sprintf('%s',strdate) '/'];
    mkdir(script_rundir);
    model_dir                           = ['/home/donnarumma/TESTS/CEBRA/SAPIENZA/cebra' sprintf('%s',strdate) '/'];
    mkdir(model_dir);
    % path to python scripts
    par.cebraModel.script_filename      = [cebra_codes_dir 'cebraModel.py'];
    % path to hd5 files
    par.cebraModel.script_rundir        = script_rundir;
    par.cebraModel.model_filename       = [model_dir model_name];
    % other parameteres
    par.cebraModel.max_iterations       = max_iterations;
    par.cebraModel.output_dimension     = output_dimension;
    disp(par.cebraModel);
    [~,out.cebraModel]                  = cebraModel(data_trials,par.cebraModel);
    model_filename                      = par.cebraModel.model_filename;
    output_dimension                    = par.cebraModel.output_dimension;
end
%% Step 3. cebraEncode -> project data on the manifold
strdate                             = datetime('now', 'Format', 'yyyyMMddHHmmss');
script_rundir                       = ['/tmp/cebra' sprintf('%s',strdate) '/'];
mkdir(script_rundir);
par.cebraEncode                     = cebraEncodeParams();
par.cebraEncode.model_filename      = model_filename;
par.cebraEncode.script_rundir       = script_rundir;
par.cebraEncode.script_filename     = [cebra_codes_dir 'cebraEncode.py'];
disp(par.cebraEncode)
data_trials                         = cebraEncode(data_trials,par.cebraEncode);

% Some plots
%% plot_trajectory2D
% meanData
par.meanData                    = meanDataParams;
par.meanData.ifclass            = true;
par.meanData.trialTypes         = [data_trials.trialType];
par.meanData.InField            = par.cebraEncode.manifold_field;
par.meanData.OutField           = par.cebraEncode.manifold_field;
par.meanData.P                  = 80;
par.meanData.opt                = [0,0,1];
data_trials_means               = meanData(data_trials, par.meanData);%

cmaps                           = linspecer(24);            % 8 directions x 3 conditions
cmaps                           = cmaps([17:24,1:8,9:16],:);% reset color as desired 
timeField                       = data_trials(1).(['time' par.cebraEncode.manifold_field]);
[~,i0]                          = min(abs(timeField));      % select   0  ms % TargetOn
istart                          = find(timeField>=-80,1);   % select -80  ms 
iend                            = find(timeField>=250,1);   % select 250  ms
par.plot_trajectory2D           = plot_trajectory2DParams;
par.plot_trajectory2D.keep      = 9:16;                     % which conditions (e.g. obs)
par.plot_trajectory2D.wd        = [1,2];                    % which components
par.plot_trajectory2D.axlab     = 'x';
par.plot_trajectory2D.cmaps     = cmaps;
par.plot_trajectory2D.cmapslight= lightCmaps(par.plot_trajectory2D.cmaps);
par.plot_trajectory2D.explained = [];
par.plot_trajectory2D.InField   = par.cebraEncode.manifold_field;  % take manifold projections to show
par.plot_trajectory2D.istart    = istart;
par.plot_trajectory2D.center    = i0;
par.plot_trajectory2D.iend      = iend;
hfg.plot_trajectory2D           = plot_trajectory2D(data_trials_means, par.plot_trajectory2D);

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
par.plot_EachDimVsTime.hfig         = figure('visible',ifplot);
par.plot_EachDimVsTime.cmaps        = cmaps;
par.plot_EachDimVsTime.cmapslight   = lightCmaps(par.plot_EachDimVsTime.cmaps);
par.plot_EachDimVsTime.InField      = par.cebraEncode.manifold_field;  % take manifold projection to show;
par.plot_EachDimVsTime.ylabel       = '$pca$';
par.plot_EachDimVsTime.keep         = 1:nConditions_m;
par.plot_EachDimVsTime.nCols        = 4;
par.plot_EachDimVsTime.legplot      = 2;
% GRAPH SUBPLOT TITLES
explained                           = [];%out.pcaCompute.explained;
channels                            = 1:output_dimension;
nChannels                           = length(channels); %% number of graphs
str                                 = cell(nChannels,1);
for ichannel=1:nChannels
    str{ichannel} = sprintf('$$\\tilde{\\mathbf x}_{%d,t}$$',ichannel);
    if ~isempty(explained)
        str{ichannel} = sprintf('%s (%2.1f)',str{ichannel},explained(ichannel));
    end
end
par.plot_EachDimVsTime.titles   = str;

par.meanData                    = meanDataParams;
par.meanData.ifclass            = true;
par.meanData.trialTypes         = [data_trials_conds.trialType];
par.meanData.InField            = par.cebraEncode.manifold_field;
par.meanData.OutField           = par.cebraEncode.manifold_field;
par.meanData.P                  = 68;
par.meanData.opt                = [1,0,0];
data_trials_means               = meanData(data_trials_conds, par.meanData);%

hfg.pcas                        = plot_EachDimVsTime(data_trials_means,par.plot_EachDimVsTime);

%% plot pcas dirs on conditions with plot_Each_mean - showing mean and confidence intervals
% remapTypes
par.remapTypes                  = remapTypesParams;
par.remapTypes.selection        = { [ 1, 9,17] , [2,10,18], [3,11,19], [4,12,20], ... 
                                    [ 5,13,21] , [6,14,22], [7,15,23], [8,16,24] };
par.remapTypes.names            = { 'Dir1', 'Dir2', 'Dir3', 'Dir4', ...
                                    'Dir5', 'Dir6', 'Dir7', 'Dir8'};
data_trials_dirs                = remapTypes(data_trials,par.remapTypes);
% plot_Each_mean
conditions_m                    = unique([data_trials_dirs.trialType]);
nConditions_m                   = length(conditions_m);
cmaps                           = linspecer(nConditions_m); % 
%cmaps                           = cmaps([2,1,3],:);         % color map for joint, obs and solo

% plot_EachDimVsTime
par.plot_EachDimVsTime              = plot_EachDimVsTimeParams;
par.plot_EachDimVsTime.hfig         = figure('visible',ifplot);
par.plot_EachDimVsTime.cmaps        = cmaps;
par.plot_EachDimVsTime.cmapslight   = lightCmaps(par.plot_EachDimVsTime.cmaps);
par.plot_EachDimVsTime.InField      = par.cebraEncode.manifold_field;  % take manifold projection to show;
par.plot_EachDimVsTime.ylabel       = '$cebra$';
par.plot_EachDimVsTime.keep         = 1:nConditions_m;
par.plot_EachDimVsTime.nCols        = 4;
par.plot_EachDimVsTime.legplot      = 2;
% GRAPH SUBPLOT TITLES
explained                           = [];%out.pcaCompute.explained;
channels                            = 1:output_dimension;
nChannels                           = length(channels); %% number of graphs
str                                 = cell(nChannels,1);
for ichannel=1:nChannels
    str{ichannel} = sprintf('$$\\tilde{\\mathbf x}_{%d,t}$$',ichannel);
    if ~isempty(explained)
        str{ichannel} = sprintf('%s (%2.1f)',str{ichannel},explained(ichannel));
    end
end
par.plot_EachDimVsTime.titles   = str;

par.meanData                    = meanDataParams;
par.meanData.ifclass            = true;
par.meanData.trialTypes         = [data_trials_dirs.trialType];
par.meanData.InField            = par.cebraEncode.manifold_field;
par.meanData.OutField           = par.cebraEncode.manifold_field;
par.meanData.P                  = 68;
par.meanData.opt                = [1,0,0];
data_trials_means               = meanData(data_trials_dirs, par.meanData);%

hfg.pcas                        = plot_EachDimVsTime(data_trials_means,par.plot_EachDimVsTime);

%% plot_scatterGradient
data_trials_conds_scatter       = data_trials_conds;
InField                         = par.cebraEncode.manifold_field;
nTrials                         = length(data_trials_conds_scatter);
for iTrial=1:nTrials
    time    = data_trials_conds_scatter(iTrial).(['time' InField]);
    nTimes  = length(time);
    tType   = data_trials_conds_scatter(iTrial).trialType;
    data_trials_conds_scatter(iTrial).repTrialType =repmat(tType,1,nTimes);
    data_trials_conds_scatter(iTrial).behavior     =repmat(tType,1,nTimes);
    names   = {'Joint', 'Obs', 'Solo'}; 
    data_trials_conds_scatter(iTrial).repTrialName=names(data_trials_conds_scatter(iTrial).repTrialType);
end
par.plot_scatterGradient                = plot_scatterGradientParams();
par.plot_scatterGradient.hfig           = figure('visible',ifplot);
par.plot_scatterGradient.InField        = InField;

par.plot_scatterGradient.InGradient     = ['time' InField];
par.plot_scatterGradient.lats           = [1,2,3];              % directions to be plot     
%par.plot_scatterGradient.lats            = [2,3,1];              % directions to be plot     
% start gradient color for each class
par.plot_scatterGradient.cmaps           = linspecer(3);
par.plot_scatterGradient.cmapslight      = lightCmaps(par.plot_scatterGradient.cmaps);
% 
hfg.plot_scatterGradient    = plot_scatterGradient(data_trials_conds_scatter,par.plot_scatterGradient);

%% print hfg plots - if not ifplot, save fig results
if ~ifplot
    par.hfigPrint               = hfigPrintParams();
    par.hfigPrint.pdf_file      = [model_dir mfilename];
    par.hfigPrint.save_dir      = model_dir; 
    hfigPrint(hfg,par.hfigPrint)
end
