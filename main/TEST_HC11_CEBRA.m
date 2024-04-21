% function TEST_HC11_CEBRA
clear;

% par.irng                            = 10;                  % for reproducibility
% rng(par.irng);
%% Step 0. load raw data
hc11_dir                            = '~/DATA/HC11/';
behavior_filename                   = [hc11_dir 'behavior.hd5'];
neural_filename                     = [hc11_dir 'neural.hd5'];
group_field                         = 'data';
behavior_field                      = 'behavior';
manifold_field                      = 'manifold';
neural_field                        = 'neural';
% nChannels   x nTimes
neural                              = h5read(neural_filename,   ['/' group_field '/' neural_field]);
% nComponents x nTimes
behavior                            = h5read(behavior_filename, ['/' group_field '/' behavior_field]);
%% Step 1. Arrange Trials - prepare data in TrialBox format
trialNames                          = {'Left','Right'};         % label names
repTrialType                        = behavior(2,:)+1;          % left 1, right 2;
data_trials(1).behavior             = behavior;
data_trials(1).neural               = neural;
data_trials(1).repTrialType         = repTrialType;
data_trials(1).repTrialName         = trialNames(data_trials(1).repTrialType);
n_ms                                = 25; % sample time in ms
n_s                                 = n_ms/1000;
nTimes                              = size(neural,2);
data_trials(1).timeneural           = linspace(0,nTimes*n_s,nTimes);
%% Step 2. perform cebra model
strdate                             = datetime('now', 'Format', 'yyyyMMddHHmmss');
script_rundir                       = ['/tmp/cebra' sprintf('%s',strdate) '/'];
mkdir(script_rundir);
cebra_codes_dir                     = '~/tools/TrialBox/pywraps/';
par.cebraModel                      = cebraModelParams;
% path to python scripts
par.cebraModel.script_filename      = [cebra_codes_dir 'cebraModel.py'];
% path to hd5 files
par.cebraModel.script_rundir        = script_rundir;
par.cebraModel.model_filename       = [script_rundir 'cebra_model.pkl'];
% other parameteres
par.cebraModel.max_iterations       = 10000;
par.cebraModel.output_dimension     = 3;
disp(par.cebraModel);
[~,out.cebraModel]                  = cebraModel(data_trials,par.cebraModel);
%% Step 3. cebraEncode -> project data on the manifold
strdate                             = datetime('now', 'Format', 'yyyyMMddHHmmss');
script_rundir                       = ['/tmp/cebra' sprintf('%s',strdate) '/'];
mkdir(script_rundir);
par.cebraEncode                     = cebraEncodeParams();
par.cebraEncode.model_filename      = par.cebraModel.model_filename;
par.cebraEncode.script_rundir       = script_rundir;
par.cebraEncode.script_filename     = [cebra_codes_dir 'cebraEncode.py'];
disp(par.cebraEncode)
data_trials                         = cebraEncode(data_trials,par.cebraEncode);

%% plot_scatterGradient
gfield                              = 'gradients';
nTrials                             = length(data_trials);
for iTrial=1:nTrials
    data_trials(iTrial).(gfield)    = data_trials(iTrial).(par.cebraModel.behavior_field)(1,:);
end
par.plot_scatterGradient            = plot_scatterGradientParams();
par.plot_scatterGradient.InField    = par.cebraEncode.manifold_field;%'cebra';
par.plot_scatterGradient.InGradient = gfield;
par.plot_scatterGradient.lats       = [2,3,1];             % directions to be plot     
par.plot_scatterGradient.reverse    = [false,false,true];  % reverse directions axis
% start gradient color for each class
par.plot_scatterGradient.cmapslight = [[1.0,0.0,1.0]; ...  % left  start from magenta
                                       [1.0,1.0,0.4]];     % right start from yellow
% end gradient color for each class
par.plot_scatterGradient.cmaps      = [[0.0,1.0,1.0]; ...  % left  end in cyan
                                            [0.0,0.5,0.4]];     % right end in green
par.plot_scatterGradinnt.label      = 'm';
hfg.plot_scatterGradient            = plot_scatterGradient(data_trials,par.plot_scatterGradient);
title('cebra','interpreter','none')