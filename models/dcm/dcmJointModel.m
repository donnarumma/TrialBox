function DCM = dcmJointModel(data_trials, pms, varargin)
% dcmJointModel  Fit one joint-monkey DCM from trial-aligned spectral data.
%
% DCM = dcmJointModel(data_trials)
% DCM = dcmJointModel(data_trials, pms)
% DCM = dcmJointModel(data_trials, pms, 'Name', value, ...)
%
% PURPOSE
%   Build the trial-by-condition design matrix used by SPM/DCM, create the
%   corresponding CSD model structure through csdModel, run the variational
%   inversion, and package the fitted session into one DCM struct.
%
% PARAMETERS
%   Use dcmJointModelParams for defaults. Parameters can be provided as a
%   pms struct, as name-value pairs, or both with precedence:
%   defaults < pms struct < name-value pairs.

executionTimer = tic;
fprintf('Function: %s ', mfilename);

if nargin < 2
    pms = struct();
end

pms = resolveFunctionParams(mfilename, pms, varargin{:});
inputFieldName = char(string(pms.InField));

assert(~isempty(data_trials), '%s received an empty data_trials array.', mfilename);
assert(isfield(data_trials, inputFieldName), ...
    '%s expected field %s in data_trials.', mfilename, inputFieldName);
assert(isfield(data_trials, 'klabels'), ...
    '%s expected field klabels in data_trials.', mfilename);

[~, allConditionLabels] = getJointMonkeysLabels(1:24);
selectedConditionIndices = pms.isdemo(:)';

if isempty(selectedConditionIndices)
    selectedConditionIndices = 1:numel(allConditionLabels);
end

if any(selectedConditionIndices < 1) || any(selectedConditionIndices > numel(allConditionLabels))
    error('%s received invalid condition indices in pms.isdemo.', mfilename);
end

nInputTrials = numel(data_trials);
xY = struct();
xY.y = cell(1, nInputTrials);
xU = struct();
xU.X = nan(nInputTrials, numel(allConditionLabels));

for trialIndex = 1:nInputTrials
    xY.y{trialIndex} = data_trials(trialIndex).(inputFieldName);
    xU.X(trialIndex, :) = data_trials(trialIndex).klabels;
end

nChannels = size(xY.y{1}, 2);
if ~isempty(pms.L)
    nSources = size(pms.L, 2);
else
    nSources = nChannels;
end

selectedTrialMask = any(xU.X(:, selectedConditionIndices), 2);
selectedConditionMask = sum(xU.X(selectedTrialMask, :), 1) > 0;

xU.X = xU.X(selectedTrialMask, selectedConditionMask);
xY.y = xY.y(selectedTrialMask);
selectedDataTrials = data_trials(selectedTrialMask);
selectedConditionLabels = allConditionLabels(selectedConditionMask);
xU.name = selectedConditionLabels;

if isempty(selectedDataTrials)
    error('%s did not find any trial matching the selected condition indices.', mfilename);
end

nConditions = size(xU.X, 2);

modelBuilderPms = pms;
modelBuilderPms.nConditions = nConditions;
modelBuilderPms.nSources = nSources;
modelBuilderPms.nChannels = nChannels;
modelBuilderPms.whichmodel = pms.whichmodel;
modelBuilderPms.customMode = local_resolve_custom_model(pms.custom_model);
modelBuilderPms.donlfp = pms.donlfp;
modelBuilderPms.Hz = pms.Hz;
modelBuilderPms.description = pms.description;

M_CSD = csdModel(modelBuilderPms);
spectralPredictionFunction = str2func(M_CSD.IS);

[posteriorMean, posteriorCovariance, logNoisePrecision, freeEnergy, ~, ~, ~, freeEnergyHistory] = ...
    DONNARUMMA_spm_nlsi_GN(M_CSD, xU, xY);

DCM = struct();
DCM.M = M_CSD;
DCM.Hz = M_CSD.Hz;
DCM.Ep = posteriorMean;
DCM.Cp = posteriorCovariance;
DCM.Ce = exp(-logNoisePrecision);

predictedChannelCsd = spectralPredictionFunction(posteriorMean, M_CSD, xU);
DCM.Hc = predictedChannelCsd;
DCM.Rc = spm_unvec(spm_vec(xY.y) - spm_vec(predictedChannelCsd), predictedChannelCsd);

sourceSpaceModel = local_build_source_space_model(M_CSD);
sourceSpacePosterior = posteriorMean;
sourceSpacePosterior.L = sourceSpaceModel.pE.L;
sourceSpacePosterior.b = sourceSpacePosterior.b - 32;
sourceSpacePosterior.c = sourceSpacePosterior.c - 32;

[predictedSourceCsd, ~, directedTransferFunctions] = spm_csd_mtf(sourceSpacePosterior, sourceSpaceModel, xU);
[crossCovarianceFunctions, peristimulusTime] = spm_csd2ccf(predictedSourceCsd, M_CSD.Hz);
[coherenceFunctions, frequencySpecificDelay] = spm_csd2coh(predictedSourceCsd, M_CSD.Hz);

DCM.dtf = directedTransferFunctions;
DCM.ccf = crossCovarianceFunctions;
DCM.coh = coherenceFunctions;
DCM.fsd = frequencySpecificDelay;
DCM.pst = peristimulusTime;
DCM.Hs = predictedSourceCsd;

priorVector = full(spm_vec(M_CSD.pE));
warningState = warning('off', 'SPM:negativeVariance');
DCM.Pp = spm_unvec(1 - spm_Ncdf(priorVector, abs(spm_vec(posteriorMean)), diag(posteriorCovariance)), posteriorMean);
DCM.Vp = spm_unvec(full(diag(posteriorCovariance)), posteriorMean);
warning(warningState);

DCM.F_History = freeEnergyHistory;
DCM.xU = xU;
DCM.xY = struct();
DCM.xY.csd = xY.y;
DCM.xY.lfp = local_collect_lfp_trials(selectedDataTrials);
DCM.xY.U = [];
DCM.xY.Hz = M_CSD.Hz;
DCM.xY.name = local_resolve_channel_names(nChannels);
DCM.Chamber = selectedDataTrials(1).Chamber;
DCM.F = freeEnergy;
DCM.selectedTrialMask = selectedTrialMask;
DCM.selectedConditionMask = selectedConditionMask;
DCM.par = pms;

if abs(DCM.F - DCM.F_History(end)) > max(eps(abs(DCM.F)), eps(abs(DCM.F_History(end))))
    warning('%s final Free Energy does not match the last history element.', mfilename);
end

executionTime = toc(executionTimer);
fprintf('| Time Elapsed: %.2f s\n', executionTime);
DCM.exectime = executionTime;
end


function customModelHandle = local_resolve_custom_model(customModelSpec)
if isempty(customModelSpec)
    customModelHandle = [];
    return;
end

if isa(customModelSpec, 'function_handle')
    customModelHandle = customModelSpec;
    return;
end

if isstring(customModelSpec)
    if strlength(customModelSpec) == 0
        customModelHandle = [];
        return;
    end
    customModelSpec = char(customModelSpec);
end

if ischar(customModelSpec)
    if isempty(strtrim(customModelSpec))
        customModelHandle = [];
        return;
    end
    customModelHandle = str2func(customModelSpec);
    return;
end

error('%s expected custom_model to be empty, char, string, or function handle.', mfilename);
end


function sourceSpaceModel = local_build_source_space_model(M_CSD)
if isfield(M_CSD, 'U')
    sourceSpaceModel = rmfield(M_CSD, 'U');
else
    sourceSpaceModel = M_CSD;
end

sourceSpaceModel.dipfit.type = 'LFP';
sourceSpaceModel.U = 1;
end


function lfpTrials = local_collect_lfp_trials(data_trials)
nTrials = numel(data_trials);
lfpTrials = cell(1, nTrials);

for trialIndex = 1:nTrials
    lfpTrials{trialIndex} = data_trials(trialIndex).LFP';
end
end


function channelNames = local_resolve_channel_names(nChannels)
if nChannels == 2
    channelNames = {'monkey S', 'monkey K'};
    return;
end

channelNames = cell(1, nChannels);
for channelIndex = 1:nChannels
    channelNames{channelIndex} = sprintf('channel %d', channelIndex);
end
end
