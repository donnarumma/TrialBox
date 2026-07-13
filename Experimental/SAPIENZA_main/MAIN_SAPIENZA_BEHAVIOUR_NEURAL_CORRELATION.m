%% MAIN_SAPIENZA_BEHAVIOUR_NEURAL_CORRELATION
% Relate RAW behavioural readouts to session-level DCM neural readouts.
%
% The analysis builds two outputs:
%   1) a long condition table: chamber/session/model/direction/condition
%   2) a joint-action table: chamber/session/model/direction with behaviour
%      asymmetries and DCM crossed-influence readouts.
%
% Positive SAdvantage_* means monkey S performs better than monkey K for that
% metric polarity. Positive NeuralBalanceSToKMinusKToS means the DCM readout
% favours S to K over K to S.

clear;
close all;

currentScriptDir = fileparts(mfilename('fullpath'));
trialBoxRoot = fileparts(fileparts(currentScriptDir));
neuroConnStatRoot = '/home/donnarumma/OneDrive/tools/NeuroConnStat';
if ~exist(neuroConnStatRoot, 'dir')
    neuroConnStatRoot = '/media/Data/OneDrive/tools/NeuroConnStat';
end
spmRoot = '/media/Data/OneDrive/tools/spm12';
if ~exist(spmRoot, 'dir')
    spmRoot = '/home/donnarumma/OneDrive/tools/spm12';
end

addpath(fullfile(trialBoxRoot, 'Experimental', 'Experimental_utility', 'Utility_Behaviour'));
addpath(fullfile(trialBoxRoot, 'utilities', 'pms'));
if exist(spmRoot, 'dir')
    addpath(spmRoot);
    addpath(fullfile(spmRoot, 'toolbox', 'dcm_meeg'));
end
addpath(genpath(fullfile(neuroConnStatRoot, 'DCMspace')));

analysisPms = struct();
analysisPms.behaviourRoot = '/TESTS/SAPIENZA/BEHAVIOUR';
analysisPms.dcmRoot = '/TESTS/SAPIENZA/DCM';
analysisPms.outputRoot = fullfile(analysisPms.behaviourRoot, 'NeuralCorrelation');
analysisPms.chamberNames = {'Frontal', 'Parietal'};
analysisPms.modelNames = {'M2', 'M4', 'M7', 'M8', 'M9', 'M11', 'M12'};
analysisPms.primaryModelNames = {'M2', 'M4', 'M7', 'M8'};
analysisPms.baselineModelNames = {'M9', 'M11', 'M12'};
analysisPms.behaviourMetrics = {'RT', 'PV', 'PVT', 'AE_abs', 'CMT', 'EC', 'ExT'};
analysisPms.keyNeuralTargets = {'NeuralBalanceSToKMinusKToS', 'NeuralJointGainBalance', ...
    'JointWeightSToK', 'JointWeightKToS', 'JointGainSToK', 'JointGainKToS'};
analysisPms.conditionLabels = {'Act S / Obs K', 'Obs S / Act K', 'Act S / Act K'};
analysisPms.correlationType = 'Spearman';
analysisPms.minimumSamplesForCorrelation = 5;
analysisPms.directionalPlotPairs = { ...
    'SAdvantage_RT', 'NeuralBalanceSToKMinusKToS'; ...
    'SAdvantage_RT', 'NeuralJointGainBalance'; ...
    'SAdvantage_AE_abs', 'NeuralJointGainBalance'; ...
    'SAdvantage_EC', 'NeuralJointGainBalance'; ...
    'SAdvantage_ExT', 'NeuralJointGainBalance'; ...
    'KJointBenefit_RT', 'NeuralJointGainBalance'; ...
    'KJointBenefit_RT', 'JointWeightSToK'; ...
    'SJointBenefit_RT', 'JointWeightKToS'; ...
    'KJointBenefit_RT', 'JointGainSToK'; ...
    'SJointBenefit_RT', 'JointGainKToS'; ...
    'KJointBenefit_CMT', 'JointWeightSToK'; ...
    'KJointBenefit_EC', 'JointGainSToK'};
analysisPms.focusModelRows = { ...
    'Frontal',  'M7', 'Frontal PEB-best trajectory model'; ...
    'Frontal',  'M8', 'Frontal model-space competitor'; ...
    'Parietal', 'M2', 'Parietal PEB/BMS-best direct-coupling model'; ...
    'Parietal', 'M7', 'Parietal trajectory comparator'};
analysisPms.defaultLFPMethod = 'median';
analysisPms.lfpBestRows = { ...
    'Frontal',  'SK030', 'Maximum'; ...
    'Frontal',  'SK009', 'Maximum'; ...
    'Frontal',  'SK029', 'Maximum'; ...
    'Frontal',  'SK028', 'Maximum'; ...
    'Frontal',  'SK026', 'mean'; ...
    'Parietal', 'SK054', 'Maximum'; ...
    'Parietal', 'SK064', 'Maximum'};
analysisPms.comparePms = struct( ...
    'CondPairs', [1 2; 1 3; 2 3], ...
    'WeightScale', 1, ...
    'Correction', 'none', ...
    'Alpha', 0.15, ...
    'SummaryOptions', struct( ...
        'includeBaselineInA', true, ...
        'includeBaselineInB', true, ...
        'includeD', true));

if ~exist(analysisPms.outputRoot, 'dir')
    mkdir(analysisPms.outputRoot);
end

conditionRows = cell(0, numel(local_condition_table_columns(analysisPms)));
jointRows = cell(0, numel(local_joint_table_columns(analysisPms)));
missingRows = cell(0, 5);

for chamberCounter = 1:numel(analysisPms.chamberNames)
    chamberName = analysisPms.chamberNames{chamberCounter};
    sessionNames = local_find_behaviour_sessions(analysisPms.behaviourRoot, chamberName);

    fprintf('Behaviour-neural correlation: %s (%d behaviour sessions found)\n', ...
        chamberName, numel(sessionNames));

    for sessionCounter = 1:numel(sessionNames)
        sessionName = sessionNames{sessionCounter};
        behaviourPath = local_find_latest_behaviour_file( ...
            analysisPms.behaviourRoot, chamberName, sessionName);

        if isempty(behaviourPath)
            missingRows(end + 1, :) = {chamberName, sessionName, '', '', 'missing behaviour'}; %#ok<SAGROW>
            continue;
        end

        loadedBehaviour = load(behaviourPath, 'Data');
        behaviourData = loadedBehaviour.Data.Behav;

        for modelCounter = 1:numel(analysisPms.modelNames)
            modelName = analysisPms.modelNames{modelCounter};
            expectedLFPMethod = local_expected_lfp_method(analysisPms, chamberName, sessionName);
            dcmPath = local_find_lfp_best_dcm( ...
                analysisPms.dcmRoot, modelName, chamberName, sessionName, expectedLFPMethod);

            if isempty(dcmPath)
                missingRows(end + 1, :) = { ...
                    chamberName, sessionName, modelName, '', ...
                    sprintf('missing DCM for LFP-best method %s', expectedLFPMethod)}; %#ok<SAGROW>
                continue;
            end

            try
                loadedDcm = load(dcmPath, 'DCM');
                dcmStruct = loadedDcm.DCM;
                neuralReadout = local_extract_neural_readout(dcmStruct, modelName, analysisPms);
                freeEnergy = local_get_free_energy(dcmStruct);
                lfpMethod = local_parse_lfp_method(dcmPath);

                [newConditionRows, newJointRows] = local_build_session_rows( ...
                    behaviourData, neuralReadout, chamberName, sessionName, modelName, ...
                    lfpMethod, freeEnergy, behaviourPath, dcmPath, analysisPms);

                conditionRows = [conditionRows; newConditionRows]; %#ok<AGROW>
                jointRows = [jointRows; newJointRows]; %#ok<AGROW>
            catch analysisError
                missingRows(end + 1, :) = {chamberName, sessionName, modelName, dcmPath, analysisError.message}; %#ok<SAGROW>
            end
        end
    end
end

conditionTable = cell2table(conditionRows, 'VariableNames', local_condition_table_columns(analysisPms));
jointTable = cell2table(jointRows, 'VariableNames', local_joint_table_columns(analysisPms));
missingTable = cell2table(missingRows, ...
    'VariableNames', {'Chamber', 'Session', 'Model', 'Path', 'Reason'});
correlationTable = local_build_correlation_table(jointTable, analysisPms);
directionalCorrelationTable = local_build_directional_correlation_table(jointTable, analysisPms);
directionalAllTargetCorrelationTable = local_build_directional_all_target_correlation_table( ...
    jointTable, analysisPms);
directionalJointGainBalanceCorrelationTable = local_extract_target_correlation_table( ...
    directionalAllTargetCorrelationTable, 'NeuralJointGainBalance');
behaviourSummaryTable = local_build_behaviour_summary(jointTable, analysisPms);
jointGainBalanceCorrelationTable = local_extract_target_correlation_table( ...
    correlationTable, 'NeuralJointGainBalance');
primaryCorrelationTable = local_filter_table_models(correlationTable, analysisPms.primaryModelNames);
primaryDirectionalAllTargetCorrelationTable = local_filter_table_models( ...
    directionalAllTargetCorrelationTable, analysisPms.primaryModelNames);

writetable(conditionTable, fullfile(analysisPms.outputRoot, 'BehaviourNeural_conditionRows.csv'));
writetable(jointTable, fullfile(analysisPms.outputRoot, 'BehaviourNeural_jointRows.csv'));
writetable(correlationTable, fullfile(analysisPms.outputRoot, 'BehaviourNeural_correlations.csv'));
writetable(jointGainBalanceCorrelationTable, ...
    fullfile(analysisPms.outputRoot, 'BehaviourNeural_jointGainBalanceCorrelations.csv'));
writetable(directionalCorrelationTable, fullfile(analysisPms.outputRoot, 'BehaviourNeural_directionalCorrelations.csv'));
writetable(directionalAllTargetCorrelationTable, ...
    fullfile(analysisPms.outputRoot, 'BehaviourNeural_directionalAllTargetCorrelations.csv'));
writetable(directionalJointGainBalanceCorrelationTable, ...
    fullfile(analysisPms.outputRoot, 'BehaviourNeural_directionalJointGainBalanceCorrelations.csv'));
writetable(primaryCorrelationTable, ...
    fullfile(analysisPms.outputRoot, 'BehaviourNeural_primaryModelCorrelations.csv'));
writetable(primaryDirectionalAllTargetCorrelationTable, ...
    fullfile(analysisPms.outputRoot, 'BehaviourNeural_directionalPrimaryModelCorrelations.csv'));
writetable(behaviourSummaryTable, fullfile(analysisPms.outputRoot, 'BehaviourNeural_behaviourSummary.csv'));
writetable(missingTable, fullfile(analysisPms.outputRoot, 'BehaviourNeural_missingInputs.csv'));

save(fullfile(analysisPms.outputRoot, 'BehaviourNeural_correlation.mat'), ...
    'analysisPms', 'conditionTable', 'jointTable', 'correlationTable', ...
    'jointGainBalanceCorrelationTable', 'directionalCorrelationTable', ...
    'directionalAllTargetCorrelationTable', ...
    'directionalJointGainBalanceCorrelationTable', 'behaviourSummaryTable', ...
    'primaryCorrelationTable', 'primaryDirectionalAllTargetCorrelationTable', ...
    'missingTable', '-v7.3');

local_plot_joint_gain_balance_metric_grid(jointGainBalanceCorrelationTable, analysisPms);
local_plot_key_correlations(jointTable, correlationTable, analysisPms);
local_plot_focus_model_correlations(jointTable, analysisPms);
local_plot_directional_correlations(directionalCorrelationTable, analysisPms);
local_plot_m7_rt_directional_stability( ...
    directionalAllTargetCorrelationTable, correlationTable, analysisPms);
local_plot_m7_nonrt_directional_bars( ...
    directionalAllTargetCorrelationTable, correlationTable, analysisPms);
local_save_directional_mean_abs_report( ...
    directionalAllTargetCorrelationTable, correlationTable, analysisPms);
local_plot_best_receiver_directional_bars( ...
    directionalAllTargetCorrelationTable, correlationTable, analysisPms);
local_plot_k_receiver_directional_bars( ...
    directionalAllTargetCorrelationTable, correlationTable, analysisPms);
local_plot_receiver_benefit_correlations(jointTable, analysisPms);
local_save_top_correlation_report(correlationTable, directionalAllTargetCorrelationTable, ...
    jointGainBalanceCorrelationTable, analysisPms);

fprintf('Behaviour-neural analysis saved in %s\n', analysisPms.outputRoot);


function sessionNames = local_find_behaviour_sessions(behaviourRoot, chamberName)
sessionDirInfo = dir(fullfile(behaviourRoot, chamberName, 'SK*'));
sessionDirInfo = sessionDirInfo([sessionDirInfo.isdir]);
sessionNames = {sessionDirInfo.name};
sessionNames = sort(sessionNames);
end


function behaviourPath = local_find_latest_behaviour_file(behaviourRoot, chamberName, sessionName)
filePattern = sprintf('Behaviour_RAW_%s%s_D*_C*.mat', chamberName, sessionName);
fileInfo = dir(fullfile(behaviourRoot, chamberName, sessionName, filePattern));
behaviourPath = local_latest_file_path(fileInfo);
end


function expectedMethod = local_expected_lfp_method(analysisPms, chamberName, sessionName)
expectedMethod = analysisPms.defaultLFPMethod;

for lfpRowCounter = 1:size(analysisPms.lfpBestRows, 1)
    rowChamber = analysisPms.lfpBestRows{lfpRowCounter, 1};
    rowSession = analysisPms.lfpBestRows{lfpRowCounter, 2};
    rowMethod = analysisPms.lfpBestRows{lfpRowCounter, 3};

    if strcmpi(rowChamber, chamberName) && strcmpi(rowSession, sessionName)
        expectedMethod = rowMethod;
        return;
    end
end
end


function dcmPath = local_find_lfp_best_dcm(dcmRoot, modelName, chamberName, sessionName, expectedLFPMethod)
filePattern = sprintf('LFP_*_%s_%s%s_D*_C*_*', modelName, chamberName, sessionName);
fileInfo = dir(fullfile(dcmRoot, modelName, chamberName, sessionName, filePattern));
if isempty(fileInfo)
    dcmPath = '';
    return;
end

fileMethod = strings(numel(fileInfo), 1);
for fileCounter = 1:numel(fileInfo)
    fileMethod(fileCounter) = local_parse_lfp_method(fileInfo(fileCounter).name);
end

methodMask = strcmpi(fileMethod, expectedLFPMethod);
dcmPath = local_latest_file_path(fileInfo(methodMask));
end


function latestPath = local_latest_file_path(fileInfo)
latestPath = '';
if isempty(fileInfo)
    return;
end

[~, sortOrder] = sort([fileInfo.datenum], 'descend');
latestFile = fileInfo(sortOrder(1));
latestPath = fullfile(latestFile.folder, latestFile.name);
end


function neuralReadout = local_extract_neural_readout(dcmStruct, modelName, analysisPms)
switch upper(modelName)
    case 'M2'
        comparisonResult = dcm_compare(dcmStruct, analysisPms.comparePms);
        summaryData = dcm_extract_summary( ...
            comparisonResult, ...
            'Correction', analysisPms.comparePms.Correction, ...
            'Alpha', analysisPms.comparePms.Alpha);

        directionCount = numel(summaryData.meta.direction_ids);
        conditionCount = height(summaryData.Conn_weights_by_direction{1, 1});
        weightSToK = nan(directionCount, conditionCount);
        weightKToS = nan(directionCount, conditionCount);

        for directionCounter = 1:directionCount
            weightSToK(directionCounter, :) = summaryData.Conn_weights_by_direction{directionCounter, 1}.Weight(:)';
            weightKToS(directionCounter, :) = summaryData.Conn_weights_by_direction{directionCounter, 2}.Weight(:)';
        end

        neuralReadout = struct();
        neuralReadout.weightSToK = weightSToK;
        neuralReadout.weightKToS = weightKToS;
        neuralReadout.directionLabels = summaryData.meta.direction_labels;
        neuralReadout.weightSToKLabel = 'W1->2';
        neuralReadout.weightKToSLabel = 'W2->1';

    otherwise
        pathwayData = dcm_extract_summary_traj2D(dcmStruct);
        neuralReadout = struct();
        neuralReadout.weightSToK = pathwayData.pathwayStrength12;
        neuralReadout.weightKToS = pathwayData.pathwayStrength21;
        neuralReadout.directionLabels = pathwayData.directionLabels;
        neuralReadout.weightSToKLabel = 'TrajS->AreaK';
        neuralReadout.weightKToSLabel = 'TrajK->AreaS';
end
end


function freeEnergy = local_get_free_energy(dcmStruct)
freeEnergy = NaN;
if isfield(dcmStruct, 'F_History') && ~isempty(dcmStruct.F_History)
    freeEnergy = dcmStruct.F_History(end);
elseif isfield(dcmStruct, 'F') && ~isempty(dcmStruct.F)
    freeEnergy = dcmStruct.F;
end
end


function lfpMethod = local_parse_lfp_method(dcmPath)
lfpMethod = "";
methodTokens = regexp(dcmPath, 'LFP_([^_/]+)_', 'tokens', 'once');
if ~isempty(methodTokens)
    lfpMethod = string(methodTokens{1});
end
end


function [conditionRows, jointRows] = local_build_session_rows( ...
    behaviourData, neuralReadout, chamberName, sessionName, modelName, ...
    lfpMethod, freeEnergy, behaviourPath, dcmPath, analysisPms)

conditionRows = {};
jointRows = {};

for directionCounter = 1:8
    sourceActSToK = neuralReadout.weightSToK(directionCounter, 1);
    sourceActKToS = neuralReadout.weightKToS(directionCounter, 2);
    jointSToK = neuralReadout.weightSToK(directionCounter, 3);
    jointKToS = neuralReadout.weightKToS(directionCounter, 3);
    jointGainSToK = jointSToK - sourceActSToK;
    jointGainKToS = jointKToS - sourceActKToS;
    neuralBalance = jointSToK - jointKToS;
    neuralJointGainBalance = jointGainSToK - jointGainKToS;

    for conditionCode = 1:3
        behaviourIndex = (conditionCode - 1) * 8 + directionCounter;
        rowValues = { ...
            chamberName, sessionName, modelName, string(lfpMethod), freeEnergy, ...
            directionCounter, conditionCode, analysisPms.conditionLabels{conditionCode}, ...
            neuralReadout.weightSToK(directionCounter, conditionCode), ...
            neuralReadout.weightKToS(directionCounter, conditionCode), ...
            neuralReadout.weightSToK(directionCounter, conditionCode) - neuralReadout.weightKToS(directionCounter, conditionCode), ...
            jointGainSToK, jointGainKToS, neuralJointGainBalance, ...
            string(neuralReadout.weightSToKLabel), string(neuralReadout.weightKToSLabel), ...
            string(behaviourPath), string(dcmPath)};

        for metricCounter = 1:numel(analysisPms.behaviourMetrics)
            metricName = analysisPms.behaviourMetrics{metricCounter};
            rowValues = [rowValues, { ...
                behaviourData.S(behaviourIndex).(metricName), ...
                behaviourData.K(behaviourIndex).(metricName)}]; %#ok<AGROW>
        end

        conditionRows(end + 1, :) = rowValues; %#ok<AGROW>
    end

    jointValues = { ...
        chamberName, sessionName, modelName, string(lfpMethod), freeEnergy, directionCounter, ...
        jointSToK, jointKToS, neuralBalance, jointGainSToK, jointGainKToS, ...
        neuralJointGainBalance, string(neuralReadout.weightSToKLabel), ...
        string(neuralReadout.weightKToSLabel), string(behaviourPath), string(dcmPath)};

    for metricCounter = 1:numel(analysisPms.behaviourMetrics)
        metricName = analysisPms.behaviourMetrics{metricCounter};
        sActValue = behaviourData.S(directionCounter).(metricName);
        kActValue = behaviourData.K(8 + directionCounter).(metricName);
        sJointValue = behaviourData.S(16 + directionCounter).(metricName);
        kJointValue = behaviourData.K(16 + directionCounter).(metricName);
        sAdvantage = local_s_advantage(metricName, sJointValue, kJointValue);
        sJointBenefit = local_joint_benefit(metricName, sActValue, sJointValue);
        kJointBenefit = local_joint_benefit(metricName, kActValue, kJointValue);

        jointValues = [jointValues, { ...
            sActValue, kActValue, sJointValue, kJointValue, ...
            sAdvantage, sJointBenefit, kJointBenefit}]; %#ok<AGROW>
    end

    jointRows(end + 1, :) = jointValues; %#ok<AGROW>
end
end


function sAdvantage = local_s_advantage(metricName, sJointValue, kJointValue)
if local_is_lower_better(metricName)
    sAdvantage = kJointValue - sJointValue;
else
    sAdvantage = sJointValue - kJointValue;
end
end


function jointBenefit = local_joint_benefit(metricName, actValue, jointValue)
if local_is_lower_better(metricName)
    jointBenefit = actValue - jointValue;
else
    jointBenefit = jointValue - actValue;
end
end


function isLowerBetter = local_is_lower_better(metricName)
isLowerBetter = ismember(metricName, {'RT', 'PVT', 'AE_abs', 'CMT', 'ExT'});
end


function variableNames = local_condition_table_columns(analysisPms)
variableNames = {'Chamber', 'Session', 'Model', 'LFPMethod', 'FreeEnergy', ...
    'Direction', 'ConditionCode', 'ConditionName', ...
    'WeightSToK', 'WeightKToS', 'WeightBalanceSToKMinusKToS', ...
    'JointGainSToK', 'JointGainKToS', 'NeuralJointGainBalance', ...
    'WeightSToKLabel', 'WeightKToSLabel', 'BehaviourPath', 'DCMPath'};

for metricCounter = 1:numel(analysisPms.behaviourMetrics)
    metricName = analysisPms.behaviourMetrics{metricCounter};
    variableNames = [variableNames, { ...
        sprintf('S_%s', metricName), ...
        sprintf('K_%s', metricName)}]; %#ok<AGROW>
end
end


function variableNames = local_joint_table_columns(analysisPms)
variableNames = {'Chamber', 'Session', 'Model', 'LFPMethod', 'FreeEnergy', ...
    'Direction', 'JointWeightSToK', 'JointWeightKToS', ...
    'NeuralBalanceSToKMinusKToS', 'JointGainSToK', 'JointGainKToS', ...
    'NeuralJointGainBalance', 'WeightSToKLabel', 'WeightKToSLabel', ...
    'BehaviourPath', 'DCMPath'};

for metricCounter = 1:numel(analysisPms.behaviourMetrics)
    metricName = analysisPms.behaviourMetrics{metricCounter};
    variableNames = [variableNames, { ...
        sprintf('SAct_%s', metricName), ...
        sprintf('KAct_%s', metricName), ...
        sprintf('SJoint_%s', metricName), ...
        sprintf('KJoint_%s', metricName), ...
        sprintf('SAdvantage_%s', metricName), ...
        sprintf('SJointBenefit_%s', metricName), ...
        sprintf('KJointBenefit_%s', metricName)}]; %#ok<AGROW>
end
end


function correlationTable = local_build_correlation_table(jointTable, analysisPms)
if isempty(jointTable)
    correlationTable = table();
    return;
end

targetNames = analysisPms.keyNeuralTargets;
predictorNames = {};
for metricCounter = 1:numel(analysisPms.behaviourMetrics)
    metricName = analysisPms.behaviourMetrics{metricCounter};
    predictorNames = [predictorNames, { ...
        sprintf('SAdvantage_%s', metricName), ...
        sprintf('SJointBenefit_%s', metricName), ...
        sprintf('KJointBenefit_%s', metricName)}]; %#ok<AGROW>
end

correlationRows = {};
chamberNames = unique(jointTable.Chamber, 'stable');
modelNames = unique(jointTable.Model, 'stable');

for chamberCounter = 1:numel(chamberNames)
    chamberName = chamberNames{chamberCounter};
    for modelCounter = 1:numel(modelNames)
        modelName = modelNames{modelCounter};
        groupMask = strcmp(jointTable.Chamber, chamberName) & strcmp(jointTable.Model, modelName);

        for predictorCounter = 1:numel(predictorNames)
            predictorName = predictorNames{predictorCounter};
            predictorValues = jointTable.(predictorName)(groupMask);

            for targetCounter = 1:numel(targetNames)
                targetName = targetNames{targetCounter};
                targetValues = jointTable.(targetName)(groupMask);
                validSamples = isfinite(predictorValues) & isfinite(targetValues);
                sampleCount = sum(validSamples);
                rhoValue = NaN;
                probabilityValue = NaN;

                if sampleCount >= analysisPms.minimumSamplesForCorrelation
                    [rhoValue, probabilityValue] = corr( ...
                        predictorValues(validSamples), ...
                        targetValues(validSamples), ...
                        'Type', analysisPms.correlationType, ...
                        'Rows', 'complete');
                end

                correlationRows(end + 1, :) = { ...
                    chamberName, modelName, predictorName, targetName, ...
                    analysisPms.correlationType, sampleCount, rhoValue, probabilityValue}; %#ok<AGROW>
            end
        end
    end
end

correlationTable = cell2table(correlationRows, ...
    'VariableNames', {'Chamber', 'Model', 'Predictor', 'Target', ...
    'CorrelationType', 'N', 'Rho', 'PValue'});
end


function targetCorrelationTable = local_extract_target_correlation_table(correlationTable, targetName)
if isempty(correlationTable)
    targetCorrelationTable = table();
    return;
end

targetMask = strcmp(correlationTable.Target, targetName);
targetCorrelationTable = correlationTable(targetMask, :);
if isempty(targetCorrelationTable)
    return;
end

sortProbability = targetCorrelationTable.PValue;
sortProbability(~isfinite(sortProbability)) = Inf;
targetCorrelationTable.SortProbability = sortProbability;
targetCorrelationTable.SortAbsRho = -abs(targetCorrelationTable.Rho);
targetCorrelationTable = sortrows(targetCorrelationTable, {'SortProbability', 'SortAbsRho'});
targetCorrelationTable.SortProbability = [];
targetCorrelationTable.SortAbsRho = [];
end


function directionalCorrelationTable = local_build_directional_correlation_table(jointTable, analysisPms)
if isempty(jointTable)
    directionalCorrelationTable = table();
    return;
end

predictorNames = unique(analysisPms.directionalPlotPairs(:, 1), 'stable');
targetNames = unique(analysisPms.directionalPlotPairs(:, 2), 'stable');

directionalRows = {};
chamberNames = unique(jointTable.Chamber, 'stable');
modelNames = unique(jointTable.Model, 'stable');
directionIds = unique(jointTable.Direction, 'stable');

for chamberCounter = 1:numel(chamberNames)
    chamberName = chamberNames{chamberCounter};
    for modelCounter = 1:numel(modelNames)
        modelName = modelNames{modelCounter};
        for directionCounter = 1:numel(directionIds)
            directionId = directionIds(directionCounter);
            groupMask = strcmp(jointTable.Chamber, chamberName) & ...
                strcmp(jointTable.Model, modelName) & ...
                jointTable.Direction == directionId;

            for predictorCounter = 1:numel(predictorNames)
                predictorName = predictorNames{predictorCounter};
                predictorValues = jointTable.(predictorName)(groupMask);

                for targetCounter = 1:numel(targetNames)
                    targetName = targetNames{targetCounter};
                    targetValues = jointTable.(targetName)(groupMask);
                    validSamples = isfinite(predictorValues) & isfinite(targetValues);
                    sampleCount = sum(validSamples);
                    rhoValue = NaN;
                    probabilityValue = NaN;

                    if sampleCount >= analysisPms.minimumSamplesForCorrelation
                        [rhoValue, probabilityValue] = corr( ...
                            predictorValues(validSamples), ...
                            targetValues(validSamples), ...
                            'Type', analysisPms.correlationType, ...
                            'Rows', 'complete');
                    end

                    directionalRows(end + 1, :) = { ...
                        chamberName, modelName, directionId, predictorName, targetName, ...
                        analysisPms.correlationType, sampleCount, rhoValue, probabilityValue}; %#ok<AGROW>
                end
            end
        end
    end
end

directionalCorrelationTable = cell2table(directionalRows, ...
    'VariableNames', {'Chamber', 'Model', 'Direction', 'Predictor', 'Target', ...
    'CorrelationType', 'N', 'Rho', 'PValue'});
end


function directionalTargetCorrelationTable = local_build_directional_target_correlation_table( ...
    jointTable, analysisPms, targetName)
if isempty(jointTable)
    directionalTargetCorrelationTable = table();
    return;
end

predictorNames = {};
for metricCounter = 1:numel(analysisPms.behaviourMetrics)
    metricName = analysisPms.behaviourMetrics{metricCounter};
    predictorNames = [predictorNames, { ...
        sprintf('SAdvantage_%s', metricName), ...
        sprintf('SJointBenefit_%s', metricName), ...
        sprintf('KJointBenefit_%s', metricName)}]; %#ok<AGROW>
end

directionalRows = {};
chamberNames = unique(jointTable.Chamber, 'stable');
modelNames = unique(jointTable.Model, 'stable');
directionIds = unique(jointTable.Direction, 'stable');

for chamberCounter = 1:numel(chamberNames)
    chamberName = chamberNames{chamberCounter};
    for modelCounter = 1:numel(modelNames)
        modelName = modelNames{modelCounter};
        for directionCounter = 1:numel(directionIds)
            directionId = directionIds(directionCounter);
            groupMask = strcmp(jointTable.Chamber, chamberName) & ...
                strcmp(jointTable.Model, modelName) & ...
                jointTable.Direction == directionId;
            targetValues = jointTable.(targetName)(groupMask);

            for predictorCounter = 1:numel(predictorNames)
                predictorName = predictorNames{predictorCounter};
                predictorValues = jointTable.(predictorName)(groupMask);
                validSamples = isfinite(predictorValues) & isfinite(targetValues);
                sampleCount = sum(validSamples);
                rhoValue = NaN;
                probabilityValue = NaN;

                if sampleCount >= analysisPms.minimumSamplesForCorrelation
                    [rhoValue, probabilityValue] = corr( ...
                        predictorValues(validSamples), ...
                        targetValues(validSamples), ...
                        'Type', analysisPms.correlationType, ...
                        'Rows', 'complete');
                end

                directionalRows(end + 1, :) = { ...
                    chamberName, modelName, directionId, predictorName, targetName, ...
                    analysisPms.correlationType, sampleCount, rhoValue, probabilityValue}; %#ok<AGROW>
            end
        end
    end
end

directionalTargetCorrelationTable = cell2table(directionalRows, ...
    'VariableNames', {'Chamber', 'Model', 'Direction', 'Predictor', 'Target', ...
    'CorrelationType', 'N', 'Rho', 'PValue'});
end


function directionalAllTargetCorrelationTable = local_build_directional_all_target_correlation_table( ...
    jointTable, analysisPms)
if isempty(jointTable)
    directionalAllTargetCorrelationTable = table();
    return;
end

targetTables = cell(numel(analysisPms.keyNeuralTargets), 1);
for targetCounter = 1:numel(analysisPms.keyNeuralTargets)
    targetName = analysisPms.keyNeuralTargets{targetCounter};
    targetTables{targetCounter} = local_build_directional_target_correlation_table( ...
        jointTable, analysisPms, targetName);
end

directionalAllTargetCorrelationTable = vertcat(targetTables{:});
end


function behaviourSummaryTable = local_build_behaviour_summary(jointTable, analysisPms)
summaryRows = {};
if isempty(jointTable)
    behaviourSummaryTable = table();
    return;
end

chamberNames = unique(jointTable.Chamber, 'stable');

for chamberCounter = 1:numel(chamberNames)
    chamberName = chamberNames{chamberCounter};
    chamberMask = strcmp(jointTable.Chamber, chamberName);
    sessionDirectionMask = chamberMask & strcmp(jointTable.Model, jointTable.Model(find(chamberMask, 1, 'first')));

    for metricCounter = 1:numel(analysisPms.behaviourMetrics)
        metricName = analysisPms.behaviourMetrics{metricCounter};
        sAdvantageName = sprintf('SAdvantage_%s', metricName);
        sBenefitName = sprintf('SJointBenefit_%s', metricName);
        kBenefitName = sprintf('KJointBenefit_%s', metricName);

        summaryRows(end + 1, :) = { ...
            chamberName, metricName, ...
            mean(jointTable.(sAdvantageName)(sessionDirectionMask), 'omitnan'), ...
            std(jointTable.(sAdvantageName)(sessionDirectionMask), 'omitnan'), ...
            mean(jointTable.(sBenefitName)(sessionDirectionMask), 'omitnan'), ...
            std(jointTable.(sBenefitName)(sessionDirectionMask), 'omitnan'), ...
            mean(jointTable.(kBenefitName)(sessionDirectionMask), 'omitnan'), ...
            std(jointTable.(kBenefitName)(sessionDirectionMask), 'omitnan')}; %#ok<AGROW>
    end
end

behaviourSummaryTable = cell2table(summaryRows, ...
    'VariableNames', {'Chamber', 'Metric', ...
    'MeanSAdvantage', 'SdSAdvantage', ...
    'MeanSJointBenefit', 'SdSJointBenefit', ...
    'MeanKJointBenefit', 'SdKJointBenefit'});
end


function local_plot_joint_gain_balance_metric_grid(jointGainBalanceCorrelationTable, analysisPms)
if isempty(jointGainBalanceCorrelationTable)
    return;
end

plotDir = fullfile(analysisPms.outputRoot, 'plots');
if ~exist(plotDir, 'dir')
    mkdir(plotDir);
end

predictorPrefixes = {'SAdvantage', 'SJointBenefit', 'KJointBenefit'};
predictorLabels = {'S advantage', 'S joint benefit', 'K joint benefit'};
metricNames = analysisPms.behaviourMetrics;
modelNames = unique(jointGainBalanceCorrelationTable.Model, 'stable');

maxAbsRho = max(abs(jointGainBalanceCorrelationTable.Rho), [], 'omitnan');
if ~isfinite(maxAbsRho) || maxAbsRho == 0
    maxAbsRho = 0.1;
end
colorLimit = min(1, max(0.1, maxAbsRho));

hFigure = figure('Color', 'w', 'Visible', 'off', ...
    'Name', 'Behaviour metrics vs NeuralJointGainBalance', ...
    'Position', [100 100 1360 900]);
tiledlayout(numel(predictorPrefixes), numel(analysisPms.chamberNames), ...
    'TileSpacing', 'compact', 'Padding', 'compact');

for predictorCounter = 1:numel(predictorPrefixes)
    predictorPrefix = predictorPrefixes{predictorCounter};
    predictorLabel = predictorLabels{predictorCounter};

    for chamberCounter = 1:numel(analysisPms.chamberNames)
        chamberName = analysisPms.chamberNames{chamberCounter};
        correlationMatrix = nan(numel(modelNames), numel(metricNames));
        probabilityMatrix = nan(numel(modelNames), numel(metricNames));

        for modelCounter = 1:numel(modelNames)
            modelName = modelNames{modelCounter};
            for metricCounter = 1:numel(metricNames)
                metricName = metricNames{metricCounter};
                predictorName = sprintf('%s_%s', predictorPrefix, metricName);
                valueMask = strcmp(jointGainBalanceCorrelationTable.Chamber, chamberName) & ...
                    strcmp(jointGainBalanceCorrelationTable.Model, modelName) & ...
                    strcmp(jointGainBalanceCorrelationTable.Predictor, predictorName);

                if any(valueMask)
                    rowIndex = find(valueMask, 1, 'first');
                    correlationMatrix(modelCounter, metricCounter) = jointGainBalanceCorrelationTable.Rho(rowIndex);
                    probabilityMatrix(modelCounter, metricCounter) = jointGainBalanceCorrelationTable.PValue(rowIndex);
                end
            end
        end

        nexttile;
        imageHandle = imagesc(correlationMatrix, [-colorLimit colorLimit]);
        set(imageHandle, 'AlphaData', isfinite(correlationMatrix));
        set(gca, 'Color', [0.92 0.92 0.92]);
        colormap(local_diverging_colormap(256));
        colorbar;
        xticks(1:numel(metricNames));
        xticklabels(metricNames);
        xtickangle(35);
        yticks(1:numel(modelNames));
        yticklabels(modelNames);
        xlabel('Behaviour metric');
        ylabel('Model');
        title(sprintf('%s: %s vs G_{bal}', chamberName, predictorLabel), ...
            'Interpreter', 'tex');

        for modelCounter = 1:numel(modelNames)
            for metricCounter = 1:numel(metricNames)
                correlationValue = correlationMatrix(modelCounter, metricCounter);
                probabilityValue = probabilityMatrix(modelCounter, metricCounter);
                if ~isfinite(correlationValue)
                    continue;
                end

                significanceMarker = '';
                if isfinite(probabilityValue) && probabilityValue < 0.05
                    significanceMarker = '*';
                end

                textColor = [0 0 0];
                if abs(correlationValue) > 0.55 * colorLimit
                    textColor = [1 1 1];
                end

                text(metricCounter, modelCounter, sprintf('%.2f%s', correlationValue, significanceMarker), ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', ...
                    'FontSize', 7, ...
                    'Color', textColor);
            end
        end
    end
end

sgtitle('Behaviour correlations with Eq. 70 / NeuralJointGainBalance (G_{bal})', ...
    'Interpreter', 'tex');
figureStem = 'BehaviourNeural_Eq70_JointGainBalance_AllMetrics';
savefig(hFigure, fullfile(plotDir, [figureStem '.fig']));
exportgraphics(hFigure, fullfile(plotDir, [figureStem '.pdf']), 'ContentType', 'vector');
close(hFigure);
end


function local_plot_focus_model_correlations(jointTable, analysisPms)
if isempty(jointTable)
    return;
end

plotDir = fullfile(analysisPms.outputRoot, 'plots', 'focus_models');
if ~exist(plotDir, 'dir')
    mkdir(plotDir);
end

for pairCounter = 1:size(analysisPms.directionalPlotPairs, 1)
    predictorName = analysisPms.directionalPlotPairs{pairCounter, 1};
    targetName = analysisPms.directionalPlotPairs{pairCounter, 2};

    hFigure = figure('Color', 'w', 'Visible', 'off', ...
        'Name', sprintf('Focus %s vs %s', predictorName, targetName), ...
        'Position', [100 100 900 900]);
    tiledlayout(size(analysisPms.focusModelRows, 1), 1, ...
        'TileSpacing', 'compact', 'Padding', 'compact');

    for focusCounter = 1:size(analysisPms.focusModelRows, 1)
        chamberName = analysisPms.focusModelRows{focusCounter, 1};
        modelName = analysisPms.focusModelRows{focusCounter, 2};
        modelNote = analysisPms.focusModelRows{focusCounter, 3};
        groupMask = strcmp(jointTable.Chamber, chamberName) & strcmp(jointTable.Model, modelName);
        predictorValues = jointTable.(predictorName)(groupMask);
        targetValues = jointTable.(targetName)(groupMask);
        directionValues = jointTable.Direction(groupMask);
        validSamples = isfinite(predictorValues) & isfinite(targetValues);

        nexttile;
        hold on;
        scatter(predictorValues(validSamples), targetValues(validSamples), ...
            42, directionValues(validSamples), 'filled', ...
            'MarkerFaceAlpha', 0.72, 'MarkerEdgeColor', [0.2 0.2 0.2]);
        xline(0, ':', 'Color', [0.35 0.35 0.35]);
        yline(0, ':', 'Color', [0.35 0.35 0.35]);
        grid on;
        colormap(turbo(8));
        caxis([1 8]);

        rhoValue = NaN;
        probabilityValue = NaN;
        sampleCount = sum(validSamples);
        if sampleCount >= analysisPms.minimumSamplesForCorrelation
            [rhoValue, probabilityValue] = corr( ...
                predictorValues(validSamples), targetValues(validSamples), ...
                'Type', analysisPms.correlationType, ...
                'Rows', 'complete');
        end

        title(sprintf('%s %s: %s | rho = %.3f, p = %.2g, N = %d', ...
            chamberName, modelName, modelNote, rhoValue, probabilityValue, sampleCount), ...
            'Interpreter', 'none');
        xlabel(local_pretty_correlation_label(predictorName), 'Interpreter', 'tex');
        ylabel(local_pretty_correlation_label(targetName), 'Interpreter', 'tex');
        colorbar('Ticks', 1:8, 'TickLabels', compose('D%d', 1:8));
    end

    sgtitle(sprintf('Best-model focus: %s vs %s', ...
        local_pretty_correlation_label(predictorName), ...
        local_pretty_correlation_label(targetName)), ...
        'Interpreter', 'tex');

    figureStem = sprintf('BehaviourNeural_FocusModels_%s_vs_%s', predictorName, targetName);
    savefig(hFigure, fullfile(plotDir, [figureStem '.fig']));
    exportgraphics(hFigure, fullfile(plotDir, [figureStem '.pdf']), 'ContentType', 'vector');
    close(hFigure);
end
end


function local_plot_directional_correlations(directionalCorrelationTable, analysisPms)
if isempty(directionalCorrelationTable)
    return;
end

plotDir = fullfile(analysisPms.outputRoot, 'plots', 'directional_correlations');
if ~exist(plotDir, 'dir')
    mkdir(plotDir);
end

for pairCounter = 1:size(analysisPms.directionalPlotPairs, 1)
    predictorName = analysisPms.directionalPlotPairs{pairCounter, 1};
    targetName = analysisPms.directionalPlotPairs{pairCounter, 2};
    pairMask = strcmp(directionalCorrelationTable.Predictor, predictorName) & ...
        strcmp(directionalCorrelationTable.Target, targetName);
    pairTable = directionalCorrelationTable(pairMask, :);

    hFigure = figure('Color', 'w', 'Visible', 'off', ...
        'Name', sprintf('Directional %s vs %s', predictorName, targetName), ...
        'Position', [100 100 1180 520]);
    tiledlayout(1, numel(analysisPms.chamberNames), ...
        'TileSpacing', 'compact', 'Padding', 'compact');

    chamberMaxAbsRho = max(abs(pairTable.Rho), [], 'omitnan');
    if ~isfinite(chamberMaxAbsRho) || chamberMaxAbsRho == 0
        chamberMaxAbsRho = 0.1;
    end
    colorLimit = min(1, max(0.1, chamberMaxAbsRho));

    for chamberCounter = 1:numel(analysisPms.chamberNames)
        chamberName = analysisPms.chamberNames{chamberCounter};
        chamberMask = strcmp(pairTable.Chamber, chamberName);
        chamberTable = pairTable(chamberMask, :);
        modelNames = unique(chamberTable.Model, 'stable');
        directionIds = 1:8;
        correlationMatrix = nan(numel(modelNames), numel(directionIds));

        for modelCounter = 1:numel(modelNames)
            modelName = modelNames{modelCounter};
            for directionCounter = 1:numel(directionIds)
                directionId = directionIds(directionCounter);
                valueMask = strcmp(chamberTable.Model, modelName) & chamberTable.Direction == directionId;
                if any(valueMask)
                    correlationMatrix(modelCounter, directionCounter) = chamberTable.Rho(find(valueMask, 1, 'first'));
                end
            end
        end

        nexttile;
        imageHandle = imagesc(correlationMatrix, [-colorLimit colorLimit]);
        set(imageHandle, 'AlphaData', isfinite(correlationMatrix));
        set(gca, 'Color', [0.92 0.92 0.92]);
        colormap(local_diverging_colormap(256));
        colorbar;
        xticks(directionIds);
        xticklabels(compose('D%d', directionIds));
        yticks(1:numel(modelNames));
        yticklabels(modelNames);
        xlabel('Direction');
        ylabel('Model');
        title(chamberName, 'Interpreter', 'none');
    end

    sgtitle(sprintf('Directional Spearman rho: %s vs %s', ...
        local_pretty_correlation_label(predictorName), ...
        local_pretty_correlation_label(targetName)), ...
        'Interpreter', 'tex');

    figureStem = sprintf('BehaviourNeural_Directional_%s_vs_%s', predictorName, targetName);
    savefig(hFigure, fullfile(plotDir, [figureStem '.fig']));
    exportgraphics(hFigure, fullfile(plotDir, [figureStem '.pdf']), 'ContentType', 'vector');
    close(hFigure);
end
end


function local_plot_key_correlations(jointTable, correlationTable, analysisPms)
if isempty(jointTable)
    return;
end

plotDir = fullfile(analysisPms.outputRoot, 'plots');
if ~exist(plotDir, 'dir')
    mkdir(plotDir);
end

keyPredictors = {'SAdvantage_RT', 'SAdvantage_AE_abs', 'SJointBenefit_RT', 'KJointBenefit_RT'};
keyTargets = {'NeuralBalanceSToKMinusKToS', 'NeuralJointGainBalance'};

for predictorCounter = 1:numel(keyPredictors)
    predictorName = keyPredictors{predictorCounter};
    for targetCounter = 1:numel(keyTargets)
        targetName = keyTargets{targetCounter};
        hFigure = figure('Color', 'w', 'Visible', 'off', ...
            'Name', sprintf('%s vs %s', predictorName, targetName), ...
            'Position', [100 100 1100 520]);
        tiledlayout(1, numel(analysisPms.chamberNames), ...
            'TileSpacing', 'compact', 'Padding', 'compact');

        for chamberCounter = 1:numel(analysisPms.chamberNames)
            chamberName = analysisPms.chamberNames{chamberCounter};
            nexttile;
            hold on;
            modelNames = unique(jointTable.Model(strcmp(jointTable.Chamber, chamberName)), 'stable');
            colorOrder = lines(max(1, numel(modelNames)));

            for modelCounter = 1:numel(modelNames)
                modelName = modelNames{modelCounter};
                groupMask = strcmp(jointTable.Chamber, chamberName) & strcmp(jointTable.Model, modelName);
                scatter(jointTable.(predictorName)(groupMask), jointTable.(targetName)(groupMask), ...
                    28, colorOrder(modelCounter, :), 'filled', ...
                    'MarkerFaceAlpha', 0.65, ...
                    'DisplayName', modelName);
            end

            xlabel(local_pretty_correlation_label(predictorName), 'Interpreter', 'tex');
            ylabel(local_pretty_correlation_label(targetName), 'Interpreter', 'tex');
            title(chamberName, 'Interpreter', 'none');
            grid on;
            legend('Location', 'best');
        end

        titleText = sprintf('%s vs %s', ...
            local_pretty_correlation_label(predictorName), ...
            local_pretty_correlation_label(targetName));
        sgtitle(titleText, 'Interpreter', 'tex');
        figureStem = sprintf('BehaviourNeural_%s_vs_%s', predictorName, targetName);
        savefig(hFigure, fullfile(plotDir, [figureStem '.fig']));
        exportgraphics(hFigure, fullfile(plotDir, [figureStem '.pdf']), 'ContentType', 'vector');
        close(hFigure);
    end
end

if ~isempty(correlationTable)
    sortedCorrelationTable = sortrows(correlationTable, 'PValue', 'ascend');
    writetable(sortedCorrelationTable, fullfile(analysisPms.outputRoot, 'BehaviourNeural_correlations_sorted.csv'));
end
end


function local_plot_m7_rt_directional_stability( ...
    directionalAllTargetCorrelationTable, correlationTable, analysisPms)
if isempty(directionalAllTargetCorrelationTable) || isempty(correlationTable)
    return;
end

plotDir = fullfile(analysisPms.outputRoot, 'plots', 'stable_results');
if ~exist(plotDir, 'dir')
    mkdir(plotDir);
end

modelName = 'M7';
predictorName = 'SAdvantage_RT';
targetName = 'NeuralBalanceSToKMinusKToS';
figureStem = 'BehaviourNeural_M7_RT_NeuralBalance_DirectionalStability';
summaryRows = {};

hFigure = figure('Color', 'w', 'Visible', 'off', ...
    'Name', 'M7 RT directional stability', ...
    'Position', [100 100 1160 520]);
tiledlayout(1, numel(analysisPms.chamberNames), ...
    'TileSpacing', 'compact', 'Padding', 'compact');

for chamberCounter = 1:numel(analysisPms.chamberNames)
    chamberName = analysisPms.chamberNames{chamberCounter};
    directionMask = strcmp(directionalAllTargetCorrelationTable.Chamber, chamberName) & ...
        strcmp(directionalAllTargetCorrelationTable.Model, modelName) & ...
        strcmp(directionalAllTargetCorrelationTable.Predictor, predictorName) & ...
        strcmp(directionalAllTargetCorrelationTable.Target, targetName);
    chamberDirectionalTable = directionalAllTargetCorrelationTable(directionMask, :);
    chamberDirectionalTable = sortrows(chamberDirectionalTable, 'Direction');

    globalMask = strcmp(correlationTable.Chamber, chamberName) & ...
        strcmp(correlationTable.Model, modelName) & ...
        strcmp(correlationTable.Predictor, predictorName) & ...
        strcmp(correlationTable.Target, targetName);
    chamberGlobalTable = correlationTable(globalMask, :);

    if isempty(chamberDirectionalTable) || isempty(chamberGlobalTable)
        nexttile;
        title(sprintf('%s: missing data', chamberName), 'Interpreter', 'none');
        axis off;
        continue;
    end

    globalRho = chamberGlobalTable.Rho(1);
    globalPValue = chamberGlobalTable.PValue(1);
    globalSampleCount = chamberGlobalTable.N(1);
    directionIds = chamberDirectionalTable.Direction;
    directionRhos = chamberDirectionalTable.Rho;
    directionPValues = chamberDirectionalTable.PValue;
    directionSampleCounts = chamberDirectionalTable.N;

    for directionCounter = 1:numel(directionIds)
        summaryRows(end + 1, :) = { ...
            chamberName, modelName, predictorName, targetName, ...
            directionIds(directionCounter), directionRhos(directionCounter), ...
            directionPValues(directionCounter), directionSampleCounts(directionCounter), ...
            globalRho, globalPValue, globalSampleCount}; %#ok<AGROW>
    end

    nexttile;
    barHandle = bar(directionIds, directionRhos, ...
        'FaceColor', [0.19 0.36 0.70], ...
        'EdgeColor', 'none');
    hold on;
    yline(0, ':', 'Color', [0.35 0.35 0.35]);
    yline(globalRho, '--', 'Color', [0.65 0.18 0.18], 'LineWidth', 1.4);
    grid on;
    xticks(directionIds);
    xlabel('Direction');
    ylabel('Spearman rho');
    title(sprintf('%s M7: S RT advantage vs N_{bal}', chamberName), ...
        'Interpreter', 'tex');
    subtitle(sprintf('global rho = %.3f, p = %.2g, N = %d', ...
        globalRho, globalPValue, globalSampleCount), ...
        'Interpreter', 'none');

    for directionCounter = 1:numel(directionIds)
        rhoValue = directionRhos(directionCounter);
        if rhoValue < 0
            textVerticalAlignment = 'top';
            textYOffset = -0.025;
        else
            textVerticalAlignment = 'bottom';
            textYOffset = 0.025;
        end
        text(directionIds(directionCounter), rhoValue + textYOffset, ...
            sprintf('%.2f', rhoValue), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', textVerticalAlignment, ...
            'FontSize', 9);
    end

    minYLimit = min([directionRhos(:); globalRho; -0.05]);
    maxYLimit = max([directionRhos(:); globalRho; 0.05]);
    ylim([minYLimit - 0.08, maxYLimit + 0.08]);
    set(gca, 'FontSize', 10);
end

sgtitle('Direction-resolved stable RT effect in ModelM7', ...
    'Interpreter', 'none');
summaryTable = cell2table(summaryRows, ...
    'VariableNames', {'Chamber', 'Model', 'Predictor', 'Target', ...
    'Direction', 'DirectionalRho', 'DirectionalPValue', 'DirectionalN', ...
    'GlobalRho', 'GlobalPValue', 'GlobalN'});
writetable(summaryTable, fullfile(analysisPms.outputRoot, ...
    [figureStem '.csv']));
savefig(hFigure, fullfile(plotDir, [figureStem '.fig']));
exportgraphics(hFigure, fullfile(plotDir, [figureStem '.pdf']), ...
    'ContentType', 'vector');
close(hFigure);
end


function local_plot_m7_nonrt_directional_bars( ...
    directionalAllTargetCorrelationTable, correlationTable, analysisPms)
if isempty(directionalAllTargetCorrelationTable) || isempty(correlationTable)
    return;
end

plotDir = fullfile(analysisPms.outputRoot, 'plots', 'stable_results');
if ~exist(plotDir, 'dir')
    mkdir(plotDir);
end

modelName = 'M7';
targetName = 'NeuralJointGainBalance';
figureStem = 'BehaviourNeural_M7_NonRT_JointGainBalance_DirectionalBars';
barSpecs = { ...
    'SAdvantage_AE_abs', 'S AE advantage vs G_{bal}'; ...
    'SAdvantage_EC', 'S EC advantage vs G_{bal}'; ...
    'SAdvantage_ExT', 'S ExT advantage vs G_{bal}'};
summaryRows = {};

hFigure = figure('Color', 'w', 'Visible', 'off', ...
    'Name', 'M7 non-RT directional bars', ...
    'Position', [100 100 1180 980]);
tiledlayout(size(barSpecs, 1), numel(analysisPms.chamberNames), ...
    'TileSpacing', 'compact', 'Padding', 'compact');

for predictorCounter = 1:size(barSpecs, 1)
    predictorName = barSpecs{predictorCounter, 1};
    predictorLabel = barSpecs{predictorCounter, 2};

    for chamberCounter = 1:numel(analysisPms.chamberNames)
        chamberName = analysisPms.chamberNames{chamberCounter};
        directionMask = strcmp(directionalAllTargetCorrelationTable.Chamber, chamberName) & ...
            strcmp(directionalAllTargetCorrelationTable.Model, modelName) & ...
            strcmp(directionalAllTargetCorrelationTable.Predictor, predictorName) & ...
            strcmp(directionalAllTargetCorrelationTable.Target, targetName);
        chamberDirectionalTable = directionalAllTargetCorrelationTable(directionMask, :);
        chamberDirectionalTable = sortrows(chamberDirectionalTable, 'Direction');

        globalMask = strcmp(correlationTable.Chamber, chamberName) & ...
            strcmp(correlationTable.Model, modelName) & ...
            strcmp(correlationTable.Predictor, predictorName) & ...
            strcmp(correlationTable.Target, targetName);
        chamberGlobalTable = correlationTable(globalMask, :);

        nexttile;
        if isempty(chamberDirectionalTable) || isempty(chamberGlobalTable)
            title(sprintf('%s: missing %s', chamberName, predictorName), ...
                'Interpreter', 'none');
            axis off;
            continue;
        end

        globalRho = chamberGlobalTable.Rho(1);
        globalPValue = chamberGlobalTable.PValue(1);
        globalSampleCount = chamberGlobalTable.N(1);
        directionIds = chamberDirectionalTable.Direction;
        directionRhos = chamberDirectionalTable.Rho;
        directionPValues = chamberDirectionalTable.PValue;
        directionSampleCounts = chamberDirectionalTable.N;

        for directionCounter = 1:numel(directionIds)
            summaryRows(end + 1, :) = { ...
                chamberName, modelName, predictorName, targetName, ...
                directionIds(directionCounter), directionRhos(directionCounter), ...
                directionPValues(directionCounter), directionSampleCounts(directionCounter), ...
                globalRho, globalPValue, globalSampleCount}; %#ok<AGROW>
        end

        barHandle = bar(directionIds, directionRhos, ...
            'FaceColor', 'flat', ...
            'EdgeColor', 'none');
        for directionCounter = 1:numel(directionIds)
            if directionRhos(directionCounter) >= 0
                barHandle.CData(directionCounter, :) = [0.78 0.18 0.18];
            else
                barHandle.CData(directionCounter, :) = [0.16 0.34 0.75];
            end
        end
        hold on;
        yline(0, ':', 'Color', [0.35 0.35 0.35]);
        yline(globalRho, '--', 'Color', [0.20 0.20 0.20], 'LineWidth', 1.3);
        grid on;
        xticks(directionIds);
        xlabel('Direction');
        ylabel('Spearman rho');
        title(sprintf('%s M7: %s', chamberName, predictorLabel), ...
            'Interpreter', 'tex');
        subtitle(sprintf('global rho = %.3f, p = %.2g, N = %d', ...
            globalRho, globalPValue, globalSampleCount), ...
            'Interpreter', 'none');

        for directionCounter = 1:numel(directionIds)
            rhoValue = directionRhos(directionCounter);
            if rhoValue < 0
                textVerticalAlignment = 'top';
                textYOffset = -0.025;
            else
                textVerticalAlignment = 'bottom';
                textYOffset = 0.025;
            end
            text(directionIds(directionCounter), rhoValue + textYOffset, ...
                sprintf('%.2f', rhoValue), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', textVerticalAlignment, ...
                'FontSize', 8);
        end

        minYLimit = min([directionRhos(:); globalRho; -0.05]);
        maxYLimit = max([directionRhos(:); globalRho; 0.05]);
        ylim([minYLimit - 0.08, maxYLimit + 0.08]);
        set(gca, 'FontSize', 9);
    end
end

sgtitle('Direction-resolved non-RT behaviour effects in ModelM7', ...
    'Interpreter', 'none');
summaryTable = cell2table(summaryRows, ...
    'VariableNames', {'Chamber', 'Model', 'Predictor', 'Target', ...
    'Direction', 'DirectionalRho', 'DirectionalPValue', 'DirectionalN', ...
    'GlobalRho', 'GlobalPValue', 'GlobalN'});
writetable(summaryTable, fullfile(analysisPms.outputRoot, ...
    [figureStem '.csv']));
savefig(hFigure, fullfile(plotDir, [figureStem '.fig']));
exportgraphics(hFigure, fullfile(plotDir, [figureStem '.pdf']), ...
    'ContentType', 'vector');
close(hFigure);
end


function local_save_directional_mean_abs_report( ...
    directionalAllTargetCorrelationTable, correlationTable, analysisPms)
if isempty(directionalAllTargetCorrelationTable)
    return;
end

reportDir = fullfile(analysisPms.outputRoot, 'DirectionalMeanAbs');
if ~exist(reportDir, 'dir')
    mkdir(reportDir);
end

primaryDirectionalTable = local_filter_table_models( ...
    directionalAllTargetCorrelationTable, analysisPms.primaryModelNames);
chamberNames = unique(primaryDirectionalTable.Chamber, 'stable');
modelNames = unique(primaryDirectionalTable.Model, 'stable');
predictorNames = unique(primaryDirectionalTable.Predictor, 'stable');
targetNames = unique(primaryDirectionalTable.Target, 'stable');
summaryRows = {};

for chamberCounter = 1:numel(chamberNames)
    chamberName = chamberNames{chamberCounter};
    for modelCounter = 1:numel(modelNames)
        modelName = modelNames{modelCounter};
        for predictorCounter = 1:numel(predictorNames)
            predictorName = predictorNames{predictorCounter};
            for targetCounter = 1:numel(targetNames)
                targetName = targetNames{targetCounter};
                directionMask = strcmp(primaryDirectionalTable.Chamber, chamberName) & ...
                    strcmp(primaryDirectionalTable.Model, modelName) & ...
                    strcmp(primaryDirectionalTable.Predictor, predictorName) & ...
                    strcmp(primaryDirectionalTable.Target, targetName);
                directionTable = primaryDirectionalTable(directionMask, :);
                if isempty(directionTable)
                    continue;
                end

                directionTable = sortrows(directionTable, 'Direction');
                finiteMask = isfinite(directionTable.Rho);
                directionTable = directionTable(finiteMask, :);
                directionCount = height(directionTable);
                if directionCount < 8
                    continue;
                end

                directionRhos = directionTable.Rho;
                directionAbsRhos = abs(directionRhos);
                positiveCount = sum(directionRhos > 0);
                negativeCount = sum(directionRhos < 0);
                majoritySignFraction = max(positiveCount, negativeCount) / directionCount;
                sameSignAllDirections = positiveCount == directionCount || ...
                    negativeCount == directionCount;
                globalMask = strcmp(correlationTable.Chamber, chamberName) & ...
                    strcmp(correlationTable.Model, modelName) & ...
                    strcmp(correlationTable.Predictor, predictorName) & ...
                    strcmp(correlationTable.Target, targetName);
                globalTable = correlationTable(globalMask, :);
                globalRho = NaN;
                globalPValue = NaN;
                globalSampleCount = NaN;
                if ~isempty(globalTable)
                    globalRho = globalTable.Rho(1);
                    globalPValue = globalTable.PValue(1);
                    globalSampleCount = globalTable.N(1);
                end

                directionSummary = strings(directionCount, 1);
                for directionCounter = 1:directionCount
                    directionSummary(directionCounter) = sprintf('D%d=%+.2f', ...
                        directionTable.Direction(directionCounter), ...
                        directionRhos(directionCounter));
                end

                summaryRows(end + 1, :) = { ...
                    chamberName, modelName, predictorName, targetName, ...
                    directionCount, mean(directionAbsRhos, 'omitnan'), ...
                    mean(directionRhos, 'omitnan'), ...
                    median(directionAbsRhos, 'omitnan'), ...
                    min(directionAbsRhos), max(directionAbsRhos), ...
                    majoritySignFraction, sameSignAllDirections, ...
                    positiveCount, negativeCount, ...
                    globalRho, globalPValue, globalSampleCount, ...
                    strjoin(directionSummary, ' ')}; %#ok<AGROW>
            end
        end
    end
end

if isempty(summaryRows)
    return;
end

summaryTable = cell2table(summaryRows, ...
    'VariableNames', {'Chamber', 'Model', 'Predictor', 'Target', ...
    'DirectionCount', 'MeanAbsRho', 'SignedMeanRho', 'MedianAbsRho', ...
    'MinAbsRho', 'MaxAbsRho', 'MajoritySignFraction', ...
    'SameSignAllDirections', 'PositiveDirectionCount', ...
    'NegativeDirectionCount', 'GlobalRho', 'GlobalPValue', 'GlobalN', ...
    'DirectionRhos'});
summaryTable.NegativeMeanAbsRho = -summaryTable.MeanAbsRho;
summaryTable.NegativeMedianAbsRho = -summaryTable.MedianAbsRho;
summaryTable = sortrows(summaryTable, ...
    {'Chamber', 'NegativeMeanAbsRho', 'NegativeMedianAbsRho'});
summaryTable.NegativeMeanAbsRho = [];
summaryTable.NegativeMedianAbsRho = [];

writetable(summaryTable, fullfile(reportDir, ...
    'BehaviourNeural_directionalMeanAbsCorrelations.csv'));

topCount = min(60, height(summaryTable));
topTable = summaryTable(1:topCount, :);
writetable(topTable, fullfile(reportDir, ...
    'BehaviourNeural_topDirectionalMeanAbsCorrelations.csv'));
end


function local_plot_best_receiver_directional_bars( ...
    directionalAllTargetCorrelationTable, correlationTable, analysisPms)
if isempty(directionalAllTargetCorrelationTable) || isempty(correlationTable)
    return;
end

plotDir = fullfile(analysisPms.outputRoot, 'plots', 'stable_results');
if ~exist(plotDir, 'dir')
    mkdir(plotDir);
end

figureStem = 'BehaviourNeural_BestReceiverBenefit_DirectionalBars';
bestSpecs = { ...
    'Frontal',  'S', 'M7', 'SJointBenefit_PV',  'JointWeightKToS', ...
        'B_S^{(PV)} vs I_S'; ...
    'Frontal',  'K', 'M2', 'KJointBenefit_ExT', 'JointGainSToK', ...
        'B_K^{(ExT)} vs \Delta I_K'; ...
    'Parietal', 'S', 'M7', 'SJointBenefit_RT',  'JointWeightKToS', ...
        'B_S^{(RT)} vs I_S'; ...
    'Parietal', 'K', 'M2', 'KJointBenefit_CMT', 'JointWeightSToK', ...
        'B_K^{(CMT)} vs I_K'};
summaryRows = {};

hFigure = figure('Color', 'w', 'Visible', 'off', ...
    'Name', 'Best receiver benefit directional bars', ...
    'Position', [100 100 1180 760]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for monkeyCounter = 1:2
    if monkeyCounter == 1
        monkeyName = 'S';
    else
        monkeyName = 'K';
    end

    for chamberCounter = 1:numel(analysisPms.chamberNames)
        chamberName = analysisPms.chamberNames{chamberCounter};
        specMask = strcmp(bestSpecs(:, 1), chamberName) & ...
            strcmp(bestSpecs(:, 2), monkeyName);
        specRow = bestSpecs(specMask, :);
        nexttile;

        if isempty(specRow)
            title(sprintf('%s %s: missing spec', chamberName, monkeyName), ...
                'Interpreter', 'none');
            axis off;
            continue;
        end

        modelName = specRow{1, 3};
        predictorName = specRow{1, 4};
        targetName = specRow{1, 5};
        plotLabel = specRow{1, 6};

        directionMask = strcmp(directionalAllTargetCorrelationTable.Chamber, chamberName) & ...
            strcmp(directionalAllTargetCorrelationTable.Model, modelName) & ...
            strcmp(directionalAllTargetCorrelationTable.Predictor, predictorName) & ...
            strcmp(directionalAllTargetCorrelationTable.Target, targetName);
        directionTable = directionalAllTargetCorrelationTable(directionMask, :);
        directionTable = sortrows(directionTable, 'Direction');

        globalMask = strcmp(correlationTable.Chamber, chamberName) & ...
            strcmp(correlationTable.Model, modelName) & ...
            strcmp(correlationTable.Predictor, predictorName) & ...
            strcmp(correlationTable.Target, targetName);
        globalTable = correlationTable(globalMask, :);

        if isempty(directionTable) || isempty(globalTable)
            title(sprintf('%s %s: missing data', chamberName, monkeyName), ...
                'Interpreter', 'none');
            axis off;
            continue;
        end

        directionIds = directionTable.Direction;
        directionRhos = directionTable.Rho;
        directionPValues = directionTable.PValue;
        directionSampleCounts = directionTable.N;
        globalRho = globalTable.Rho(1);
        globalPValue = globalTable.PValue(1);
        globalSampleCount = globalTable.N(1);
        meanAbsRho = mean(abs(directionRhos), 'omitnan');
        signedMeanRho = mean(directionRhos, 'omitnan');
        positiveCount = sum(directionRhos > 0);
        negativeCount = sum(directionRhos < 0);
        majoritySignFraction = max(positiveCount, negativeCount) / numel(directionRhos);

        for directionCounter = 1:numel(directionIds)
            summaryRows(end + 1, :) = { ...
                chamberName, monkeyName, modelName, predictorName, targetName, ...
                directionIds(directionCounter), directionRhos(directionCounter), ...
                directionPValues(directionCounter), directionSampleCounts(directionCounter), ...
                meanAbsRho, signedMeanRho, majoritySignFraction, ...
                globalRho, globalPValue, globalSampleCount}; %#ok<AGROW>
        end

        barHandle = bar(directionIds, directionRhos, ...
            'FaceColor', 'flat', ...
            'EdgeColor', 'none');
        for directionCounter = 1:numel(directionIds)
            if directionRhos(directionCounter) >= 0
                barHandle.CData(directionCounter, :) = [0.78 0.18 0.18];
            else
                barHandle.CData(directionCounter, :) = [0.16 0.34 0.75];
            end
        end
        hold on;
        yline(0, ':', 'Color', [0.35 0.35 0.35]);
        yline(globalRho, '--', 'Color', [0.20 0.20 0.20], 'LineWidth', 1.3);
        grid on;
        xticks(directionIds);
        xlabel('Direction');
        ylabel('Spearman rho');
        title(sprintf('%s %s (%s): %s', chamberName, monkeyName, modelName, plotLabel), ...
            'Interpreter', 'tex');
        subtitle(sprintf('mean |rho|=%.3f, signed mean=%.3f, global rho=%.3f, p=%.2g', ...
            meanAbsRho, signedMeanRho, globalRho, globalPValue), ...
            'Interpreter', 'none');

        for directionCounter = 1:numel(directionIds)
            rhoValue = directionRhos(directionCounter);
            if rhoValue < 0
                textVerticalAlignment = 'top';
                textYOffset = -0.025;
            else
                textVerticalAlignment = 'bottom';
                textYOffset = 0.025;
            end
            text(directionIds(directionCounter), rhoValue + textYOffset, ...
                sprintf('%.2f', rhoValue), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', textVerticalAlignment, ...
                'FontSize', 8);
        end

        minYLimit = min([directionRhos(:); globalRho; -0.05]);
        maxYLimit = max([directionRhos(:); globalRho; 0.05]);
        ylim([minYLimit - 0.08, maxYLimit + 0.08]);
        set(gca, 'FontSize', 8);
    end
end

sgtitle('Best receiver-specific behaviour-neural directional correlations', ...
    'Interpreter', 'none');
summaryTable = cell2table(summaryRows, ...
    'VariableNames', {'Chamber', 'Monkey', 'Model', 'Predictor', 'Target', ...
    'Direction', 'DirectionalRho', 'DirectionalPValue', 'DirectionalN', ...
    'MeanAbsRho', 'SignedMeanRho', 'MajoritySignFraction', ...
    'GlobalRho', 'GlobalPValue', 'GlobalN'});
writetable(summaryTable, fullfile(analysisPms.outputRoot, ...
    [figureStem '.csv']));
savefig(hFigure, fullfile(plotDir, [figureStem '.fig']));
exportgraphics(hFigure, fullfile(plotDir, [figureStem '.pdf']), ...
    'ContentType', 'vector');
close(hFigure);
end


function local_plot_k_receiver_directional_bars( ...
    directionalAllTargetCorrelationTable, correlationTable, analysisPms)
if isempty(directionalAllTargetCorrelationTable) || isempty(correlationTable)
    return;
end

plotDir = fullfile(analysisPms.outputRoot, 'plots', 'stable_results');
if ~exist(plotDir, 'dir')
    mkdir(plotDir);
end

figureStem = 'BehaviourNeural_BestKReceiverBenefit_DirectionalBars';
barSpecs = { ...
    'Frontal',  'M2', 'KJointBenefit_CMT', 'JointWeightSToK', ...
        'Frontal K: B_K^{(CMT)} vs I_K'; ...
    'Frontal',  'M2', 'KJointBenefit_EC',  'JointGainSToK', ...
        'Frontal K: B_K^{(EC)} vs \Delta I_K'; ...
    'Parietal', 'M2', 'KJointBenefit_CMT', 'JointWeightSToK', ...
        'Parietal K: B_K^{(CMT)} vs I_K'; ...
    'Parietal', 'M2', 'KJointBenefit_EC',  'JointGainSToK', ...
        'Parietal K: B_K^{(EC)} vs \Delta I_K'};
summaryRows = {};

hFigure = figure('Color', 'w', 'Visible', 'off', ...
    'Name', 'Best K receiver benefit directional bars', ...
    'Position', [100 100 1180 760]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for specCounter = 1:size(barSpecs, 1)
    chamberName = barSpecs{specCounter, 1};
    modelName = barSpecs{specCounter, 2};
    predictorName = barSpecs{specCounter, 3};
    targetName = barSpecs{specCounter, 4};
    plotLabel = barSpecs{specCounter, 5};

    directionMask = strcmp(directionalAllTargetCorrelationTable.Chamber, chamberName) & ...
        strcmp(directionalAllTargetCorrelationTable.Model, modelName) & ...
        strcmp(directionalAllTargetCorrelationTable.Predictor, predictorName) & ...
        strcmp(directionalAllTargetCorrelationTable.Target, targetName);
    directionTable = directionalAllTargetCorrelationTable(directionMask, :);
    directionTable = sortrows(directionTable, 'Direction');

    globalMask = strcmp(correlationTable.Chamber, chamberName) & ...
        strcmp(correlationTable.Model, modelName) & ...
        strcmp(correlationTable.Predictor, predictorName) & ...
        strcmp(correlationTable.Target, targetName);
    globalTable = correlationTable(globalMask, :);

    nexttile;
    if isempty(directionTable) || isempty(globalTable)
        title(sprintf('%s: missing %s', chamberName, predictorName), ...
            'Interpreter', 'none');
        axis off;
        continue;
    end

    directionIds = directionTable.Direction;
    directionRhos = directionTable.Rho;
    directionPValues = directionTable.PValue;
    directionSampleCounts = directionTable.N;
    globalRho = globalTable.Rho(1);
    globalPValue = globalTable.PValue(1);
    globalSampleCount = globalTable.N(1);
    meanAbsRho = mean(abs(directionRhos), 'omitnan');
    signedMeanRho = mean(directionRhos, 'omitnan');
    positiveCount = sum(directionRhos > 0);
    negativeCount = sum(directionRhos < 0);
    majoritySignFraction = max(positiveCount, negativeCount) / numel(directionRhos);

    for directionCounter = 1:numel(directionIds)
        summaryRows(end + 1, :) = { ...
            chamberName, 'K', modelName, predictorName, targetName, ...
            directionIds(directionCounter), directionRhos(directionCounter), ...
            directionPValues(directionCounter), directionSampleCounts(directionCounter), ...
            meanAbsRho, signedMeanRho, majoritySignFraction, ...
            globalRho, globalPValue, globalSampleCount}; %#ok<AGROW>
    end

    barHandle = bar(directionIds, directionRhos, ...
        'FaceColor', 'flat', ...
        'EdgeColor', 'none');
    for directionCounter = 1:numel(directionIds)
        if directionRhos(directionCounter) >= 0
            barHandle.CData(directionCounter, :) = [0.78 0.18 0.18];
        else
            barHandle.CData(directionCounter, :) = [0.16 0.34 0.75];
        end
    end
    hold on;
    yline(0, ':', 'Color', [0.35 0.35 0.35]);
    yline(globalRho, '--', 'Color', [0.20 0.20 0.20], 'LineWidth', 1.3);
    grid on;
    xticks(directionIds);
    xlabel('Direction');
    ylabel('Spearman rho');
    title(sprintf('%s (%s): %s', chamberName, modelName, plotLabel), ...
        'Interpreter', 'tex');
    subtitle(sprintf('mean |rho|=%.3f, signed mean=%.3f, global rho=%.3f, p=%.2g', ...
        meanAbsRho, signedMeanRho, globalRho, globalPValue), ...
        'Interpreter', 'none');

    for directionCounter = 1:numel(directionIds)
        rhoValue = directionRhos(directionCounter);
        if rhoValue < 0
            textVerticalAlignment = 'top';
            textYOffset = -0.025;
        else
            textVerticalAlignment = 'bottom';
            textYOffset = 0.025;
        end
        text(directionIds(directionCounter), rhoValue + textYOffset, ...
            sprintf('%.2f', rhoValue), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', textVerticalAlignment, ...
            'FontSize', 8);
    end

    minYLimit = min([directionRhos(:); globalRho; -0.05]);
    maxYLimit = max([directionRhos(:); globalRho; 0.05]);
    ylim([minYLimit - 0.08, maxYLimit + 0.08]);
    set(gca, 'FontSize', 8);
end

sgtitle('K receiver-specific behaviour-neural directional correlations', ...
    'Interpreter', 'none');
summaryTable = cell2table(summaryRows, ...
    'VariableNames', {'Chamber', 'Monkey', 'Model', 'Predictor', 'Target', ...
    'Direction', 'DirectionalRho', 'DirectionalPValue', 'DirectionalN', ...
    'MeanAbsRho', 'SignedMeanRho', 'MajoritySignFraction', ...
    'GlobalRho', 'GlobalPValue', 'GlobalN'});
writetable(summaryTable, fullfile(analysisPms.outputRoot, ...
    [figureStem '.csv']));
savefig(hFigure, fullfile(plotDir, [figureStem '.fig']));
exportgraphics(hFigure, fullfile(plotDir, [figureStem '.pdf']), ...
    'ContentType', 'vector');
close(hFigure);
end


function local_plot_receiver_benefit_correlations(jointTable, analysisPms)
if isempty(jointTable)
    return;
end

plotDir = fullfile(analysisPms.outputRoot, 'plots', 'receiver_benefit');
if ~exist(plotDir, 'dir')
    mkdir(plotDir);
end

figureSpecs = { ...
    'JointWeight', 'joint incoming weight', { ...
        'KJointBenefit_RT', 'JointWeightSToK', 'B_K^{(RT)} vs I_K'; ...
        'SJointBenefit_RT', 'JointWeightKToS', 'B_S^{(RT)} vs I_S'}; ...
    'JointGain', 'joint-minus-source-active incoming gain', { ...
        'KJointBenefit_RT', 'JointGainSToK', 'B_K^{(RT)} vs \Delta I_K'; ...
        'SJointBenefit_RT', 'JointGainKToS', 'B_S^{(RT)} vs \Delta I_S'}};

for figureCounter = 1:size(figureSpecs, 1)
    figureStemSuffix = figureSpecs{figureCounter, 1};
    figureTitle = figureSpecs{figureCounter, 2};
    receiverPairs = figureSpecs{figureCounter, 3};

    hFigure = figure('Color', 'w', 'Visible', 'off', ...
        'Name', sprintf('Receiver benefit %s', figureTitle), ...
        'Position', [100 100 1180 780]);
    tiledlayout(size(receiverPairs, 1), numel(analysisPms.chamberNames), ...
        'TileSpacing', 'compact', 'Padding', 'compact');

    for pairCounter = 1:size(receiverPairs, 1)
        predictorName = receiverPairs{pairCounter, 1};
        targetName = receiverPairs{pairCounter, 2};
        pairTitle = receiverPairs{pairCounter, 3};

        for chamberCounter = 1:numel(analysisPms.chamberNames)
            chamberName = analysisPms.chamberNames{chamberCounter};
            nexttile;
            hold on;

            chamberMask = strcmp(jointTable.Chamber, chamberName);
            modelNames = unique(jointTable.Model(chamberMask), 'stable');
            colorOrder = lines(max(1, numel(modelNames)));

            for modelCounter = 1:numel(modelNames)
                modelName = modelNames{modelCounter};
                modelMask = chamberMask & strcmp(jointTable.Model, modelName);
                scatter(jointTable.(predictorName)(modelMask), jointTable.(targetName)(modelMask), ...
                    28, colorOrder(modelCounter, :), 'filled', ...
                    'MarkerFaceAlpha', 0.60, ...
                    'DisplayName', modelName);
            end

            predictorValues = jointTable.(predictorName)(chamberMask);
            targetValues = jointTable.(targetName)(chamberMask);
            validSamples = isfinite(predictorValues) & isfinite(targetValues);
            sampleCount = sum(validSamples);
            rhoValue = NaN;
            probabilityValue = NaN;
            if sampleCount >= analysisPms.minimumSamplesForCorrelation
                [rhoValue, probabilityValue] = corr( ...
                    predictorValues(validSamples), ...
                    targetValues(validSamples), ...
                    'Type', analysisPms.correlationType, ...
                    'Rows', 'complete');
            end

            xline(0, ':', 'Color', [0.35 0.35 0.35]);
            yline(0, ':', 'Color', [0.35 0.35 0.35]);
            grid on;
            xlabel(local_pretty_correlation_label(predictorName), 'Interpreter', 'tex');
            ylabel(local_pretty_correlation_label(targetName), 'Interpreter', 'tex');
            title(sprintf('%s: %s | rho = %.3f, p = %.2g, N = %d', ...
                chamberName, pairTitle, rhoValue, probabilityValue, sampleCount), ...
                'Interpreter', 'tex');
            if pairCounter == 1 && chamberCounter == numel(analysisPms.chamberNames)
                legend('Location', 'bestoutside');
            end
        end
    end

    sgtitle(sprintf('Receiver benefit vs neural %s', figureTitle), ...
        'Interpreter', 'none');
    figureStem = sprintf('BehaviourNeural_ReceiverBenefit_RT_%s', figureStemSuffix);
    savefig(hFigure, fullfile(plotDir, [figureStem '.fig']));
    exportgraphics(hFigure, fullfile(plotDir, [figureStem '.pdf']), 'ContentType', 'vector');
    close(hFigure);
end
end


function local_save_top_correlation_report(correlationTable, directionalCorrelationTable, ...
    jointGainBalanceCorrelationTable, analysisPms)

reportDir = fullfile(analysisPms.outputRoot, 'TopCorrelations');
plotDir = fullfile(reportDir, 'plots');
if ~exist(plotDir, 'dir')
    mkdir(plotDir);
end

topCount = 20;
primaryCorrelationTable = local_filter_table_models(correlationTable, analysisPms.primaryModelNames);
primaryJointGainBalanceCorrelationTable = local_filter_table_models( ...
    jointGainBalanceCorrelationTable, analysisPms.primaryModelNames);
primaryDirectionalCorrelationTable = local_filter_table_models( ...
    directionalCorrelationTable, analysisPms.primaryModelNames);

topEq70Table = local_select_top_correlation_rows(primaryJointGainBalanceCorrelationTable, topCount);
topGlobalTable = local_select_top_correlation_rows(primaryCorrelationTable, topCount);
topDirectionalTable = local_select_top_correlation_rows(primaryDirectionalCorrelationTable, topCount);
topEq70AllModelsTable = local_select_top_correlation_rows(jointGainBalanceCorrelationTable, topCount);
topGlobalAllModelsTable = local_select_top_correlation_rows(correlationTable, topCount);
topDirectionalAllModelsTable = local_select_top_correlation_rows(directionalCorrelationTable, topCount);

writetable(topEq70Table, fullfile(reportDir, 'BehaviourNeural_topEq70Correlations.csv'));
writetable(topGlobalTable, fullfile(reportDir, 'BehaviourNeural_topGlobalCorrelations.csv'));
writetable(topDirectionalTable, fullfile(reportDir, 'BehaviourNeural_topDirectionalCorrelations.csv'));
writetable(topEq70AllModelsTable, fullfile(reportDir, 'BehaviourNeural_topEq70AllModelsCorrelations.csv'));
writetable(topGlobalAllModelsTable, fullfile(reportDir, 'BehaviourNeural_topGlobalAllModelsCorrelations.csv'));
writetable(topDirectionalAllModelsTable, fullfile(reportDir, 'BehaviourNeural_topDirectionalAllModelsCorrelations.csv'));

local_plot_top_correlation_bars(topEq70Table, ...
    fullfile(plotDir, 'BehaviourNeural_TopEq70Correlations'), ...
    'Top primary-model global correlations with Eq. 70 / G_{bal}');
local_plot_top_correlation_bars(topGlobalTable, ...
    fullfile(plotDir, 'BehaviourNeural_TopGlobalCorrelations'), ...
    'Top primary-model global correlations across neural readouts');
local_plot_top_correlation_bars(topDirectionalTable, ...
    fullfile(plotDir, 'BehaviourNeural_TopDirectionalCorrelations'), ...
    'Top primary-model direction-specific correlations across neural readouts');
end


function filteredTable = local_filter_table_models(inputTable, modelNames)
if isempty(inputTable)
    filteredTable = inputTable;
    return;
end

modelMask = ismember(inputTable.Model, modelNames);
filteredTable = inputTable(modelMask, :);
end


function topTable = local_select_top_correlation_rows(inputTable, topCount)
if isempty(inputTable)
    topTable = inputTable;
    return;
end

validRows = isfinite(inputTable.Rho) & isfinite(inputTable.PValue);
topTable = inputTable(validRows, :);
if isempty(topTable)
    return;
end

topTable.AbsRho = abs(topTable.Rho);
topTable.NegAbsRho = -topTable.AbsRho;
topTable = sortrows(topTable, {'PValue', 'NegAbsRho'});
topTable.NegAbsRho = [];

if height(topTable) > topCount
    topTable = topTable(1:topCount, :);
end
end


function local_plot_top_correlation_bars(topTable, figureStem, figureTitle)
if isempty(topTable)
    return;
end

rowCount = height(topTable);
rhoValues = flipud(topTable.Rho(:));
probabilityValues = flipud(topTable.PValue(:));
sampleCounts = flipud(topTable.N(:));
absRhoValues = abs(rhoValues);
labelStrings = strings(rowCount, 1);

for rowCounter = 1:rowCount
    tableRow = topTable(rowCount - rowCounter + 1, :);
    labelStrings(rowCounter) = local_top_correlation_label(tableRow);
end

hFigure = figure('Color', 'w', 'Visible', 'off', ...
    'Name', figureTitle, ...
    'Position', [100 100 1380 820]);
hAxes = axes(hFigure);
hold(hAxes, 'on');

barHandle = barh(hAxes, 1:rowCount, rhoValues, ...
    'FaceColor', 'flat', ...
    'EdgeColor', 'none');
for rowCounter = 1:rowCount
    if rhoValues(rowCounter) >= 0
        barHandle.CData(rowCounter, :) = [0.78 0.18 0.18];
    else
        barHandle.CData(rowCounter, :) = [0.16 0.34 0.75];
    end
end

xline(hAxes, 0, ':', 'Color', [0.25 0.25 0.25]);
grid(hAxes, 'on');
set(hAxes, ...
    'YTick', 1:rowCount, ...
    'YTickLabel', labelStrings, ...
    'TickLabelInterpreter', 'tex', ...
    'YLim', [0 rowCount + 1]);
xlabel(hAxes, 'Spearman rho');
title(hAxes, figureTitle, 'Interpreter', 'tex');

for rowCounter = 1:rowCount
    labelOffset = 0.015 * max(1, max(absRhoValues));
    if rhoValues(rowCounter) >= 0
        textX = rhoValues(rowCounter) + labelOffset;
        horizontalAlignment = 'left';
    else
        textX = rhoValues(rowCounter) - labelOffset;
        horizontalAlignment = 'right';
    end

    text(hAxes, textX, rowCounter, ...
        sprintf('p=%.2g, N=%d', probabilityValues(rowCounter), sampleCounts(rowCounter)), ...
        'HorizontalAlignment', horizontalAlignment, ...
        'VerticalAlignment', 'middle', ...
        'FontSize', 8);
end

maxAbsRho = max(absRhoValues);
if ~isfinite(maxAbsRho) || maxAbsRho == 0
    maxAbsRho = 0.1;
end
xlim(hAxes, [-maxAbsRho maxAbsRho] * 1.22);

savefig(hFigure, [figureStem '.fig']);
exportgraphics(hFigure, [figureStem '.pdf'], 'ContentType', 'vector');
close(hFigure);
end


function labelString = local_top_correlation_label(tableRow)
chamberName = char(string(tableRow.Chamber));
modelName = char(string(tableRow.Model));
predictorName = char(string(tableRow.Predictor));
targetName = char(string(tableRow.Target));

if any(strcmp(tableRow.Properties.VariableNames, 'Direction'))
    directionText = sprintf(' D%d', tableRow.Direction);
else
    directionText = '';
end

labelString = sprintf('%s %s%s: %s \\rightarrow %s', ...
    chamberName, modelName, directionText, ...
    local_pretty_correlation_label(predictorName), ...
    local_pretty_correlation_label(targetName));
end


function labelString = local_pretty_correlation_label(rawName)
rawName = char(string(rawName));

switch rawName
    case 'NeuralBalanceSToKMinusKToS'
        labelString = 'N_{bal}';
        return;
    case 'NeuralJointGainBalance'
        labelString = 'G_{bal}';
        return;
    case 'JointWeightSToK'
        labelString = 'I_K';
        return;
    case 'JointWeightKToS'
        labelString = 'I_S';
        return;
    case 'JointGainSToK'
        labelString = '\Delta I_K';
        return;
    case 'JointGainKToS'
        labelString = '\Delta I_S';
        return;
end

tokenParts = regexp(rawName, '^(SAdvantage|SJointBenefit|KJointBenefit)_(.+)$', ...
    'tokens', 'once');
if ~isempty(tokenParts)
    metricName = tokenParts{2};
    switch tokenParts{1}
        case 'SAdvantage'
            labelString = sprintf('A_S^{(%s)}', metricName);
        case 'SJointBenefit'
            labelString = sprintf('B_S^{(%s)}', metricName);
        case 'KJointBenefit'
            labelString = sprintf('B_K^{(%s)}', metricName);
    end
    return;
end

labelString = strrep(rawName, '_', '\_');
end


function colorMap = local_diverging_colormap(colorCount)
halfColorCount = colorCount / 2;
blueToWhite = [linspace(0.16, 1, halfColorCount)', ...
    linspace(0.34, 1, halfColorCount)', ...
    linspace(0.75, 1, halfColorCount)'];
whiteToRed = [linspace(1, 0.78, halfColorCount)', ...
    linspace(1, 0.18, halfColorCount)', ...
    linspace(1, 0.18, halfColorCount)'];
colorMap = [blueToWhite; whiteToRed];
end
