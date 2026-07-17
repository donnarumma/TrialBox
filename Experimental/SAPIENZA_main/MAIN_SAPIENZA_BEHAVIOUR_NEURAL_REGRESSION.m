%% MAIN_SAPIENZA_BEHAVIOUR_NEURAL_REGRESSION
% Leave-one-session-out regressions for continuous behaviour-neural readouts.
%
% This complements the median-split classification checks.  Here the target
% remains continuous, so the key metric is predictive R^2 on held-out sessions.

clear;
close all;

analysisPms = struct();
analysisPms.inputRoot = '/TESTS/SAPIENZA/BEHAVIOUR/NeuralCorrelation';
analysisPms.outputRoot = '/TESTS/SAPIENZA/BEHAVIOUR/NeuralRegression';
analysisPms.visible = 'off';
analysisPms.minimumSamples = 16;
analysisPms.minimumSessions = 4;
analysisPms.ridgeLambdas = [0.01 0.1 1 10 100];
analysisPms.defaultRidgeLambda = 1;
analysisPms.modelRows = cell2table({ ...
    'Frontal',  'M2', 'Frontal M2 direct-coupling comparator'; ...
    'Frontal',  'M4', 'Frontal M4 trajectory comparator'; ...
    'Frontal',  'M7', 'Frontal M7 primary trajectory model'; ...
    'Frontal',  'M8', 'Frontal M8 trajectory/null competitor'; ...
    'Parietal', 'M2', 'Parietal M2 primary direct-coupling model'; ...
    'Parietal', 'M4', 'Parietal M4 trajectory comparator'; ...
    'Parietal', 'M7', 'Parietal M7 trajectory comparator'; ...
    'Parietal', 'M8', 'Parietal M8 trajectory/null competitor'}, ...
    'VariableNames', {'Chamber', 'Model', 'Description'});
analysisPms.targetRows = cell2table({ ...
    'S_RT_Advantage',       'SAdvantage_RT',        'S RT advantage in joint action'; ...
    'S_RT_Benefit',         'SJointBenefit_RT',     'S RT joint benefit'; ...
    'K_RT_Benefit',         'KJointBenefit_RT',     'K RT joint benefit'; ...
    'S_AE_Benefit',         'SJointBenefit_AE_abs', 'S absolute-angular-error joint benefit'; ...
    'K_AE_Benefit',         'KJointBenefit_AE_abs', 'K absolute-angular-error joint benefit'; ...
    'S_AE_Advantage',       'SAdvantage_AE_abs',    'S absolute-angular-error advantage in joint action'}, ...
    'VariableNames', {'TargetId', 'TargetVariable', 'Description'});
analysisPms.featureRows = local_feature_rows();

if ~exist(analysisPms.outputRoot, 'dir')
    mkdir(analysisPms.outputRoot);
end

plotRoot = fullfile(analysisPms.outputRoot, 'plots');
pdfRoot = fullfile(plotRoot, 'PDF');
figRoot = fullfile(plotRoot, 'FIG');
if ~exist(pdfRoot, 'dir')
    mkdir(pdfRoot);
end
if ~exist(figRoot, 'dir')
    mkdir(figRoot);
end

jointRowsPath = fullfile(analysisPms.inputRoot, 'BehaviourNeural_jointRows.csv');
if ~exist(jointRowsPath, 'file')
    error('%s:MissingInput', mfilename, 'Missing %s', jointRowsPath);
end

jointTable = readtable(jointRowsPath);
jointTable = local_add_direction_features(jointTable);
jointTable = local_add_interaction_features(jointTable);

[resultTable, predictionTable] = local_run_regression_grid(jointTable, analysisPms);
topResultTable = local_build_top_result_table(resultTable);

writetable(resultTable, fullfile(analysisPms.outputRoot, ...
    'BehaviourNeural_regressionResults.csv'));
writetable(predictionTable, fullfile(analysisPms.outputRoot, ...
    'BehaviourNeural_regressionPredictions.csv'));
writetable(topResultTable, fullfile(analysisPms.outputRoot, ...
    'BehaviourNeural_regressionTopResults.csv'));

local_plot_top_results(topResultTable, pdfRoot, figRoot, analysisPms);
local_plot_primary_results(resultTable, pdfRoot, figRoot, analysisPms);

save(fullfile(analysisPms.outputRoot, 'BehaviourNeural_regression.mat'), ...
    'analysisPms', 'resultTable', 'predictionTable', 'topResultTable', '-v7.3');

fprintf('Behaviour-neural regression analysis saved in %s\n', analysisPms.outputRoot);


function featureRows = local_feature_rows()
jointWeightNames = {'JointWeightSToK', 'JointWeightKToS', 'NeuralBalanceSToKMinusKToS'};
jointGainNames = {'JointGainSToK', 'JointGainKToS', 'NeuralJointGainBalance'};
directionNames = {'DirectionSin', 'DirectionCos'};
allNeuralNames = [jointWeightNames, jointGainNames];
weightInteractionNames = local_interaction_names(jointWeightNames);
gainInteractionNames = local_interaction_names(jointGainNames);
allInteractionNames = local_interaction_names(allNeuralNames);

featureRows = cell2table({ ...
    'DirectionOnly',           directionNames, ...
    'Direction-only baseline'; ...
    'JointWeights',            jointWeightNames, ...
    'Joint-action crossed weights'; ...
    'JointWeightsPlusDirection', [jointWeightNames, directionNames], ...
    'Joint weights plus direction'; ...
    'JointWeightsByDirection', [jointWeightNames, directionNames, weightInteractionNames], ...
    'Joint weights with direction interactions'; ...
    'JointGains',              jointGainNames, ...
    'Joint-minus-active gain readouts'; ...
    'JointGainsPlusDirection', [jointGainNames, directionNames], ...
    'Joint gains plus direction'; ...
    'JointGainsByDirection', [jointGainNames, directionNames, gainInteractionNames], ...
    'Joint gains with direction interactions'; ...
    'AllNeural',               allNeuralNames, ...
    'All DCM readouts'; ...
    'AllNeuralPlusDirection',  [allNeuralNames, directionNames], ...
    'All DCM readouts plus direction baseline'; ...
    'AllNeuralByDirection',    [allNeuralNames, directionNames, allInteractionNames], ...
    'All DCM readouts with direction interactions'}, ...
    'VariableNames', {'FeatureSet', 'FeatureNames', 'Description'});
end


function interactionNames = local_interaction_names(sourceNames)
interactionNames = cell(1, numel(sourceNames) .* 2);
interactionCounter = 0;
for sourceCounter = 1:numel(sourceNames)
    interactionCounter = interactionCounter + 1;
    interactionNames{interactionCounter} = [sourceNames{sourceCounter}, '_x_DirectionSin'];
    interactionCounter = interactionCounter + 1;
    interactionNames{interactionCounter} = [sourceNames{sourceCounter}, '_x_DirectionCos'];
end
end


function jointTable = local_add_direction_features(jointTable)
directionAngles = 2 .* pi .* (jointTable.Direction - 1) ./ 8;
jointTable.DirectionSin = sin(directionAngles);
jointTable.DirectionCos = cos(directionAngles);
end


function jointTable = local_add_interaction_features(jointTable)
sourceNames = {'JointWeightSToK', 'JointWeightKToS', 'NeuralBalanceSToKMinusKToS', ...
    'JointGainSToK', 'JointGainKToS', 'NeuralJointGainBalance'};

for sourceCounter = 1:numel(sourceNames)
    sourceName = sourceNames{sourceCounter};
    jointTable.([sourceName, '_x_DirectionSin']) = jointTable.(sourceName) .* jointTable.DirectionSin;
    jointTable.([sourceName, '_x_DirectionCos']) = jointTable.(sourceName) .* jointTable.DirectionCos;
end
end


function [resultTable, predictionTable] = local_run_regression_grid(jointTable, analysisPms)
resultRows = {};
predictionRows = {};

for modelRowCounter = 1:height(analysisPms.modelRows)
    modelRow = analysisPms.modelRows(modelRowCounter, :);
    modelMask = strcmp(jointTable.Chamber, modelRow.Chamber{1}) & ...
        strcmp(jointTable.Model, modelRow.Model{1});
    modelTable = jointTable(modelMask, :);

    for targetCounter = 1:height(analysisPms.targetRows)
        targetRow = analysisPms.targetRows(targetCounter, :);

        for featureCounter = 1:height(analysisPms.featureRows)
            featureRow = analysisPms.featureRows(featureCounter, :);
            featureNames = featureRow.FeatureNames{1};
            [featureMatrix, targetValues, sessionLabels, directionVector, validTable] = ...
                local_prepare_regression_data(modelTable, targetRow.TargetVariable{1}, featureNames);

            [resultRow, samplePredictions] = local_evaluate_regression( ...
                featureMatrix, targetValues, sessionLabels, directionVector, ...
                modelRow, targetRow, featureRow, validTable, analysisPms);

            resultRows(end + 1, :) = resultRow; %#ok<AGROW>
            predictionRows = [predictionRows; samplePredictions]; %#ok<AGROW>
        end
    end
end

resultTable = cell2table(resultRows, 'VariableNames', local_result_columns());
predictionTable = cell2table(predictionRows, 'VariableNames', local_prediction_columns());
end


function [featureMatrix, targetValues, sessionLabels, directionVector, validTable] = ...
    local_prepare_regression_data(modelTable, targetName, featureNames)
targetValues = modelTable.(targetName);
featureMatrix = table2array(modelTable(:, featureNames));
validMask = isfinite(targetValues) & all(isfinite(featureMatrix), 2);

validTable = modelTable(validMask, :);
featureMatrix = featureMatrix(validMask, :);
targetValues = targetValues(validMask);
sessionLabels = validTable.Session;
directionVector = validTable.Direction;
end


function [resultRow, predictionRows] = local_evaluate_regression( ...
    featureMatrix, targetValues, sessionLabels, directionVector, ...
    modelRow, targetRow, featureRow, validTable, analysisPms)

sampleCount = numel(targetValues);
sessionCount = numel(unique(sessionLabels, 'stable'));
isValidProblem = sampleCount >= analysisPms.minimumSamples && ...
    sessionCount >= analysisPms.minimumSessions && ...
    numel(unique(targetValues)) > 1;

if isValidProblem
    [predictedValues, foldLabels, selectedLambdas] = local_leave_one_session_out_regression( ...
        featureMatrix, targetValues, sessionLabels, analysisPms);
    metricStruct = local_regression_metrics(targetValues, predictedValues);
else
    predictedValues = nan(sampleCount, 1);
    foldLabels = strings(sampleCount, 1);
    selectedLambdas = nan(sampleCount, 1);
    metricStruct = local_empty_metrics();
end

featureNames = featureRow.FeatureNames{1};
resultRow = { ...
    modelRow.Chamber{1}, modelRow.Model{1}, modelRow.Description{1}, ...
    targetRow.TargetId{1}, targetRow.TargetVariable{1}, targetRow.Description{1}, ...
    featureRow.FeatureSet{1}, strjoin(featureNames, '+'), featureRow.Description{1}, ...
    sampleCount, sessionCount, metricStruct.PredictiveR2, metricStruct.PearsonR, ...
    metricStruct.MAE, metricStruct.RMSE, metricStruct.BaselineRMSE, ...
    median(selectedLambdas, 'omitnan'), isValidProblem};

predictionRows = cell(sampleCount, numel(local_prediction_columns()));
for sampleCounter = 1:sampleCount
    predictionRows(sampleCounter, :) = { ...
        modelRow.Chamber{1}, modelRow.Model{1}, targetRow.TargetId{1}, ...
        targetRow.TargetVariable{1}, featureRow.FeatureSet{1}, ...
        validTable.Session{sampleCounter}, directionVector(sampleCounter), ...
        targetValues(sampleCounter), predictedValues(sampleCounter), ...
        targetValues(sampleCounter) - predictedValues(sampleCounter), ...
        selectedLambdas(sampleCounter), string(foldLabels(sampleCounter))};
end
end


function [predictedValues, foldLabels, selectedLambdas] = local_leave_one_session_out_regression( ...
    featureMatrix, targetValues, sessionLabels, analysisPms)
sampleCount = numel(targetValues);
predictedValues = nan(sampleCount, 1);
selectedLambdas = nan(sampleCount, 1);
foldLabels = strings(sampleCount, 1);
uniqueSessions = unique(sessionLabels, 'stable');

for foldCounter = 1:numel(uniqueSessions)
    testMask = strcmp(sessionLabels, uniqueSessions{foldCounter});
    trainMask = ~testMask;
    trainFeatures = featureMatrix(trainMask, :);
    testFeatures = featureMatrix(testMask, :);
    trainTargets = targetValues(trainMask);
    trainSessions = sessionLabels(trainMask);

    selectedLambda = local_choose_ridge_lambda(trainFeatures, trainTargets, trainSessions, analysisPms);
    predictedValues(testMask) = local_predict_ridge_linear( ...
        trainFeatures, trainTargets, testFeatures, selectedLambda);
    selectedLambdas(testMask) = selectedLambda;
    foldLabels(testMask) = string(uniqueSessions{foldCounter});
end
end


function selectedLambda = local_choose_ridge_lambda(trainFeatures, trainTargets, trainSessions, analysisPms)
uniqueTrainSessions = unique(trainSessions, 'stable');
if numel(uniqueTrainSessions) < 3
    selectedLambda = analysisPms.defaultRidgeLambda;
    return;
end

meanSquaredErrors = nan(numel(analysisPms.ridgeLambdas), 1);
for lambdaCounter = 1:numel(analysisPms.ridgeLambdas)
    lambdaValue = analysisPms.ridgeLambdas(lambdaCounter);
    innerPredictions = nan(numel(trainTargets), 1);

    for foldCounter = 1:numel(uniqueTrainSessions)
        innerTestMask = strcmp(trainSessions, uniqueTrainSessions{foldCounter});
        innerTrainMask = ~innerTestMask;
        innerPredictions(innerTestMask) = local_predict_ridge_linear( ...
            trainFeatures(innerTrainMask, :), trainTargets(innerTrainMask), ...
            trainFeatures(innerTestMask, :), lambdaValue);
    end

    residualValues = trainTargets - innerPredictions;
    meanSquaredErrors(lambdaCounter) = mean(residualValues .^ 2, 'omitnan');
end

[~, bestLambdaIndex] = min(meanSquaredErrors);
selectedLambda = analysisPms.ridgeLambdas(bestLambdaIndex);
end


function predictedValues = local_predict_ridge_linear(trainFeatures, trainTargets, testFeatures, lambdaValue)
[trainFeatures, testFeatures] = local_standardize_from_train(trainFeatures, testFeatures);

targetMean = mean(trainTargets, 'omitnan');
targetScale = std(trainTargets, 'omitnan');
if targetScale == 0 || ~isfinite(targetScale)
    predictedValues = repmat(targetMean, size(testFeatures, 1), 1);
    return;
end

trainTargetsStandardized = (trainTargets - targetMean) ./ targetScale;
trainDesign = [ones(size(trainFeatures, 1), 1), trainFeatures];
testDesign = [ones(size(testFeatures, 1), 1), testFeatures];
ridgePenalty = lambdaValue .* eye(size(trainDesign, 2));
ridgePenalty(1, 1) = 0;

coefficients = (trainDesign' * trainDesign + ridgePenalty) \ (trainDesign' * trainTargetsStandardized);
predictedStandardized = testDesign * coefficients;
predictedValues = predictedStandardized .* targetScale + targetMean;
end


function [trainFeatures, testFeatures] = local_standardize_from_train(trainFeatures, testFeatures)
featureMean = mean(trainFeatures, 1, 'omitnan');
featureScale = std(trainFeatures, 0, 1, 'omitnan');
featureScale(featureScale == 0 | ~isfinite(featureScale)) = 1;
trainFeatures = (trainFeatures - featureMean) ./ featureScale;
testFeatures = (testFeatures - featureMean) ./ featureScale;
end


function metricStruct = local_regression_metrics(targetValues, predictedValues)
validMask = isfinite(targetValues) & isfinite(predictedValues);
targetValues = targetValues(validMask);
predictedValues = predictedValues(validMask);

if isempty(targetValues)
    metricStruct = local_empty_metrics();
    return;
end

residualValues = targetValues - predictedValues;
totalSumSquares = sum((targetValues - mean(targetValues, 'omitnan')) .^ 2);
residualSumSquares = sum(residualValues .^ 2);
if totalSumSquares == 0
    predictiveR2 = NaN;
else
    predictiveR2 = 1 - residualSumSquares ./ totalSumSquares;
end

metricStruct = struct();
metricStruct.PredictiveR2 = predictiveR2;
metricStruct.PearsonR = local_pearson_correlation(targetValues, predictedValues);
metricStruct.MAE = mean(abs(residualValues), 'omitnan');
metricStruct.RMSE = sqrt(mean(residualValues .^ 2, 'omitnan'));
metricStruct.BaselineRMSE = sqrt(mean((targetValues - mean(targetValues, 'omitnan')) .^ 2, 'omitnan'));
end


function metricStruct = local_empty_metrics()
metricStruct = struct();
metricStruct.PredictiveR2 = NaN;
metricStruct.PearsonR = NaN;
metricStruct.MAE = NaN;
metricStruct.RMSE = NaN;
metricStruct.BaselineRMSE = NaN;
end


function correlationValue = local_pearson_correlation(firstValues, secondValues)
validMask = isfinite(firstValues) & isfinite(secondValues);
firstValues = firstValues(validMask);
secondValues = secondValues(validMask);

if numel(firstValues) < 3 || std(firstValues) == 0 || std(secondValues) == 0
    correlationValue = NaN;
    return;
end

firstValues = firstValues - mean(firstValues);
secondValues = secondValues - mean(secondValues);
correlationValue = sum(firstValues .* secondValues) ./ ...
    sqrt(sum(firstValues .^ 2) .* sum(secondValues .^ 2));
end


function topResultTable = local_build_top_result_table(resultTable)
validMask = resultTable.IsValidProblem & isfinite(resultTable.PredictiveR2);
topResultTable = resultTable(validMask, :);
topResultTable = sortrows(topResultTable, {'PredictiveR2', 'PearsonR'}, {'descend', 'descend'});
if height(topResultTable) > 40
    topResultTable = topResultTable(1:40, :);
end
end


function local_plot_top_results(topResultTable, pdfRoot, figRoot, analysisPms)
if isempty(topResultTable)
    return;
end

plotTable = topResultTable(1:min(20, height(topResultTable)), :);
figureHandle = figure('Visible', analysisPms.visible, 'Color', 'w', ...
    'Position', [100 100 950 720]);
barh(plotTable.PredictiveR2, 'FaceColor', [0.20 0.45 0.70], 'EdgeColor', 'none');
set(gca, 'YDir', 'reverse');
set(gca, 'YTick', 1:height(plotTable));
yticklabels(local_result_labels(plotTable));
xline(0, 'k-', 'LineWidth', 1);
xlabel('Leave-one-session-out predictive R^2');
title('Top continuous behaviour-neural regression checks');
grid on;
box off;

local_save_figure(figureHandle, pdfRoot, figRoot, 'BehaviourNeural_Regression_TopResults');
close(figureHandle);
end


function local_plot_primary_results(resultTable, pdfRoot, figRoot, analysisPms)
primaryRows = cell2table({ ...
    'Frontal',  'M7', 'S_AE_Benefit',   'Frontal M7: S |AE| benefit'; ...
    'Frontal',  'M7', 'S_RT_Advantage', 'Frontal M7: S RT advantage'; ...
    'Parietal', 'M2', 'S_AE_Benefit',   'Parietal M2: S |AE| benefit'; ...
    'Parietal', 'M2', 'S_RT_Advantage', 'Parietal M2: S RT advantage'}, ...
    'VariableNames', {'Chamber', 'Model', 'TargetId', 'TitleText'});

figureHandle = figure('Visible', analysisPms.visible, 'Color', 'w', ...
    'Position', [100 100 1050 760]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for primaryCounter = 1:height(primaryRows)
    nexttile;
    primaryRow = primaryRows(primaryCounter, :);
    plotMask = strcmp(resultTable.Chamber, primaryRow.Chamber{1}) & ...
        strcmp(resultTable.Model, primaryRow.Model{1}) & ...
        strcmp(resultTable.TargetId, primaryRow.TargetId{1}) & ...
        resultTable.IsValidProblem;
    plotTable = resultTable(plotMask, :);
    plotTable = sortrows(plotTable, 'PredictiveR2', 'descend');

    barh(plotTable.PredictiveR2, 'FaceColor', [0.28 0.55 0.38], 'EdgeColor', 'none');
    set(gca, 'YDir', 'reverse');
    set(gca, 'YTick', 1:height(plotTable));
    yticklabels(plotTable.FeatureSet);
    xline(0, 'k-', 'LineWidth', 1);
    xlabel('Predictive R^2');
    title(primaryRow.TitleText{1}, 'Interpreter', 'none');
    grid on;
    box off;
end

local_save_figure(figureHandle, pdfRoot, figRoot, 'BehaviourNeural_Regression_PrimaryModels');
close(figureHandle);
end


function resultLabels = local_result_labels(resultTable)
resultLabels = strings(height(resultTable), 1);
for rowCounter = 1:height(resultTable)
    resultLabels(rowCounter) = sprintf('%s %s | %s | %s', ...
        resultTable.Chamber{rowCounter}, resultTable.Model{rowCounter}, ...
        resultTable.TargetId{rowCounter}, resultTable.FeatureSet{rowCounter});
end
end


function local_save_figure(figureHandle, pdfRoot, figRoot, fileStem)
pdfPath = fullfile(pdfRoot, [fileStem, '.pdf']);
figPath = fullfile(figRoot, [fileStem, '.fig']);
set(figureHandle, 'PaperPositionMode', 'auto');
print(figureHandle, pdfPath, '-dpdf', '-bestfit');
savefig(figureHandle, figPath);
end


function columnNames = local_result_columns()
columnNames = {'Chamber', 'Model', 'ModelDescription', ...
    'TargetId', 'TargetVariable', 'TargetDescription', ...
    'FeatureSet', 'FeatureNames', 'FeatureDescription', ...
    'N', 'SessionCount', 'PredictiveR2', 'PearsonR', ...
    'MAE', 'RMSE', 'BaselineRMSE', 'MedianSelectedLambda', 'IsValidProblem'};
end


function columnNames = local_prediction_columns()
columnNames = {'Chamber', 'Model', 'TargetId', 'TargetVariable', 'FeatureSet', ...
    'Session', 'Direction', 'TargetValue', 'PredictedValue', 'ResidualValue', ...
    'SelectedLambda', 'FoldSession'};
end
