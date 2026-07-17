%% MAIN_SAPIENZA_BEHAVIOUR_NEURAL_DIRECTIONAL_PREDICTION
% Direction-wise continuous regression and high/low classification checks.
%
% This script mirrors the direction-resolved correlation plots more closely
% than the pooled regression/classification scripts.  For each fixed direction,
% one sample is one session, and validation leaves one session out.

clear;
close all;

analysisPms = struct();
analysisPms.inputRoot = '/TESTS/SAPIENZA/BEHAVIOUR/NeuralCorrelation';
analysisPms.outputRoot = '/TESTS/SAPIENZA/BEHAVIOUR/NeuralDirectionalPrediction';
analysisPms.visible = 'off';
analysisPms.minimumSamples = 12;
analysisPms.minimumSamplesPerClass = 5;
analysisPms.ridgeLambda = 1;
analysisPms.maximumLogisticIterations = 80;
analysisPms.logisticTolerance = 1e-7;
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
analysisPms.featureRows = cell2table({ ...
    'JointWeights', {'JointWeightSToK', 'JointWeightKToS', 'NeuralBalanceSToKMinusKToS'}, ...
    'Joint-action crossed weights'; ...
    'JointGains', {'JointGainSToK', 'JointGainKToS', 'NeuralJointGainBalance'}, ...
    'Joint-minus-active gain readouts'; ...
    'AllNeural', {'JointWeightSToK', 'JointWeightKToS', 'NeuralBalanceSToKMinusKToS', ...
                  'JointGainSToK', 'JointGainKToS', 'NeuralJointGainBalance'}, ...
    'All DCM readouts'}, ...
    'VariableNames', {'FeatureSet', 'FeatureNames', 'Description'});

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

[regressionResultTable, regressionPredictionTable, classificationResultTable, classificationPredictionTable] = ...
    local_run_directional_prediction_grid(jointTable, analysisPms);

topRegressionTable = local_top_regression_table(regressionResultTable);
topClassificationTable = local_top_classification_table(classificationResultTable);

writetable(regressionResultTable, fullfile(analysisPms.outputRoot, ...
    'BehaviourNeural_directionalRegressionResults.csv'));
writetable(regressionPredictionTable, fullfile(analysisPms.outputRoot, ...
    'BehaviourNeural_directionalRegressionPredictions.csv'));
writetable(topRegressionTable, fullfile(analysisPms.outputRoot, ...
    'BehaviourNeural_directionalRegressionTopResults.csv'));
writetable(classificationResultTable, fullfile(analysisPms.outputRoot, ...
    'BehaviourNeural_directionalClassificationResults.csv'));
writetable(classificationPredictionTable, fullfile(analysisPms.outputRoot, ...
    'BehaviourNeural_directionalClassificationPredictions.csv'));
writetable(topClassificationTable, fullfile(analysisPms.outputRoot, ...
    'BehaviourNeural_directionalClassificationTopResults.csv'));

local_plot_top_regression(topRegressionTable, pdfRoot, figRoot, analysisPms);
local_plot_top_classification(topClassificationTable, pdfRoot, figRoot, analysisPms);
local_plot_primary_direction_panels(regressionResultTable, classificationResultTable, ...
    pdfRoot, figRoot, analysisPms);

save(fullfile(analysisPms.outputRoot, 'BehaviourNeural_directionalPrediction.mat'), ...
    'analysisPms', 'regressionResultTable', 'regressionPredictionTable', ...
    'classificationResultTable', 'classificationPredictionTable', ...
    'topRegressionTable', 'topClassificationTable', '-v7.3');

fprintf('Behaviour-neural directional prediction analysis saved in %s\n', analysisPms.outputRoot);


function [regressionResultTable, regressionPredictionTable, classificationResultTable, classificationPredictionTable] = ...
    local_run_directional_prediction_grid(jointTable, analysisPms)
regressionResultRows = {};
regressionPredictionRows = {};
classificationResultRows = {};
classificationPredictionRows = {};

for modelRowCounter = 1:height(analysisPms.modelRows)
    modelRow = analysisPms.modelRows(modelRowCounter, :);
    modelMask = strcmp(jointTable.Chamber, modelRow.Chamber{1}) & ...
        strcmp(jointTable.Model, modelRow.Model{1});
    modelTable = jointTable(modelMask, :);

    for targetCounter = 1:height(analysisPms.targetRows)
        targetRow = analysisPms.targetRows(targetCounter, :);

        for featureCounter = 1:height(analysisPms.featureRows)
            featureRow = analysisPms.featureRows(featureCounter, :);

            for directionCounter = 1:8
                directionTable = modelTable(modelTable.Direction == directionCounter, :);
                featureNames = featureRow.FeatureNames{1};
                [featureMatrix, targetValues, sessionLabels, validTable] = ...
                    local_prepare_data(directionTable, targetRow.TargetVariable{1}, featureNames);

                [regressionRow, regressionPredictions] = local_evaluate_directional_regression( ...
                    featureMatrix, targetValues, sessionLabels, modelRow, targetRow, ...
                    featureRow, directionCounter, validTable, analysisPms);
                [classificationRow, classificationPredictions] = local_evaluate_directional_classification( ...
                    featureMatrix, targetValues, sessionLabels, modelRow, targetRow, ...
                    featureRow, directionCounter, validTable, analysisPms);

                regressionResultRows(end + 1, :) = regressionRow; %#ok<AGROW>
                regressionPredictionRows = [regressionPredictionRows; regressionPredictions]; %#ok<AGROW>
                classificationResultRows(end + 1, :) = classificationRow; %#ok<AGROW>
                classificationPredictionRows = [classificationPredictionRows; classificationPredictions]; %#ok<AGROW>
            end
        end
    end
end

regressionResultTable = cell2table(regressionResultRows, ...
    'VariableNames', local_regression_result_columns());
regressionPredictionTable = cell2table(regressionPredictionRows, ...
    'VariableNames', local_regression_prediction_columns());
classificationResultTable = cell2table(classificationResultRows, ...
    'VariableNames', local_classification_result_columns());
classificationPredictionTable = cell2table(classificationPredictionRows, ...
    'VariableNames', local_classification_prediction_columns());
end


function [featureMatrix, targetValues, sessionLabels, validTable] = ...
    local_prepare_data(directionTable, targetName, featureNames)
targetValues = directionTable.(targetName);
featureMatrix = table2array(directionTable(:, featureNames));
validMask = isfinite(targetValues) & all(isfinite(featureMatrix), 2);

validTable = directionTable(validMask, :);
featureMatrix = featureMatrix(validMask, :);
targetValues = targetValues(validMask);
sessionLabels = validTable.Session;
end


function [resultRow, predictionRows] = local_evaluate_directional_regression( ...
    featureMatrix, targetValues, sessionLabels, modelRow, targetRow, featureRow, ...
    directionNumber, validTable, analysisPms)
sampleCount = numel(targetValues);
sessionCount = numel(unique(sessionLabels, 'stable'));
isValidProblem = sampleCount >= analysisPms.minimumSamples && ...
    sessionCount >= analysisPms.minimumSamples && numel(unique(targetValues)) > 1;

if isValidProblem
    [predictedValues, baselinePredictedValues, foldLabels] = local_leave_one_session_out_regression( ...
        featureMatrix, targetValues, sessionLabels, analysisPms.ridgeLambda);
    metricStruct = local_regression_metrics(targetValues, predictedValues, baselinePredictedValues);
else
    predictedValues = nan(sampleCount, 1);
    baselinePredictedValues = nan(sampleCount, 1);
    foldLabels = strings(sampleCount, 1);
    metricStruct = local_empty_regression_metrics();
end

featureNames = featureRow.FeatureNames{1};
resultRow = { ...
    modelRow.Chamber{1}, modelRow.Model{1}, modelRow.Description{1}, ...
    targetRow.TargetId{1}, targetRow.TargetVariable{1}, targetRow.Description{1}, ...
    featureRow.FeatureSet{1}, strjoin(featureNames, '+'), featureRow.Description{1}, ...
    directionNumber, sampleCount, sessionCount, ...
    metricStruct.PredictiveR2, metricStruct.PearsonR, metricStruct.MAE, ...
    metricStruct.RMSE, metricStruct.BaselineRMSE, isValidProblem};

predictionRows = cell(sampleCount, numel(local_regression_prediction_columns()));
for sampleCounter = 1:sampleCount
    predictionRows(sampleCounter, :) = { ...
        modelRow.Chamber{1}, modelRow.Model{1}, targetRow.TargetId{1}, ...
        targetRow.TargetVariable{1}, featureRow.FeatureSet{1}, directionNumber, ...
        validTable.Session{sampleCounter}, targetValues(sampleCounter), ...
        predictedValues(sampleCounter), baselinePredictedValues(sampleCounter), ...
        targetValues(sampleCounter) - predictedValues(sampleCounter), ...
        string(foldLabels(sampleCounter))};
end
end


function [predictedValues, baselinePredictedValues, foldLabels] = local_leave_one_session_out_regression( ...
    featureMatrix, targetValues, sessionLabels, ridgeLambda)
sampleCount = numel(targetValues);
predictedValues = nan(sampleCount, 1);
baselinePredictedValues = nan(sampleCount, 1);
foldLabels = strings(sampleCount, 1);
uniqueSessions = unique(sessionLabels, 'stable');

for foldCounter = 1:numel(uniqueSessions)
    testMask = strcmp(sessionLabels, uniqueSessions{foldCounter});
    trainMask = ~testMask;
    trainFeatures = featureMatrix(trainMask, :);
    testFeatures = featureMatrix(testMask, :);
    trainTargets = targetValues(trainMask);

    predictedValues(testMask) = local_predict_ridge_linear( ...
        trainFeatures, trainTargets, testFeatures, ridgeLambda);
    baselinePredictedValues(testMask) = mean(trainTargets, 'omitnan');
    foldLabels(testMask) = string(uniqueSessions{foldCounter});
end
end


function predictedValues = local_predict_ridge_linear(trainFeatures, trainTargets, testFeatures, ridgeLambda)
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
ridgePenalty = ridgeLambda .* eye(size(trainDesign, 2));
ridgePenalty(1, 1) = 0;
coefficients = (trainDesign' * trainDesign + ridgePenalty) \ (trainDesign' * trainTargetsStandardized);
predictedValues = (testDesign * coefficients) .* targetScale + targetMean;
end


function metricStruct = local_regression_metrics(targetValues, predictedValues, baselinePredictedValues)
validMask = isfinite(targetValues) & isfinite(predictedValues) & isfinite(baselinePredictedValues);
targetValues = targetValues(validMask);
predictedValues = predictedValues(validMask);
baselinePredictedValues = baselinePredictedValues(validMask);

if isempty(targetValues)
    metricStruct = local_empty_regression_metrics();
    return;
end

residualValues = targetValues - predictedValues;
baselineResidualValues = targetValues - baselinePredictedValues;
baselineSumSquares = sum(baselineResidualValues .^ 2);
if baselineSumSquares == 0
    predictiveR2 = NaN;
else
    predictiveR2 = 1 - sum(residualValues .^ 2) ./ baselineSumSquares;
end

metricStruct = struct();
metricStruct.PredictiveR2 = predictiveR2;
metricStruct.PearsonR = local_pearson_correlation(targetValues, predictedValues);
metricStruct.MAE = mean(abs(residualValues), 'omitnan');
metricStruct.RMSE = sqrt(mean(residualValues .^ 2, 'omitnan'));
metricStruct.BaselineRMSE = sqrt(mean(baselineResidualValues .^ 2, 'omitnan'));
end


function metricStruct = local_empty_regression_metrics()
metricStruct = struct();
metricStruct.PredictiveR2 = NaN;
metricStruct.PearsonR = NaN;
metricStruct.MAE = NaN;
metricStruct.RMSE = NaN;
metricStruct.BaselineRMSE = NaN;
end


function [resultRow, predictionRows] = local_evaluate_directional_classification( ...
    featureMatrix, targetValues, sessionLabels, modelRow, targetRow, featureRow, ...
    directionNumber, validTable, analysisPms)
sampleCount = numel(targetValues);
sessionCount = numel(unique(sessionLabels, 'stable'));
classThreshold = median(targetValues, 'omitnan');
validClassMask = targetValues ~= classThreshold;
labelVector = targetValues(validClassMask) > classThreshold;
classFeatureMatrix = featureMatrix(validClassMask, :);
classTargetValues = targetValues(validClassMask);
classSessionLabels = sessionLabels(validClassMask);
classTable = validTable(validClassMask, :);
positiveCount = sum(labelVector);
negativeCount = numel(labelVector) - positiveCount;
isValidProblem = numel(labelVector) >= analysisPms.minimumSamples && ...
    positiveCount >= analysisPms.minimumSamplesPerClass && ...
    negativeCount >= analysisPms.minimumSamplesPerClass && ...
    numel(unique(classSessionLabels, 'stable')) >= analysisPms.minimumSamples;

if isValidProblem
    [predictedLabels, predictedProbability, foldLabels] = local_leave_one_session_out_classification( ...
        classFeatureMatrix, labelVector, classSessionLabels, analysisPms);
    metricStruct = local_classification_metrics(labelVector, predictedLabels, predictedProbability);
else
    predictedLabels = false(numel(labelVector), 1);
    predictedProbability = nan(numel(labelVector), 1);
    foldLabels = strings(numel(labelVector), 1);
    metricStruct = local_empty_classification_metrics();
end

featureNames = featureRow.FeatureNames{1};
resultRow = { ...
    modelRow.Chamber{1}, modelRow.Model{1}, modelRow.Description{1}, ...
    targetRow.TargetId{1}, targetRow.TargetVariable{1}, targetRow.Description{1}, ...
    featureRow.FeatureSet{1}, strjoin(featureNames, '+'), featureRow.Description{1}, ...
    directionNumber, sampleCount, sessionCount, classThreshold, ...
    positiveCount, negativeCount, metricStruct.Accuracy, metricStruct.BalancedAccuracy, ...
    metricStruct.Sensitivity, metricStruct.Specificity, metricStruct.AUC, isValidProblem};

predictionRows = cell(numel(labelVector), numel(local_classification_prediction_columns()));
for sampleCounter = 1:numel(labelVector)
    predictionRows(sampleCounter, :) = { ...
        modelRow.Chamber{1}, modelRow.Model{1}, targetRow.TargetId{1}, ...
        targetRow.TargetVariable{1}, featureRow.FeatureSet{1}, directionNumber, ...
        classTable.Session{sampleCounter}, classTargetValues(sampleCounter), ...
        classThreshold, labelVector(sampleCounter), predictedLabels(sampleCounter), ...
        predictedProbability(sampleCounter), string(foldLabels(sampleCounter))};
end
end


function [predictedLabels, predictedProbability, foldLabels] = local_leave_one_session_out_classification( ...
    featureMatrix, labelVector, sessionLabels, analysisPms)
sampleCount = numel(labelVector);
predictedLabels = false(sampleCount, 1);
predictedProbability = nan(sampleCount, 1);
foldLabels = strings(sampleCount, 1);
uniqueSessions = unique(sessionLabels, 'stable');

for foldCounter = 1:numel(uniqueSessions)
    testMask = strcmp(sessionLabels, uniqueSessions{foldCounter});
    trainMask = ~testMask;
    trainFeatures = featureMatrix(trainMask, :);
    testFeatures = featureMatrix(testMask, :);
    trainLabels = labelVector(trainMask);

    if numel(unique(trainLabels)) < 2
        probability = repmat(mean(trainLabels), sum(testMask), 1);
    else
        probability = local_predict_ridge_logistic( ...
            trainFeatures, trainLabels, testFeatures, analysisPms);
    end

    predictedProbability(testMask) = probability;
    predictedLabels(testMask) = probability >= 0.5;
    foldLabels(testMask) = string(uniqueSessions{foldCounter});
end
end


function probability = local_predict_ridge_logistic(trainFeatures, trainLabels, testFeatures, analysisPms)
[trainFeatures, testFeatures] = local_standardize_from_train(trainFeatures, testFeatures);
trainDesign = [ones(size(trainFeatures, 1), 1), trainFeatures];
testDesign = [ones(size(testFeatures, 1), 1), testFeatures];
coefficientCount = size(trainDesign, 2);
coefficients = zeros(coefficientCount, 1);
ridgeWeights = analysisPms.ridgeLambda .* ones(coefficientCount, 1);
ridgeWeights(1) = 0;

for iterationCounter = 1:analysisPms.maximumLogisticIterations
    linearPredictor = trainDesign * coefficients;
    fittedProbability = local_sigmoid(linearPredictor);
    gradient = trainDesign' * (fittedProbability - trainLabels) + ridgeWeights .* coefficients;
    weightVector = max(fittedProbability .* (1 - fittedProbability), 1e-6);
    hessian = trainDesign' * (trainDesign .* weightVector) + diag(ridgeWeights);
    updateStep = hessian \ gradient;
    coefficients = coefficients - updateStep;

    if norm(updateStep) < analysisPms.logisticTolerance
        break;
    end
end

probability = local_sigmoid(testDesign * coefficients);
end


function probability = local_sigmoid(linearPredictor)
probability = 1 ./ (1 + exp(-max(min(linearPredictor, 50), -50)));
end


function [trainFeatures, testFeatures] = local_standardize_from_train(trainFeatures, testFeatures)
featureMean = mean(trainFeatures, 1, 'omitnan');
featureScale = std(trainFeatures, 0, 1, 'omitnan');
featureScale(featureScale == 0 | ~isfinite(featureScale)) = 1;
trainFeatures = (trainFeatures - featureMean) ./ featureScale;
testFeatures = (testFeatures - featureMean) ./ featureScale;
end


function metricStruct = local_classification_metrics(labelVector, predictedLabels, predictedProbability)
truePositive = sum(labelVector & predictedLabels);
trueNegative = sum(~labelVector & ~predictedLabels);
falsePositive = sum(~labelVector & predictedLabels);
falseNegative = sum(labelVector & ~predictedLabels);

sensitivity = truePositive ./ max(truePositive + falseNegative, 1);
specificity = trueNegative ./ max(trueNegative + falsePositive, 1);

metricStruct = struct();
metricStruct.Accuracy = mean(labelVector == predictedLabels);
metricStruct.BalancedAccuracy = 0.5 .* (sensitivity + specificity);
metricStruct.Sensitivity = sensitivity;
metricStruct.Specificity = specificity;
metricStruct.AUC = local_auc(labelVector, predictedProbability);
end


function metricStruct = local_empty_classification_metrics()
metricStruct = struct();
metricStruct.Accuracy = NaN;
metricStruct.BalancedAccuracy = NaN;
metricStruct.Sensitivity = NaN;
metricStruct.Specificity = NaN;
metricStruct.AUC = NaN;
end


function aucValue = local_auc(labelVector, scoreVector)
validMask = isfinite(scoreVector);
labelVector = labelVector(validMask);
scoreVector = scoreVector(validMask);

positiveScores = scoreVector(labelVector);
negativeScores = scoreVector(~labelVector);

if isempty(positiveScores) || isempty(negativeScores)
    aucValue = NaN;
    return;
end

comparisonMatrix = positiveScores(:) - negativeScores(:)';
aucValue = (sum(comparisonMatrix(:) > 0) + 0.5 .* sum(comparisonMatrix(:) == 0)) ./ numel(comparisonMatrix);
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


function topRegressionTable = local_top_regression_table(regressionResultTable)
validMask = regressionResultTable.IsValidProblem & isfinite(regressionResultTable.PredictiveR2);
topRegressionTable = regressionResultTable(validMask, :);
topRegressionTable = sortrows(topRegressionTable, {'PredictiveR2', 'PearsonR'}, {'descend', 'descend'});
if height(topRegressionTable) > 60
    topRegressionTable = topRegressionTable(1:60, :);
end
end


function topClassificationTable = local_top_classification_table(classificationResultTable)
validMask = classificationResultTable.IsValidProblem & isfinite(classificationResultTable.BalancedAccuracy);
topClassificationTable = classificationResultTable(validMask, :);
topClassificationTable = sortrows(topClassificationTable, {'BalancedAccuracy', 'AUC'}, {'descend', 'descend'});
if height(topClassificationTable) > 60
    topClassificationTable = topClassificationTable(1:60, :);
end
end


function local_plot_top_regression(topRegressionTable, pdfRoot, figRoot, analysisPms)
if isempty(topRegressionTable)
    return;
end

plotTable = topRegressionTable(1:min(24, height(topRegressionTable)), :);
figureHandle = figure('Visible', analysisPms.visible, 'Color', 'w', ...
    'Position', [100 100 980 760]);
barh(plotTable.PredictiveR2, 'FaceColor', [0.20 0.45 0.70], 'EdgeColor', 'none');
set(gca, 'YDir', 'reverse');
set(gca, 'YTick', 1:height(plotTable));
yticklabels(local_directional_result_labels(plotTable));
xline(0, 'k-', 'LineWidth', 1);
xlabel('Direction-wise leave-one-session-out predictive R^2');
title('Top direction-wise continuous behaviour-neural regressions');
grid on;
box off;
local_save_figure(figureHandle, pdfRoot, figRoot, 'BehaviourNeural_DirectionalRegression_TopResults');
close(figureHandle);
end


function local_plot_top_classification(topClassificationTable, pdfRoot, figRoot, analysisPms)
if isempty(topClassificationTable)
    return;
end

plotTable = topClassificationTable(1:min(24, height(topClassificationTable)), :);
figureHandle = figure('Visible', analysisPms.visible, 'Color', 'w', ...
    'Position', [100 100 980 760]);
barh(plotTable.BalancedAccuracy, 'FaceColor', [0.72 0.38 0.22], 'EdgeColor', 'none');
set(gca, 'YDir', 'reverse');
set(gca, 'YTick', 1:height(plotTable));
yticklabels(local_directional_result_labels(plotTable));
xline(0.5, 'k--', 'LineWidth', 1);
xlabel('Direction-wise leave-one-session-out balanced accuracy');
title('Top direction-wise high/low behaviour-neural classifications');
grid on;
box off;
local_save_figure(figureHandle, pdfRoot, figRoot, 'BehaviourNeural_DirectionalClassification_TopResults');
close(figureHandle);
end


function local_plot_primary_direction_panels(regressionResultTable, classificationResultTable, ...
    pdfRoot, figRoot, analysisPms)
primaryRows = cell2table({ ...
    'Frontal',  'M7', 'S_AE_Benefit',   'Frontal M7: S |AE| benefit'; ...
    'Frontal',  'M7', 'S_RT_Advantage', 'Frontal M7: S RT advantage'; ...
    'Parietal', 'M2', 'S_AE_Benefit',   'Parietal M2: S |AE| benefit'; ...
    'Parietal', 'M2', 'S_RT_Advantage', 'Parietal M2: S RT advantage'}, ...
    'VariableNames', {'Chamber', 'Model', 'TargetId', 'TitleText'});

featureSets = {'JointWeights', 'JointGains', 'AllNeural'};
featureColors = [0.20 0.45 0.70; 0.28 0.55 0.38; 0.72 0.38 0.22];

figureHandle = figure('Visible', analysisPms.visible, 'Color', 'w', ...
    'Position', [100 100 1120 940]);
tileLayout = tiledlayout(5, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

legendAxes = nexttile(tileLayout, [1 2]);
set(legendAxes, 'Visible', 'off');
hold(legendAxes, 'on');
legendHandles = gobjects(numel(featureSets), 1);
for featureCounter = 1:numel(featureSets)
    legendHandles(featureCounter) = patch(legendAxes, NaN, NaN, ...
        featureColors(featureCounter, :), 'EdgeColor', 'none');
end
legend(legendAxes, legendHandles, featureSets, ...
    'Orientation', 'horizontal', 'Interpreter', 'none', ...
    'Location', 'north');
axis(legendAxes, 'off');

for primaryCounter = 1:height(primaryRows)
    primaryRow = primaryRows(primaryCounter, :);
    regressionMask = strcmp(regressionResultTable.Chamber, primaryRow.Chamber{1}) & ...
        strcmp(regressionResultTable.Model, primaryRow.Model{1}) & ...
        strcmp(regressionResultTable.TargetId, primaryRow.TargetId{1});
    classificationMask = strcmp(classificationResultTable.Chamber, primaryRow.Chamber{1}) & ...
        strcmp(classificationResultTable.Model, primaryRow.Model{1}) & ...
        strcmp(classificationResultTable.TargetId, primaryRow.TargetId{1});

    nexttile;
    local_plot_direction_metric(regressionResultTable(regressionMask, :), featureSets, ...
        'PredictiveR2', featureColors, primaryRow.TitleText{1}, 'Predictive R^2', 0);

    nexttile;
    local_plot_direction_metric(classificationResultTable(classificationMask, :), featureSets, ...
        'BalancedAccuracy', featureColors, primaryRow.TitleText{1}, 'Balanced accuracy', 0.5);
end

local_save_figure(figureHandle, pdfRoot, figRoot, 'BehaviourNeural_DirectionalPrediction_PrimaryModels');
close(figureHandle);
end


function local_plot_direction_metric(plotTable, featureSets, metricName, featureColors, titleText, yLabelText, baselineValue)
metricMatrix = nan(8, numel(featureSets));
for featureCounter = 1:numel(featureSets)
    featureMask = strcmp(plotTable.FeatureSet, featureSets{featureCounter});
    featureTable = plotTable(featureMask, :);
    for directionCounter = 1:8
        directionMask = featureTable.Direction == directionCounter;
        if any(directionMask)
            metricMatrix(directionCounter, featureCounter) = featureTable.(metricName)(find(directionMask, 1, 'first'));
        end
    end
end

barHandle = bar(metricMatrix, 'grouped');
for featureCounter = 1:numel(featureSets)
    barHandle(featureCounter).FaceColor = featureColors(featureCounter, :);
    barHandle(featureCounter).EdgeColor = 'none';
end
xlineValues = 1:8; %#ok<NASGU>
yline(baselineValue, 'k--', 'LineWidth', 1);
set(gca, 'XTick', 1:8, 'XTickLabel', {'D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7', 'D8'});
ylabel(yLabelText);
title(titleText, 'Interpreter', 'none');
grid on;
box off;
end


function resultLabels = local_directional_result_labels(resultTable)
resultLabels = strings(height(resultTable), 1);
for rowCounter = 1:height(resultTable)
    resultLabels(rowCounter) = sprintf('%s %s D%d | %s | %s', ...
        resultTable.Chamber{rowCounter}, resultTable.Model{rowCounter}, ...
        resultTable.Direction(rowCounter), resultTable.TargetId{rowCounter}, ...
        resultTable.FeatureSet{rowCounter});
end
end


function local_save_figure(figureHandle, pdfRoot, figRoot, fileStem)
pdfPath = fullfile(pdfRoot, [fileStem, '.pdf']);
figPath = fullfile(figRoot, [fileStem, '.fig']);
set(figureHandle, 'PaperPositionMode', 'auto');
print(figureHandle, pdfPath, '-dpdf', '-bestfit');
savefig(figureHandle, figPath);
end


function columnNames = local_regression_result_columns()
columnNames = {'Chamber', 'Model', 'ModelDescription', ...
    'TargetId', 'TargetVariable', 'TargetDescription', ...
    'FeatureSet', 'FeatureNames', 'FeatureDescription', ...
    'Direction', 'N', 'SessionCount', 'PredictiveR2', 'PearsonR', ...
    'MAE', 'RMSE', 'BaselineRMSE', 'IsValidProblem'};
end


function columnNames = local_regression_prediction_columns()
columnNames = {'Chamber', 'Model', 'TargetId', 'TargetVariable', 'FeatureSet', ...
    'Direction', 'Session', 'TargetValue', 'PredictedValue', ...
    'BaselinePredictedValue', 'ResidualValue', 'FoldSession'};
end


function columnNames = local_classification_result_columns()
columnNames = {'Chamber', 'Model', 'ModelDescription', ...
    'TargetId', 'TargetVariable', 'TargetDescription', ...
    'FeatureSet', 'FeatureNames', 'FeatureDescription', ...
    'Direction', 'N', 'SessionCount', 'ClassThreshold', ...
    'PositiveCount', 'NegativeCount', 'Accuracy', 'BalancedAccuracy', ...
    'Sensitivity', 'Specificity', 'AUC', 'IsValidProblem'};
end


function columnNames = local_classification_prediction_columns()
columnNames = {'Chamber', 'Model', 'TargetId', 'TargetVariable', 'FeatureSet', ...
    'Direction', 'Session', 'TargetValue', 'ClassThreshold', ...
    'TrueLabel', 'PredictedLabel', 'PredictedProbability', 'FoldSession'};
end
