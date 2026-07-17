%% MAIN_SAPIENZA_BEHAVIOUR_NEURAL_CLASSIFICATION
% Lightweight classification tests for behaviour-neural Sapienza readouts.
%
% The goal is not to build a predictive black box.  Instead, we ask whether
% simple DCM readouts contain enough information to classify behaviourally
% meaningful signs, using leave-one-session-out validation.

clear;
close all;

analysisPms = struct();
analysisPms.inputRoot = '/TESTS/SAPIENZA/BEHAVIOUR/NeuralCorrelation';
analysisPms.outputRoot = '/TESTS/SAPIENZA/BEHAVIOUR/NeuralClassification';
analysisPms.visible = 'off';
analysisPms.minimumSamplesPerClass = 6;
analysisPms.permutationIterations = 250;
analysisPms.ridgeLambda = 1;
analysisPms.maximumLogisticIterations = 80;
analysisPms.logisticTolerance = 1e-7;
analysisPms.randomSeed = 17;
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
analysisPms.taskRows = cell2table({ ...
    'S_AE_Benefit_Positive',     'SJointBenefit_AE_abs', 'SignPositive', 'S improves |AE| in joint action'; ...
    'K_AE_Benefit_Positive',     'KJointBenefit_AE_abs', 'SignPositive', 'K improves |AE| in joint action'; ...
    'S_AE_Advantage_Positive',   'SAdvantage_AE_abs',    'SignPositive', 'S is more accurate than K in joint action'; ...
    'S_RT_Advantage_Positive',   'SAdvantage_RT',        'SignPositive', 'S is faster than K in joint action'; ...
    'S_RT_Benefit_Positive',     'SJointBenefit_RT',     'SignPositive', 'S is faster in joint than when acting alone'; ...
    'K_RT_Benefit_Positive',     'KJointBenefit_RT',     'SignPositive', 'K is faster in joint than when acting alone'; ...
    'S_AE_Benefit_High',         'SJointBenefit_AE_abs', 'MedianHigh',   'S has higher-than-median |AE| joint benefit'; ...
    'K_AE_Benefit_High',         'KJointBenefit_AE_abs', 'MedianHigh',   'K has higher-than-median |AE| joint benefit'; ...
    'S_AE_Advantage_High',       'SAdvantage_AE_abs',    'MedianHigh',   'S has higher-than-median |AE| advantage'; ...
    'S_RT_Advantage_High',       'SAdvantage_RT',        'MedianHigh',   'S has higher-than-median RT advantage'; ...
    'S_RT_Benefit_High',         'SJointBenefit_RT',     'MedianHigh',   'S has higher-than-median RT joint benefit'; ...
    'K_RT_Benefit_High',         'KJointBenefit_RT',     'MedianHigh',   'K has higher-than-median RT joint benefit'}, ...
    'VariableNames', {'TaskId', 'TargetVariable', 'ClassRule', 'PositiveClassMeaning'});
analysisPms.featureRows = local_feature_rows();

rng(analysisPms.randomSeed, 'twister');

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

[resultTable, predictionTable] = local_run_classification_grid(jointTable, analysisPms);
topResultTable = local_build_top_result_table(resultTable);

writetable(resultTable, fullfile(analysisPms.outputRoot, ...
    'BehaviourNeural_classificationResults.csv'));
writetable(predictionTable, fullfile(analysisPms.outputRoot, ...
    'BehaviourNeural_classificationPredictions.csv'));
writetable(topResultTable, fullfile(analysisPms.outputRoot, ...
    'BehaviourNeural_classificationTopResults.csv'));

local_plot_top_results(topResultTable, analysisPms, pdfRoot, figRoot);
local_plot_primary_model_results(resultTable, analysisPms, pdfRoot, figRoot);

save(fullfile(analysisPms.outputRoot, 'BehaviourNeural_classification.mat'), ...
    'analysisPms', 'resultTable', 'predictionTable', 'topResultTable', '-v7.3');

fprintf('Behaviour-neural classification analysis saved in %s\n', analysisPms.outputRoot);


function featureRows = local_feature_rows()
featureRows = cell2table({ ...
    'DirectionOnly',           {'DirectionSin', 'DirectionCos'}, ...
    'Direction-only baseline'; ...
    'JointWeights',            {'JointWeightSToK', 'JointWeightKToS', 'NeuralBalanceSToKMinusKToS'}, ...
    'Joint-action crossed weights'; ...
    'JointGains',              {'JointGainSToK', 'JointGainKToS', 'NeuralJointGainBalance'}, ...
    'Joint-minus-active gain readouts'; ...
    'AllNeural',               {'JointWeightSToK', 'JointWeightKToS', 'NeuralBalanceSToKMinusKToS', ...
                                'JointGainSToK', 'JointGainKToS', 'NeuralJointGainBalance'}, ...
    'All DCM readouts'; ...
    'AllNeuralPlusDirection',  {'JointWeightSToK', 'JointWeightKToS', 'NeuralBalanceSToKMinusKToS', ...
                                'JointGainSToK', 'JointGainKToS', 'NeuralJointGainBalance', ...
                                'DirectionSin', 'DirectionCos'}, ...
    'All DCM readouts plus direction baseline'}, ...
    'VariableNames', {'FeatureSet', 'FeatureNames', 'Description'});
end


function jointTable = local_add_direction_features(jointTable)
directionAngles = 2 .* pi .* (jointTable.Direction - 1) ./ 8;
jointTable.DirectionSin = sin(directionAngles);
jointTable.DirectionCos = cos(directionAngles);
end


function [resultTable, predictionTable] = local_run_classification_grid(jointTable, analysisPms)
resultRows = {};
predictionRows = {};

for modelRowCounter = 1:height(analysisPms.modelRows)
    modelRow = analysisPms.modelRows(modelRowCounter, :);
    modelMask = strcmp(jointTable.Chamber, modelRow.Chamber{1}) & ...
        strcmp(jointTable.Model, modelRow.Model{1});
    modelTable = jointTable(modelMask, :);

    for taskCounter = 1:height(analysisPms.taskRows)
        taskRow = analysisPms.taskRows(taskCounter, :);

        for featureCounter = 1:height(analysisPms.featureRows)
            featureRow = analysisPms.featureRows(featureCounter, :);
            featureNames = featureRow.FeatureNames{1};
            [featureMatrix, labelVector, targetValues, classThreshold, sessionLabels, directionVector, ...
                validTable] = local_prepare_classification_data( ...
                modelTable, taskRow.TargetVariable{1}, taskRow.ClassRule{1}, featureNames);

            [resultRow, samplePredictions] = local_evaluate_classifier( ...
                featureMatrix, labelVector, targetValues, classThreshold, sessionLabels, directionVector, ...
                modelRow, taskRow, featureRow, validTable, analysisPms);

            resultRows(end + 1, :) = resultRow; %#ok<AGROW>
            predictionRows = [predictionRows; samplePredictions]; %#ok<AGROW>
        end
    end
end

resultTable = cell2table(resultRows, 'VariableNames', local_result_columns());
predictionTable = cell2table(predictionRows, 'VariableNames', local_prediction_columns());
end


function [featureMatrix, labelVector, targetValues, classThreshold, sessionLabels, directionVector, ...
    validTable] = local_prepare_classification_data(modelTable, targetName, classRule, featureNames)
targetValues = modelTable.(targetName);
featureMatrix = table2array(modelTable(:, featureNames));
validMask = isfinite(targetValues) & all(isfinite(featureMatrix), 2);

switch classRule
    case 'SignPositive'
        classThreshold = 0;
        validMask = validMask & targetValues ~= classThreshold;
    case 'MedianHigh'
        classThreshold = median(targetValues(validMask), 'omitnan');
        validMask = validMask & targetValues ~= classThreshold;
    otherwise
        error('%s:UnknownClassRule', mfilename, 'Unknown class rule: %s', classRule);
end

validTable = modelTable(validMask, :);
featureMatrix = featureMatrix(validMask, :);
targetValues = targetValues(validMask);
labelVector = targetValues > classThreshold;
sessionLabels = validTable.Session;
directionVector = validTable.Direction;
end


function [resultRow, predictionRows] = local_evaluate_classifier( ...
    featureMatrix, labelVector, targetValues, classThreshold, sessionLabels, directionVector, ...
    modelRow, taskRow, featureRow, validTable, analysisPms)

sampleCount = numel(labelVector);
positiveCount = sum(labelVector);
negativeCount = sampleCount - positiveCount;
sessionCount = numel(unique(sessionLabels, 'stable'));

isValidProblem = sampleCount > 0 && ...
    positiveCount >= analysisPms.minimumSamplesPerClass && ...
    negativeCount >= analysisPms.minimumSamplesPerClass && ...
    sessionCount >= 3;

if isValidProblem
    [predictedLabels, predictedProbability, foldLabels] = local_leave_one_session_out( ...
        featureMatrix, labelVector, sessionLabels, analysisPms);
    metricStruct = local_classification_metrics(labelVector, predictedLabels, predictedProbability);
    [permutationPValue, permutationMean, permutationUpper95] = local_permutation_test( ...
        featureMatrix, labelVector, sessionLabels, analysisPms);
else
    predictedLabels = false(sampleCount, 1);
    predictedProbability = nan(sampleCount, 1);
    foldLabels = strings(sampleCount, 1);
    metricStruct = local_empty_metrics();
    permutationPValue = NaN;
    permutationMean = NaN;
    permutationUpper95 = NaN;
end

featureNames = featureRow.FeatureNames{1};
resultRow = { ...
    modelRow.Chamber{1}, modelRow.Model{1}, modelRow.Description{1}, ...
    taskRow.TaskId{1}, taskRow.TargetVariable{1}, taskRow.PositiveClassMeaning{1}, ...
    taskRow.ClassRule{1}, classThreshold, ...
    featureRow.FeatureSet{1}, strjoin(featureNames, '+'), featureRow.Description{1}, ...
    sampleCount, sessionCount, positiveCount, negativeCount, ...
    metricStruct.Accuracy, metricStruct.BalancedAccuracy, metricStruct.Sensitivity, ...
    metricStruct.Specificity, metricStruct.AUC, permutationPValue, ...
    permutationMean, permutationUpper95, isValidProblem};

predictionRows = cell(sampleCount, numel(local_prediction_columns()));
for sampleCounter = 1:sampleCount
    predictionRows(sampleCounter, :) = { ...
        modelRow.Chamber{1}, modelRow.Model{1}, taskRow.TaskId{1}, ...
        taskRow.TargetVariable{1}, taskRow.ClassRule{1}, classThreshold, featureRow.FeatureSet{1}, ...
        validTable.Session{sampleCounter}, directionVector(sampleCounter), ...
        targetValues(sampleCounter), labelVector(sampleCounter), ...
        predictedLabels(sampleCounter), predictedProbability(sampleCounter), ...
        string(foldLabels(sampleCounter))};
end
end


function [predictedLabels, predictedProbability, foldLabels] = local_leave_one_session_out( ...
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


function probability = local_predict_ridge_logistic(trainFeatures, trainLabels, ...
    testFeatures, analysisPms)
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


function [trainFeatures, testFeatures] = local_standardize_from_train(trainFeatures, testFeatures)
featureMean = mean(trainFeatures, 1, 'omitnan');
featureScale = std(trainFeatures, 0, 1, 'omitnan');
featureScale(~isfinite(featureScale) | featureScale == 0) = 1;
trainFeatures = (trainFeatures - featureMean) ./ featureScale;
testFeatures = (testFeatures - featureMean) ./ featureScale;
end


function probability = local_sigmoid(linearPredictor)
linearPredictor = max(min(linearPredictor, 35), -35);
probability = 1 ./ (1 + exp(-linearPredictor));
end


function metricStruct = local_classification_metrics(trueLabels, predictedLabels, predictedProbability)
trueLabels = logical(trueLabels);
predictedLabels = logical(predictedLabels);
truePositive = sum(predictedLabels & trueLabels);
trueNegative = sum(~predictedLabels & ~trueLabels);
falsePositive = sum(predictedLabels & ~trueLabels);
falseNegative = sum(~predictedLabels & trueLabels);

metricStruct = struct();
metricStruct.Accuracy = mean(predictedLabels == trueLabels);
metricStruct.Sensitivity = truePositive ./ max(truePositive + falseNegative, 1);
metricStruct.Specificity = trueNegative ./ max(trueNegative + falsePositive, 1);
metricStruct.BalancedAccuracy = mean([metricStruct.Sensitivity, metricStruct.Specificity]);
metricStruct.AUC = local_auc(trueLabels, predictedProbability);
end


function metricStruct = local_empty_metrics()
metricStruct = struct();
metricStruct.Accuracy = NaN;
metricStruct.Sensitivity = NaN;
metricStruct.Specificity = NaN;
metricStruct.BalancedAccuracy = NaN;
metricStruct.AUC = NaN;
end


function aucValue = local_auc(trueLabels, predictedProbability)
validMask = isfinite(predictedProbability);
trueLabels = trueLabels(validMask);
predictedProbability = predictedProbability(validMask);
positiveCount = sum(trueLabels);
negativeCount = sum(~trueLabels);
if positiveCount == 0 || negativeCount == 0
    aucValue = NaN;
    return;
end

rankValues = local_tied_ranks(predictedProbability);
positiveRankSum = sum(rankValues(trueLabels));
aucValue = (positiveRankSum - positiveCount * (positiveCount + 1) / 2) ./ ...
    (positiveCount * negativeCount);
end


function rankValues = local_tied_ranks(values)
[sortedValues, sortOrder] = sort(values(:));
rankValues = nan(size(values(:)));
startIndex = 1;
while startIndex <= numel(sortedValues)
    endIndex = startIndex;
    while endIndex < numel(sortedValues) && sortedValues(endIndex + 1) == sortedValues(startIndex)
        endIndex = endIndex + 1;
    end
    averageRank = mean(startIndex:endIndex);
    rankValues(sortOrder(startIndex:endIndex)) = averageRank;
    startIndex = endIndex + 1;
end
end


function [permutationPValue, permutationMean, permutationUpper95] = local_permutation_test( ...
    featureMatrix, labelVector, sessionLabels, analysisPms)
if analysisPms.permutationIterations <= 0
    permutationPValue = NaN;
    permutationMean = NaN;
    permutationUpper95 = NaN;
    return;
end

[observedPredictedLabels, observedPredictedProbability] = local_leave_one_session_out( ...
    featureMatrix, labelVector, sessionLabels, analysisPms);
observedMetrics = local_classification_metrics( ...
    labelVector, observedPredictedLabels, observedPredictedProbability);
nullBalancedAccuracy = nan(analysisPms.permutationIterations, 1);

for permutationCounter = 1:analysisPms.permutationIterations
    shuffledLabels = labelVector(randperm(numel(labelVector)));
    [permutedLabels, permutedProbability] = local_leave_one_session_out( ...
        featureMatrix, shuffledLabels, sessionLabels, analysisPms);
    permutedMetrics = local_classification_metrics( ...
        shuffledLabels, permutedLabels, permutedProbability);
    nullBalancedAccuracy(permutationCounter) = permutedMetrics.BalancedAccuracy;
end

permutationMean = mean(nullBalancedAccuracy, 'omitnan');
permutationUpper95 = prctile(nullBalancedAccuracy, 95);
permutationPValue = (1 + sum(nullBalancedAccuracy >= observedMetrics.BalancedAccuracy)) ./ ...
    (1 + sum(isfinite(nullBalancedAccuracy)));
end


function topResultTable = local_build_top_result_table(resultTable)
topResultTable = resultTable(resultTable.IsValidProblem, :);
if isempty(topResultTable)
    return;
end

sortPermutation = topResultTable.PermutationPValue;
sortPermutation(~isfinite(sortPermutation)) = Inf;
topResultTable.SortPermutation = sortPermutation;
topResultTable.SortBalancedAccuracy = -topResultTable.BalancedAccuracy;
topResultTable = sortrows(topResultTable, ...
    {'SortPermutation', 'SortBalancedAccuracy', 'Chamber', 'Model', 'TaskId'});
topResultTable.SortPermutation = [];
topResultTable.SortBalancedAccuracy = [];
end


function local_plot_top_results(topResultTable, analysisPms, pdfRoot, figRoot)
if isempty(topResultTable)
    return;
end

topCount = min(24, height(topResultTable));
plotTable = topResultTable(1:topCount, :);
labelText = strings(topCount, 1);
for rowCounter = 1:topCount
    labelText(rowCounter) = sprintf('%s %s | %s | %s', ...
        plotTable.Chamber{rowCounter}, plotTable.Model{rowCounter}, ...
        strrep(plotTable.TaskId{rowCounter}, '_', '\_'), ...
        plotTable.FeatureSet{rowCounter});
end

figureHandle = figure( ...
    'Color', 'w', ...
    'Visible', analysisPms.visible, ...
    'Name', 'Top behaviour-neural classification results', ...
    'Position', [80 80 1400 1100]);

barh(1:topCount, flipud(plotTable.BalancedAccuracy), ...
    'FaceColor', [0.25 0.43 0.72], ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.85);
hold on;
xline(0.5, 'k--', 'LineWidth', 1.1);
xlabel('Leave-one-session-out balanced accuracy');
yticks(1:topCount);
yticklabels(flipud(labelText));
xlim([0 1]);
grid on;
title(sprintf('Top classification tests; permutation N=%d', ...
    analysisPms.permutationIterations), 'Interpreter', 'none');
set(gca, 'FontSize', 9, 'TickLabelInterpreter', 'tex');

figureStem = 'BehaviourNeural_Classification_TopResults';
exportgraphics(figureHandle, fullfile(pdfRoot, [figureStem '.pdf']), ...
    'ContentType', 'vector');
savefig(figureHandle, fullfile(figRoot, [figureStem '.fig']));
close(figureHandle);
end


function local_plot_primary_model_results(resultTable, analysisPms, pdfRoot, figRoot)
focusRows = cell2table({ ...
    'Frontal',  'M7', 'S_AE_Benefit_High',       'Frontal M7: high S |AE| benefit'; ...
    'Frontal',  'M7', 'S_RT_Advantage_High',     'Frontal M7: high S RT advantage'; ...
    'Parietal', 'M2', 'S_AE_Benefit_High',       'Parietal M2: high S |AE| benefit'; ...
    'Parietal', 'M2', 'S_RT_Advantage_High',     'Parietal M2: high S RT advantage'}, ...
    'VariableNames', {'Chamber', 'Model', 'TaskId', 'TitleText'});

figureHandle = figure( ...
    'Color', 'w', ...
    'Visible', analysisPms.visible, ...
    'Name', 'Primary-model behaviour-neural classification', ...
    'Position', [80 80 1450 900]);
tileLayoutHandle = tiledlayout(2, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

for focusCounter = 1:height(focusRows)
    focusRow = focusRows(focusCounter, :);
    nexttile;
    focusMask = strcmp(resultTable.Chamber, focusRow.Chamber{1}) & ...
        strcmp(resultTable.Model, focusRow.Model{1}) & ...
        strcmp(resultTable.TaskId, focusRow.TaskId{1}) & ...
        resultTable.IsValidProblem;
    focusTable = resultTable(focusMask, :);
    [~, featureOrder] = ismember(analysisPms.featureRows.FeatureSet, focusTable.FeatureSet);
    featureOrder = featureOrder(featureOrder > 0);
    focusTable = focusTable(featureOrder, :);

    bar(focusTable.BalancedAccuracy, ...
        'FaceColor', [0.25 0.43 0.72], ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.85);
    hold on;
    yline(0.5, 'k--', 'LineWidth', 1.0);
    xticks(1:height(focusTable));
    xticklabels(focusTable.FeatureSet);
    xtickangle(35);
    ylim([0 1]);
    ylabel('Balanced accuracy');
    grid on;
    title(focusRow.TitleText{1}, 'Interpreter', 'none');
end

title(tileLayoutHandle, ...
    'Primary-model classification: DCM readouts versus direction-only baseline', ...
    'Interpreter', 'none');
figureStem = 'BehaviourNeural_Classification_PrimaryModels';
exportgraphics(figureHandle, fullfile(pdfRoot, [figureStem '.pdf']), ...
    'ContentType', 'vector');
savefig(figureHandle, fullfile(figRoot, [figureStem '.fig']));
close(figureHandle);
end


function variableNames = local_result_columns()
variableNames = {'Chamber', 'Model', 'ModelDescription', 'TaskId', ...
    'TargetVariable', 'PositiveClassMeaning', 'ClassRule', 'ClassThreshold', ...
    'FeatureSet', 'FeatureNames', ...
    'FeatureDescription', 'N', 'SessionCount', 'PositiveCount', ...
    'NegativeCount', 'Accuracy', 'BalancedAccuracy', 'Sensitivity', ...
    'Specificity', 'AUC', 'PermutationPValue', 'PermutationMeanBalancedAccuracy', ...
    'PermutationUpper95BalancedAccuracy', 'IsValidProblem'};
end


function variableNames = local_prediction_columns()
variableNames = {'Chamber', 'Model', 'TaskId', 'TargetVariable', ...
    'ClassRule', 'ClassThreshold', 'FeatureSet', 'Session', 'Direction', 'TargetValue', 'TrueLabel', ...
    'PredictedLabel', 'PredictedProbability', 'FoldSession'};
end
