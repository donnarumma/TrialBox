%% MAIN_SAPIENZA_BEHAVIOUR_NEURAL_RT_PRIMARY_MODEL_SUMMARY
% Primary-model RT summary for the current interpretation.
%
% Frontal is summarized with the trajectory-aware ModelM7, whereas Parietal
% is summarized with the direct-coupling ModelM2. The figure keeps the
% exploratory trajectory-aware Parietal results out of the main comparison.

clear;
close all;

analysisPms = struct();
analysisPms.inputRoot = '/TESTS/SAPIENZA/BEHAVIOUR/NeuralCorrelation';
analysisPms.outputRoot = fullfile(analysisPms.inputRoot, 'RTPattern');
analysisPms.visible = 'off';
analysisPms.directionIds = 1:8;
analysisPms.monkeySColor = [0.20 0.40 0.85];
analysisPms.monkeyKColor = [0.20 0.60 0.20];
analysisPms.balanceColor = [0.35 0.35 0.35];

plotDir = fullfile(analysisPms.outputRoot, 'plots');
pdfDir = fullfile(plotDir, 'PDF');
figDir = fullfile(plotDir, 'FIG');
if ~exist(pdfDir, 'dir')
    mkdir(pdfDir);
end
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

directionalPath = fullfile(analysisPms.inputRoot, ...
    'BehaviourNeural_directionalAllTargetCorrelations.csv');
globalPath = fullfile(analysisPms.inputRoot, ...
    'BehaviourNeural_correlations_sorted.csv');
if ~exist(directionalPath, 'file')
    error('%s:MissingInput', mfilename, 'Missing %s', directionalPath);
end
if ~exist(globalPath, 'file')
    error('%s:MissingInput', mfilename, 'Missing %s', globalPath);
end

directionalTable = readtable(directionalPath);
globalTable = readtable(globalPath);

summaryTable = local_build_summary_table(directionalTable, globalTable);
writetable(summaryTable, fullfile(analysisPms.outputRoot, ...
    'BehaviourNeural_RT_primaryModelSummary.csv'));

figureHandle = figure( ...
    'Color', 'w', ...
    'Visible', analysisPms.visible, ...
    'Name', 'Primary-model RT summary', ...
    'Position', [80 80 1250 980]);
tileLayoutHandle = tiledlayout(2, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

nexttile;
local_plot_single_target(summaryTable, analysisPms, ...
    'Frontal', 'M7', 'NeuralBalanceSToKMinusKToS', ...
    'Frontal M7', 'N_{bal}=I_K-I_S', analysisPms.balanceColor);

nexttile;
local_plot_incoming_pair(summaryTable, analysisPms, ...
    'Frontal', 'M7', 'Frontal M7 incoming weights');

nexttile;
local_plot_single_target(summaryTable, analysisPms, ...
    'Parietal', 'M2', 'NeuralBalanceSToKMinusKToS', ...
    'Parietal M2', 'N_{bal}=I_K-I_S', analysisPms.balanceColor);

nexttile;
local_plot_incoming_pair(summaryTable, analysisPms, ...
    'Parietal', 'M2', 'Parietal M2 incoming weights');

title(tileLayoutHandle, ...
    'S RT advantage and primary-model neural readouts', ...
    'Interpreter', 'tex');

figureStem = 'BehaviourNeural_RT_PrimaryModelSummary';
exportgraphics(figureHandle, fullfile(pdfDir, [figureStem '.pdf']), ...
    'ContentType', 'vector');
savefig(figureHandle, fullfile(figDir, [figureStem '.fig']));
close(figureHandle);

fprintf('Primary-model RT summary saved in %s\n', pdfDir);


function summaryTable = local_build_summary_table(directionalTable, globalTable)
focusPairs = cell2table({ ...
    'Frontal',  'M7', 'NeuralBalanceSToKMinusKToS', 'N_{bal}=I_K-I_S'; ...
    'Frontal',  'M7', 'JointWeightSToK',             'I_K=S\rightarrowK'; ...
    'Frontal',  'M7', 'JointWeightKToS',             'I_S=K\rightarrowS'; ...
    'Parietal', 'M2', 'NeuralBalanceSToKMinusKToS', 'N_{bal}=I_K-I_S'; ...
    'Parietal', 'M2', 'JointWeightSToK',             'I_K=S\rightarrowK'; ...
    'Parietal', 'M2', 'JointWeightKToS',             'I_S=K\rightarrowS'}, ...
    'VariableNames', {'Chamber', 'Model', 'NeuralReadout', 'ReadoutLabel'});

summaryRows = {};
for focusPairCounter = 1:height(focusPairs)
    focusPair = focusPairs(focusPairCounter, :);

    directionMask = strcmp(directionalTable.Chamber, focusPair.Chamber{1}) & ...
        strcmp(directionalTable.Model, focusPair.Model{1}) & ...
        strcmp(directionalTable.Predictor, 'SAdvantage_RT') & ...
        strcmp(directionalTable.Target, focusPair.NeuralReadout{1});
    pairDirectionalTable = sortrows(directionalTable(directionMask, :), 'Direction');

    globalMask = strcmp(globalTable.Chamber, focusPair.Chamber{1}) & ...
        strcmp(globalTable.Model, focusPair.Model{1}) & ...
        strcmp(globalTable.Predictor, 'SAdvantage_RT') & ...
        strcmp(globalTable.Target, focusPair.NeuralReadout{1});
    pairGlobalTable = globalTable(globalMask, :);

    if isempty(pairGlobalTable)
        globalRho = NaN;
        globalPValue = NaN;
        globalSampleCount = NaN;
    else
        globalRho = pairGlobalTable.Rho(1);
        globalPValue = pairGlobalTable.PValue(1);
        globalSampleCount = pairGlobalTable.N(1);
    end

    for directionCounter = 1:height(pairDirectionalTable)
        summaryRows(end + 1, :) = { ...
            focusPair.Chamber{1}, focusPair.Model{1}, ...
            focusPair.NeuralReadout{1}, focusPair.ReadoutLabel{1}, ...
            pairDirectionalTable.Direction(directionCounter), ...
            pairDirectionalTable.Rho(directionCounter), ...
            pairDirectionalTable.PValue(directionCounter), ...
            pairDirectionalTable.N(directionCounter), ...
            globalRho, globalPValue, globalSampleCount}; %#ok<AGROW>
    end
end

summaryTable = cell2table(summaryRows, ...
    'VariableNames', {'Chamber', 'Model', 'NeuralReadout', ...
    'ReadoutLabel', 'Direction', 'DirectionRho', 'DirectionPValue', ...
    'DirectionSampleCount', 'GlobalRho', 'GlobalPValue', ...
    'GlobalSampleCount'});
end


function local_plot_single_target(summaryTable, analysisPms, chamberName, ...
    modelName, neuralReadout, titleText, readoutLabel, barColor)
targetMask = strcmp(summaryTable.Chamber, chamberName) & ...
    strcmp(summaryTable.Model, modelName) & ...
    strcmp(summaryTable.NeuralReadout, neuralReadout);
targetTable = sortrows(summaryTable(targetMask, :), 'Direction');

bar(targetTable.Direction, targetTable.DirectionRho, ...
    'FaceColor', barColor, ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.82);
hold on;
yline(0, 'k-', 'LineWidth', 0.8);
yline(targetTable.GlobalRho(1), '--', ...
    'Color', [0.10 0.10 0.10], ...
    'LineWidth', 1.2);
local_format_axes(analysisPms);
title({titleText, sprintf('%s; global \\rho=%+.2f, p=%s', ...
    readoutLabel, targetTable.GlobalRho(1), ...
    local_pvalue_text(targetTable.GlobalPValue(1)))}, ...
    'Interpreter', 'tex');
end


function local_plot_incoming_pair(summaryTable, analysisPms, chamberName, ...
    modelName, titleText)
targetNames = {'JointWeightSToK', 'JointWeightKToS'};
rhoByDirection = nan(numel(analysisPms.directionIds), numel(targetNames));
globalRho = nan(1, numel(targetNames));
globalPValue = nan(1, numel(targetNames));

for targetCounter = 1:numel(targetNames)
    targetMask = strcmp(summaryTable.Chamber, chamberName) & ...
        strcmp(summaryTable.Model, modelName) & ...
        strcmp(summaryTable.NeuralReadout, targetNames{targetCounter});
    targetTable = sortrows(summaryTable(targetMask, :), 'Direction');
    rhoByDirection(:, targetCounter) = targetTable.DirectionRho;
    globalRho(targetCounter) = targetTable.GlobalRho(1);
    globalPValue(targetCounter) = targetTable.GlobalPValue(1);
end

barHandle = bar(analysisPms.directionIds, rhoByDirection, ...
    'grouped', ...
    'EdgeColor', 'none');
barHandle(1).FaceColor = analysisPms.monkeySColor;
barHandle(2).FaceColor = analysisPms.monkeyKColor;
barHandle(1).FaceAlpha = 0.82;
barHandle(2).FaceAlpha = 0.82;
hold on;
yline(0, 'k-', 'LineWidth', 0.8);
yline(globalRho(1), '--', 'Color', analysisPms.monkeySColor, 'LineWidth', 1.2);
yline(globalRho(2), '--', 'Color', analysisPms.monkeyKColor, 'LineWidth', 1.2);
local_format_axes(analysisPms);
legend({'I_K=S\rightarrowK', 'I_S=K\rightarrowS'}, ...
    'Location', 'best', ...
    'Interpreter', 'tex');
title({titleText, sprintf('global \\rho: I_K=%+.2f (p=%s), I_S=%+.2f (p=%s)', ...
    globalRho(1), local_pvalue_text(globalPValue(1)), ...
    globalRho(2), local_pvalue_text(globalPValue(2)))}, ...
    'Interpreter', 'tex');
end


function local_format_axes(analysisPms)
grid on;
xticks(analysisPms.directionIds);
xticklabels(compose('D%d', analysisPms.directionIds));
ylim([-0.75 0.75]);
ylabel('Spearman \rho', 'Interpreter', 'tex');
set(gca, 'FontSize', 10);
end


function pValueText = local_pvalue_text(pValue)
if ~isfinite(pValue)
    pValueText = 'NA';
elseif pValue < 1e-3
    pValueText = sprintf('%.1e', pValue);
else
    pValueText = sprintf('%.3f', pValue);
end
end
